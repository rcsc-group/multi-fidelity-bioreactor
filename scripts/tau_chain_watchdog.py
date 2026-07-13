"""Status/self-heal watchdog for the L9+L10 tau-only sweep (theta=7deg, 9 RPM conditions).

On each invocation:
  1. Checks each of 9 L9 (single-shot) and 9 L10 (checkpoint-restart chain) conditions
     for a finite tau_100_max in results.json.
  2. For L10 chains: if a segment finished but the next segment in its chain has
     neither a result nor a queued/running job, resubmits it from the checkpoint
     (self-healing against timeouts/preemption).
  3. Exits 0 when all 18 conditions are done, 1 otherwise.

Safe to run repeatedly -- never double-submits a queued job.
"""
from __future__ import annotations
import json, math, subprocess, sys
from pathlib import Path

ROOT     = Path(__file__).parents[1]
RUNS     = ROOT / "runs"
SCRATCH  = Path("/oscar/scratch/eaguerov/mpi_runs")
TEMPLATE = ROOT / "config" / "slurm_mpi_template.sh"
EXP_L9   = ROOT / "experiments" / "sweep_tau_theta7_l9"
EXP_L10  = ROOT / "experiments" / "sweep_tau_theta7_l10"

# Previously hardcoded to only 4 of the 9 actual L9 run_ids (17.5, 20.0, 22.5,
# 25.0), silently excluding 27.5-37.5 from tracking and undercounting the
# true total (13 instead of 18). Load the full list from the sweep metadata.
L9_RUN_IDS = json.loads((EXP_L9 / "_sweep_metadata.json").read_text())["run_ids"]


def _has_result(run_id: str) -> bool:
    # slurm_mpi_template.sh writes results.json to the CANONICAL runs/ dir
    # when a run's params.json has _canonical_run_dir/_experiment_dir set
    # (the normal case for tracked sweep segments); manually-rescued runs
    # (postprocessed directly against the scratch dir, bypassing the
    # template) only ever have it in SCRATCH. Check both.
    for f in (RUNS / run_id / "results.json", SCRATCH / run_id / "results.json"):
        if not f.exists():
            continue
        try:
            r = json.loads(f.read_text())
            v = r.get("tau_100_max")
            if v is not None and math.isfinite(float(v)):
                return True
        except Exception:
            pass
    return False


def _active_job_ids() -> set[str]:
    try:
        out = subprocess.run(
            ["squeue", "-u", "eaguerov", "--format=%.10i %.8T", "--noheader"],
            capture_output=True, text=True, timeout=15,
        ).stdout
    except Exception:
        return set()
    active = set()
    for line in out.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[1] in ("RUNNING", "PENDING"):
            active.add(parts[0].strip())
    return active


def _queued_run_ids(active_jids: set[str]) -> set[str]:
    ids = set()
    for jid_file in SCRATCH.glob("*/.slurm_jid"):
        rid = jid_file.parent.name
        jid = jid_file.read_text().strip()
        if jid in active_jids:
            ids.add(rid)
    return ids


def _resubmit(run_id: str, from_run_id: str) -> str | None:
    src_canon = RUNS / from_run_id
    next_canon = RUNS / run_id
    next_scratch = SCRATCH / run_id
    next_scratch.mkdir(parents=True, exist_ok=True)

    ckpt_src = src_canon / "checkpoint.dump"
    if not ckpt_src.exists():
        ckpt_src = SCRATCH / from_run_id / "checkpoint.dump"
    if not ckpt_src.exists():
        print(f"  [WARN] no checkpoint for {from_run_id}, cannot resubmit {run_id}")
        return None

    (next_scratch / "checkpoint.dump").write_bytes(ckpt_src.read_bytes())
    params = json.loads((next_canon / "params.json").read_text())
    (next_scratch / "params.json").write_text(json.dumps(params, indent=2))

    walltime = params.get("_walltime", "2-00:00:00")
    mem      = params.get("_mem", "4G")
    ntasks   = params.get("_ntasks", 16)
    mail_user = params.get("_mail_user", "")
    mail_type = params.get("_mail_type", "FAIL")

    cmd = [
        "sbatch", "--no-requeue",
        f"--time={walltime}", f"--mem-per-cpu={mem}", f"--ntasks={ntasks}",
        "--cpus-per-task=1",
    ]
    if mail_user:
        cmd += [f"--mail-type={mail_type}", f"--mail-user={mail_user}"]
    cmd += [
        f"--export=NONE,PARAMS={next_scratch}/params.json,DUMP={next_scratch}/checkpoint.dump",
        str(TEMPLATE),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] sbatch failed for {run_id}: {result.stderr.strip()}")
        return None
    jid = result.stdout.strip().split()[-1]
    (next_scratch / ".slurm_jid").write_text(jid)
    return jid


def main() -> None:
    manifest = json.loads((EXP_L10 / "_chain_manifest.json").read_text())
    active = _active_job_ids()
    inflight = _queued_run_ids(active)

    total = len(L9_RUN_IDS) + len(manifest)
    n_done = 0

    print("=== L9 (single-shot) ===")
    for rid in L9_RUN_IDS:
        done = _has_result(rid)
        n_done += done
        status = "DONE" if done else ("queued/running" if rid in inflight else "MISSING (no result, not queued)")
        print(f"  {rid}: {status}")

    print("=== L10 (checkpoint chain) ===")
    submitted = []
    for rpm, chain in manifest.items():
        final = chain[-1]
        done = _has_result(final)
        n_done += done
        if done:
            print(f"  {rpm} rpm: DONE (final={final})")
            continue
        # find current position: last segment with a result, or chain[0] if none
        cur_idx = -1
        for i, rid in enumerate(chain):
            if _has_result(rid):
                cur_idx = i
        if cur_idx + 1 < len(chain):
            nxt = chain[cur_idx + 1]
            cur = chain[cur_idx] if cur_idx >= 0 else None
            if nxt in inflight:
                print(f"  {rpm} rpm: seg {cur_idx+1}/{len(chain)} running (run={nxt})")
            elif cur_idx == -1:
                # segment 0 hasn't produced a result and isn't queued -> report only
                # (it may still be running without a .slurm_jid marker, or genuinely dead)
                print(f"  {rpm} rpm: seg 0/{len(chain)} not queued, no result (run={nxt}) -- check manually")
            else:
                jid = _resubmit(nxt, cur)
                if jid:
                    submitted.append((rpm, nxt, jid))
                    print(f"  {rpm} rpm: seg {cur_idx} done, resubmitted seg {cur_idx+1} -> job {jid}")
                else:
                    print(f"  {rpm} rpm: seg {cur_idx} done, FAILED to resubmit seg {cur_idx+1}")

    print(f"\n{n_done}/{total} conditions complete.")
    if submitted:
        print(f"Self-healed {len(submitted)} broken chain segment(s).")

    sys.exit(0 if n_done == total else 1)


if __name__ == "__main__":
    main()
