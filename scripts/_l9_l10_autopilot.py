"""Autonomous coordinator for the L9 warm-start rerun + L10 pause policy.

Runs unattended (the user checks in every few days, no terminal access in
the meantime). On each invocation:

  1. Finds any currently-RUNNING bioreactor-mpi job whose run_id belongs to
     the L10 chain manifest. If such a job disappeared from the queue since
     the last invocation (i.e. it finished), waits for its auto-chained next
     segment to appear in its own .out log, then cancels that next segment
     while it's still PENDING -- this is the "pause L10, give L9 priority"
     policy. Never touches a job that's already RUNNING.

  2. Tracks the L9 pilot warm-start job (17.5 rpm, run_id 44133566). Once it
     finishes, sanity-checks its results.json (finite, physically reasonable
     tau values) via postprocess.py. If it looks good, submits the other 7
     prepared L9 warm-start jobs (configs already written to
     /oscar/data/dharri15/eaguerov/mpi_runs/<run_id>/params.json). If the
     pilot looks wrong, STOPS and does not submit the rest -- that needs a
     human to look at it.

  3. Tracks all 9 L9 conditions (the "warm start sweep") for completion and
     reports final tau_100_max / tau_mean_max vs Kim et al. once all are done.

State is persisted to a small JSON file so repeated invocations (e.g. every
few minutes from a Monitor loop) don't repeat work or lose track of which
L10 jobs were already seen running.
"""
from __future__ import annotations
import json, math, re, subprocess, sys, time
from pathlib import Path

ROOT       = Path(__file__).parents[1]
RUNS       = ROOT / "runs"
SCRATCH    = Path("/oscar/scratch/eaguerov/mpi_runs")
PERSISTENT = Path("/oscar/data/dharri15/eaguerov/mpi_runs")
TEMPLATE   = ROOT / "config" / "slurm_mpi_template.sh"
STATE_FILE = ROOT / "experiments" / "_autopilot_state.json"

L9_META  = json.loads((ROOT / "experiments/sweep_tau_theta7_l9/_sweep_metadata.json").read_text())
L9_RUN_IDS = set(L9_META["run_ids"])
L10_MANIFEST = json.loads((ROOT / "experiments/sweep_tau_theta7_l10/_chain_manifest.json").read_text())
L10_RUN_IDS = {rid for chain in L10_MANIFEST.values() for rid in chain}

PILOT_RUN_ID = "44133566"  # 17.5 rpm, biggest omega jump from the 22.5 rpm seed
OTHER_L9_RUN_IDS = [rid for rid in L9_META["run_ids"] if rid != "8994c04a"]  # 22.5 already done

KIM_TAU = {  # RPM -> (tau_liq_max, tau_liq_mean), from docs/kimetal2024
    17.5: (0.0955, 0.00077), 20.0: (0.0914, 0.00086), 22.5: (0.2299, 0.00190),
    25.0: (0.1256, 0.00125), 27.5: (0.1467, 0.00139), 30.0: (0.1735, 0.00161),
    32.5: (0.2060, 0.00180), 35.0: (0.2753, 0.00221), 37.5: (1.1419, 0.00314),
}


def _load_state() -> dict:
    if STATE_FILE.exists():
        return json.loads(STATE_FILE.read_text())
    return {"seen_running_l10_jobs": {}, "l9_batch_submitted": False, "pilot_verdict": None}


def _save_state(s: dict) -> None:
    STATE_FILE.write_text(json.dumps(s, indent=2))


def _squeue() -> list[dict]:
    out = subprocess.run(
        ["squeue", "-u", "eaguerov", "--format=%i|%T", "--noheader"],
        capture_output=True, text=True, timeout=30,
    ).stdout
    rows = []
    for line in out.splitlines():
        parts = line.split("|")
        if len(parts) == 2:
            rows.append({"jobid": parts[0].strip(), "state": parts[1].strip()})
    return rows


def _job_run_id(jobid: str) -> str | None:
    out = subprocess.run(
        ["sacct", "-j", jobid, "--format=SubmitLine%300", "--noheader", "-X"],
        capture_output=True, text=True, timeout=30,
    ).stdout
    for tok in out.split():
        if tok.startswith("PARAMS=") or "PARAMS=" in tok:
            path = tok.split("PARAMS=")[-1].split(",")[0]
            return Path(path).parent.name
    return None


def _has_result(run_id: str) -> bool:
    # L9 was deliberately moved to PERSISTENT-only storage for the warm-start
    # rerun (see conversation history) -- SCRATCH still holds STALE pre-fix
    # coarse-sampling data for every original L9 run_id and must never be
    # checked here, or a finished-looking old file gets mistaken for a fresh
    # result (this exact bug fired once: it read 44133566's old scratch data
    # and declared the pilot done while the real pilot job was still PENDING).
    for base in (RUNS, PERSISTENT):
        f = base / run_id / "results.json"
        if f.exists():
            try:
                r = json.loads(f.read_text())
                v = r.get("tau_100_max")
                if v is not None and math.isfinite(float(v)):
                    return True
            except Exception:
                pass
    return False


def _postprocess(run_id: str) -> dict | None:
    run_dir = PERSISTENT / run_id
    if (run_dir / "shear_stress.dat").exists():
        sys.path.insert(0, str(ROOT))
        from scripts import postprocess as pp
        return pp.main(str(run_dir))
    return None


def _submit(run_id: str) -> str | None:
    params = PERSISTENT / run_id / "params.json"
    dump = PERSISTENT / run_id / "checkpoint.dump"
    cmd = [
        "sbatch", "--parsable", "--time=12:00:00", "--ntasks=16", "--cpus-per-task=1",
        "--mem-per-cpu=4G", "--mail-type=FAIL",
        "--mail-user=elvis_alexander_aguero_vera@brown.edu",
        f"--export=NONE,PARAMS={params},DUMP={dump}",
        str(TEMPLATE),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] sbatch failed for {run_id}: {result.stderr.strip()}")
        return None
    return result.stdout.strip()


def _pause_l10(state: dict) -> None:
    rows = _squeue()
    running_now = {}
    for row in rows:
        if row["state"] != "RUNNING":
            continue
        rid = _job_run_id(row["jobid"])
        if rid in L10_RUN_IDS:
            running_now[row["jobid"]] = rid

    prev = state["seen_running_l10_jobs"]
    # any job that WAS running last check but isn't anymore -> finished
    finished = [jid for jid in prev if jid not in running_now]
    for jid in finished:
        rid = prev[jid]
        out_log = ROOT / "logs" / f"slurm_{jid}.out"
        newjob = None
        for _ in range(40):  # poll up to ~10 min for postprocess+resubmit to land
            if out_log.exists():
                # slurm_mpi_template.sh prints "Submitted next segment: <rid> (job <jid>)"
                # -- NOT sbatch's own raw "Submitted batch job <jid>" stdout, which only
                # appears if you call sbatch directly yourself (as this script's own
                # _submit() does for L9). Match both so this doesn't silently miss the
                # L10 auto-chain line again (it already did once, letting an
                # auto-chained segment run uncancelled for a few minutes before being
                # caught manually).
                m = re.search(
                    r"Submitted (?:batch job (\d+)|next segment: \S+ \(job (\d+)\))",
                    out_log.read_text(),
                )
                if m:
                    newjob = m.group(1) or m.group(2)
            if newjob:
                break
            time.sleep(15)
        if newjob:
            st = subprocess.run(
                ["squeue", "-j", newjob, "--format=%T", "--noheader"],
                capture_output=True, text=True,
            ).stdout.strip()
            if st == "PENDING":
                subprocess.run(["scancel", newjob])
                print(f"PAUSED_L10: {rid} (job {jid}) finished; cancelled auto-chained "
                      f"segment (job {newjob}) to keep capacity free for L9")
            else:
                print(f"NOTE: {rid} (job {jid}) finished; next segment {newjob} "
                      f"already state={st}, left alone")
        else:
            print(f"NOTE: {rid} (job {jid}) finished; no auto-chain job found "
                  f"after 10 min (may be the chain's final segment)")

    state["seen_running_l10_jobs"] = running_now


def _drive_l9(state: dict) -> bool:
    """Returns True once all 9 L9 conditions are done."""
    if not state["l9_batch_submitted"]:
        if _has_result(PILOT_RUN_ID):
            res = _postprocess(PILOT_RUN_ID)
            tau_max = res.get("tau_100_max") if res else None
            ok = (
                tau_max is not None and math.isfinite(tau_max)
                and 0.02 < tau_max < 1.0  # sanity band around Kim's 0.0955 at 17.5 rpm
            )
            state["pilot_verdict"] = {"tau_100_max": tau_max, "ok": ok}
            if ok:
                print(f"PILOT_CONFIRMED: 17.5 rpm tau_100_max={tau_max:.4f} "
                      f"(Kim={KIM_TAU[17.5][0]}, ratio={tau_max/KIM_TAU[17.5][0]:.2f}) "
                      f"-- submitting remaining 7 L9 warm-start jobs")
                for rid in OTHER_L9_RUN_IDS:
                    if rid == PILOT_RUN_ID:
                        continue
                    jid = _submit(rid)
                    print(f"  submitted {rid} -> job {jid}")
                state["l9_batch_submitted"] = True
            else:
                print(f"PILOT_LOOKS_WRONG: tau_100_max={tau_max} -- NOT submitting "
                      f"the rest. Needs a human to look at this.")
                _save_state(state)
                return False
        else:
            return False  # pilot not done yet

    done = [rid for rid in L9_META["run_ids"] if _has_result(rid)]
    if len(done) == len(L9_META["run_ids"]):
        print("ALL_L9_DONE -- final comparison vs Kim et al.:")
        for rpm, rid in zip(L9_META["rpms"], L9_META["run_ids"]):
            r = _postprocess(rid) or {}
            tm, tmm = r.get("tau_100_max"), r.get("tau_mean_max")
            k_max, k_mean = KIM_TAU[rpm]
            print(f"  {rpm:>5} rpm  tau_100_max={tm:.4f} (ratio {tm/k_max:.2f})  "
                  f"tau_mean_max={tmm:.5f} (ratio {tmm/k_mean:.2f})")
        return True
    return False


def main() -> None:
    state = _load_state()
    _pause_l10(state)
    all_l9_done = _drive_l9(state)
    _save_state(state)
    sys.exit(0 if all_l9_done else 1)


if __name__ == "__main__":
    main()
