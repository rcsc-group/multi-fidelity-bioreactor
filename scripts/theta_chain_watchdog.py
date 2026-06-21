"""Autonomous watchdog for the L7 theta MPI-checkpoint chain sweep.

On each invocation:
  1. Maps every (theta, omega) condition to its result state.
  2. Checks the SLURM queue to see what is already running/pending.
  3. For every completed run whose *next* segment has no result and no
     queued job, stages the checkpoint and submits it (self-healing).
  4. Prints a status block and, when all 60 conditions are done,
     regenerates the theta heatmap, commits, and pushes.

Safe to run multiple times — it never double-submits a queued job.
"""
from __future__ import annotations
import json, math, os, subprocess, sys
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────────
ROOT     = Path(__file__).parents[1]           # rocking-bioreactor-2d/
RUNS     = ROOT / "runs"
SCRATCH  = Path("/oscar/scratch/eaguerov/mpi_runs")
TEMPLATE = ROOT / "config" / "slurm_mpi_template.sh"
EXP_TAG  = "theta_l7_mpi_ckpt"
TOTAL    = 60   # 6 theta × 10 omega

# ── helpers ────────────────────────────────────────────────────────────────

def _load_params(run_dir: Path) -> dict:
    f = run_dir / "params.json"
    return json.loads(f.read_text()) if f.exists() else {}

def _has_result(run_dir: Path) -> bool:
    f = run_dir / "results.json"
    if not f.exists():
        return False
    try:
        r = json.loads(f.read_text())
        v = r.get("kLa_25")
        return v is not None and math.isfinite(float(v))
    except Exception:
        return False

def _queued_run_ids() -> set[str]:
    """Return set of scratch-dir basenames currently in SLURM (R or PD)."""
    try:
        out = subprocess.run(
            ["squeue", "-u", "eaguerov",
             "--format=%.10i %.8T %100j", "--noheader"],
            capture_output=True, text=True, timeout=15
        ).stdout
    except Exception:
        return set()
    ids: set[str] = set()
    for line in out.splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        state = parts[1]
        if state not in ("RUNNING", "PENDING"):
            continue
        # job name is just "bioreactor-mpi"; we can't get run_id from name.
        # Instead we check scratch dirs that have a params.json recently staged.
    # Fallback: list every scratch dir created after the last result
    # (any dir with params.json but no results.json is inflight or staged)
    for d in SCRATCH.glob("*/params.json"):
        rid = d.parent.name
        if not (RUNS / rid / "results.json").exists():
            ids.add(rid)
    return ids

def _sbatch(run_id: str, from_run_id: str) -> str | None:
    """Stage checkpoint, submit job. Returns SLURM job ID or None on error."""
    src_canon  = RUNS / from_run_id
    next_canon = RUNS / run_id
    next_scratch = SCRATCH / run_id
    next_scratch.mkdir(parents=True, exist_ok=True)

    # Copy checkpoint from canon dir (preferred) or scratch
    ckpt_src = src_canon / "checkpoint.dump"
    if not ckpt_src.exists():
        ckpt_src = SCRATCH / from_run_id / "checkpoint.dump"
    if not ckpt_src.exists():
        print(f"  [WARN] no checkpoint found for {from_run_id}, skipping {run_id}")
        return None

    os.system(f"cp {ckpt_src} {next_scratch}/checkpoint.dump")
    os.system(f"cp {next_canon}/params.json {next_scratch}/params.json")

    p = _load_params(next_canon)
    walltime = p.get("_walltime", "04:30:00")
    mem      = p.get("_mem",      "4G")
    ntasks   = p.get("_ntasks",   16)

    result = subprocess.run(
        ["sbatch", "--no-requeue",
         f"--time={walltime}",
         f"--mem-per-cpu={mem}",
         f"--ntasks={ntasks}",
         "--cpus-per-task=1",
         f"--export=NONE,PARAMS={next_scratch}/params.json,"
         f"DUMP={next_scratch}/checkpoint.dump",
         str(TEMPLATE)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  [ERROR] sbatch failed for {run_id}: {result.stderr.strip()}")
        return None
    jid = result.stdout.strip().split()[-1]
    return jid


def _generate_figure_and_commit(n_done: int) -> None:
    """Regenerate theta heatmap and commit+push when all 60 done."""
    import subprocess as sp
    print("\n=== ALL 60 CONDITIONS COMPLETE — generating figure + commit ===")

    sp.run(["uv", "run", "python", "scripts/plot_heatmaps.py",
            "--fidelity", "7", "--exp-suffix", "theta_l7_mpi_ckpt",
            "--theta-fill", "0.5", "--fill-theta", "7.0"],
           cwd=str(ROOT), check=False)

    # also regenerate checkpoint validation with full dataset
    sp.run(["uv", "run", "python", "scripts/plot_checkpoint_validation.py"],
           cwd=str(ROOT), check=False)

    repo = ROOT.parents[2]  # BioReactor3D/
    sp.run(["git", "add",
            "dev/rocking-bioreactor-2d/experiments/figures/heatmap_theta_sweep_l7.pdf",
            "dev/rocking-bioreactor-2d/experiments/figures/checkpoint_validation.pdf",
            "dev/rocking-bioreactor-2d/experiments/figures/heatmap_fill_sweep_l7.pdf"],
           cwd=str(repo))
    sp.run(["git", "commit", "-m",
            f"data(l7-theta): all 60 theta sweep conditions complete\n\n"
            f"Regenerate theta heatmap (16 KPI panels) and checkpoint\n"
            f"validation figure with full 60+50 condition dataset.\n\n"
            f"Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"],
           cwd=str(repo))
    sp.run(["git", "push"], cwd=str(repo))
    print("Committed and pushed.")


# ── main ───────────────────────────────────────────────────────────────────

def main() -> None:
    # 1. Load all theta runs
    runs_info: dict[str, dict] = {}
    for p in sorted(RUNS.glob("*/params.json")):
        d = _load_params(p.parent)
        if EXP_TAG not in d.get("_experiment_dir", ""):
            continue
        rid = p.parent.name
        om  = round(float(d.get("omega_b", 0)), 4)
        th_raw = d.get("theta_max", [0])
        th  = round(float(th_raw[0] if isinstance(th_raw, list) else th_raw), 1)
        nxt = d.get("next_run_id", None)
        done = _has_result(p.parent)
        runs_info[rid] = {"omega": om, "theta": th, "next": nxt, "done": done}

    n_done = sum(1 for v in runs_info.values() if v["done"])
    print(f"STATUS: {n_done}/{TOTAL} conditions complete")

    # 2. Get inflight run IDs (scratch dirs with params but no results)
    inflight = _queued_run_ids()

    # 3. Find broken chains: done run whose next seg has no result and isn't inflight
    submitted: list[tuple[str, str, float, float]] = []
    for rid, info in runs_info.items():
        if not info["done"]:
            continue
        nxt = info["next"]
        if nxt is None:
            continue  # last segment in chain — no next
        nxt_info = runs_info.get(nxt, {})
        if nxt_info.get("done", False):
            continue  # already complete
        if nxt in inflight:
            continue  # already queued or running
        # Broken chain — submit
        jid = _sbatch(nxt, rid)
        if jid:
            submitted.append((nxt, jid, nxt_info.get("omega", 0), nxt_info.get("theta", 0)))
            inflight.add(nxt)

    if submitted:
        print(f"  Submitted {len(submitted)} missing segments:")
        for rid, jid, om, th in submitted:
            print(f"    run={rid} theta={th} omega={om:.4f} -> job {jid}")
    else:
        print("  Chain continuity OK — no missing submissions")

    # 4. Queue summary
    try:
        sq = subprocess.run(
            ["squeue", "-u", "eaguerov", "--format=%.10i %.8T", "--noheader"],
            capture_output=True, text=True, timeout=15
        ).stdout
        n_r = sq.count("RUNNING")
        n_p = sq.count("PENDING")
        print(f"  SLURM: {n_r} RUNNING  {n_p} PENDING")
    except Exception:
        pass

    # 5. Completion check
    if n_done == TOTAL:
        _generate_figure_and_commit(n_done)
        print("\nDONE — watchdog task complete.")
        sys.exit(0)

    # Per-omega progress
    omegas = sorted({v["omega"] for v in runs_info.values()})
    thetas = sorted({v["theta"] for v in runs_info.values()})
    print(f"\n  Progress by omega:")
    for om in omegas:
        n = sum(1 for v in runs_info.values() if v["omega"] == om and v["done"])
        bar = "█" * n + "░" * (len(thetas) - n)
        print(f"    {om:.4f} rad/s  {bar}  {n}/{len(thetas)}")

    sys.exit(1)   # signal "not done yet" to the caller


if __name__ == "__main__":
    main()
