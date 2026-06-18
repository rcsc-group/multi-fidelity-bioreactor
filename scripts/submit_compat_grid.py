"""Submit MPI checkpoint-restart compatibility grid (H12–H17).

Grid axes: fidelity (5, 6) × ranks (2, 16) × omega (same, 2×)
Each cell restarts from an existing checkpoint and must survive past t_mix
with finite oxy_liq_sum.

Pass/fail is determined by check_restart_run.py after each job completes.

Usage:
    uv run python scripts/submit_compat_grid.py [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
MPI_TEMPLATE = PROJECT_ROOT / "config" / "slurm_mpi_template.sh"
SCRATCH_BASE = Path("/oscar/scratch/eaguerov/mpi_runs")
RUNS_ROOT = PROJECT_ROOT / "runs"

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-stripped"

# Verified source checkpoints (pre-confirmed healthy runs)
CHECKPOINT_L5 = Path("/oscar/scratch/eaguerov/mpi_runs/smoke_l5_seg0/checkpoint.dump")
CHECKPOINT_L6 = Path("/oscar/scratch/eaguerov/mpi_runs/bench_f6_mpi_a/checkpoint.dump")

# t_checkpoint values matching the checkpoint dumps
T_CKPT_L5 = 97.294
T_CKPT_L6 = 147.7

# Common rocking parameters (match the source checkpoints)
ROCKING_BASE = {
    "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
    "fill_level": 0.5,
}

# 5 mix cycles → t_mix = t_checkpoint + 5*0.6081 ≈ t_checkpoint + 3.04
# t_end (relative) = 12 gives ≈9 nd time units of oxy data past t_mix
GRID = [
    {
        "id": "grid_H12",
        "hypothesis": "H12",
        "fidelity": 5,
        "ntasks": 2,
        "omega_b": 1.5708,
        "omega_b_prev": 0.0,  # same omega: no rescaling
        "t_checkpoint": T_CKPT_L5,
        "t_end": 12.0,
        "checkpoint": CHECKPOINT_L5,
        "walltime": "00:20:00",
        "notes": "l5, 2ranks, same-omega restart (baseline validation)",
    },
    {
        "id": "grid_H13",
        "hypothesis": "H13",
        "fidelity": 5,
        "ntasks": 2,
        "omega_b": 3.1416,
        "omega_b_prev": 1.5708,  # 2× omega: velocity rescaled by 0.5
        "t_checkpoint": T_CKPT_L5,
        "t_end": 12.0,
        "checkpoint": CHECKPOINT_L5,
        "walltime": "00:20:00",
        "notes": "l5, 2ranks, diff-omega (2x) restart: tests velocity rescaling",
    },
    {
        "id": "grid_H14",
        "hypothesis": "H14",
        "fidelity": 5,
        "ntasks": 16,
        "omega_b": 1.5708,
        "omega_b_prev": 0.0,
        "t_checkpoint": T_CKPT_L5,
        "t_end": 12.0,
        "checkpoint": CHECKPOINT_L5,
        "walltime": "00:20:00",
        "notes": "l5, 16ranks, same-omega restart: tests rank scaling",
    },
    {
        "id": "grid_H15",
        "hypothesis": "H15",
        "fidelity": 5,
        "ntasks": 16,
        "omega_b": 3.1416,
        "omega_b_prev": 1.5708,
        "t_checkpoint": T_CKPT_L5,
        "t_end": 12.0,
        "checkpoint": CHECKPOINT_L5,
        "walltime": "00:20:00",
        "notes": "l5, 16ranks, diff-omega: full stress test (most MPI ghosts + rescaling)",
    },
    {
        "id": "grid_H16",
        "hypothesis": "H16",
        "fidelity": 6,
        "ntasks": 2,
        "omega_b": 1.5708,
        "omega_b_prev": 0.0,
        "t_checkpoint": T_CKPT_L6,
        "t_end": 12.0,
        "checkpoint": CHECKPOINT_L6,
        "walltime": "01:00:00",
        "notes": "l6, 2ranks, same-omega: tests higher fidelity restart",
    },
    {
        "id": "grid_H17",
        "hypothesis": "H17",
        "fidelity": 6,
        "ntasks": 16,
        "omega_b": 1.5708,
        "omega_b_prev": 0.0,
        "t_checkpoint": T_CKPT_L6,
        "t_end": 12.0,
        "checkpoint": CHECKPOINT_L6,
        "walltime": "01:00:00",
        "notes": "l6, 16ranks, same-omega: full production-scale restart",
    },
]


def build_params(cell: dict) -> dict:
    p = dict(ROCKING_BASE)
    p.update({
        "fidelity": cell["fidelity"],
        "omega_b": cell["omega_b"],
        "n_mix_cycles": 5,
        "t_checkpoint": cell["t_checkpoint"],
        "t_end": cell["t_end"],
        "run_id": cell["id"],
        "_binary": BINARY,
        "_ntasks": cell["ntasks"],
        "_walltime": cell["walltime"],
        "_mem": "2G",
    })
    if cell["omega_b_prev"] > 0.0:
        p["omega_b_prev"] = cell["omega_b_prev"]
    run_dir = RUNS_ROOT / cell["id"]
    p["_canonical_run_dir"] = str(run_dir.resolve())
    return p


def stage_and_submit(cell: dict, dry_run: bool = False) -> str:
    run_id = cell["id"]
    params = build_params(cell)

    run_dir = RUNS_ROOT / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "params.json").write_text(json.dumps(params, indent=2))

    scratch_dir = SCRATCH_BASE / run_id
    scratch_dir.mkdir(parents=True, exist_ok=True)
    (scratch_dir / "params.json").write_text(json.dumps(params, indent=2))

    dump_src = cell["checkpoint"]
    dump_dst = scratch_dir / "checkpoint.dump"
    if not dump_src.exists():
        raise FileNotFoundError(f"Checkpoint not found: {dump_src}")
    shutil.copy2(dump_src, dump_dst)

    export_str = f"NONE,PARAMS={scratch_dir}/params.json,DUMP={dump_dst}"
    cmd = [
        "sbatch", "--no-requeue",
        f"--job-name=cgrid-{cell['hypothesis']}",
        f"--ntasks={cell['ntasks']}",
        f"--time={cell['walltime']}",
        "--mem-per-cpu=2G",
        "--cpus-per-task=1",
        f"--export={export_str}",
        str(MPI_TEMPLATE.resolve()),
    ]

    if dry_run:
        print(f"  [DRY] {run_id}: {' '.join(cmd)}")
        return "DRY_RUN"

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    job_id = re.search(r"(\d+)", result.stdout).group(1)
    print(f"  {run_id}: job {job_id}  ({cell['ntasks']} ranks, l{cell['fidelity']}, omega={cell['omega_b']})")
    return job_id


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    print(f"Binary: {BINARY}")
    import os
    if not os.path.exists(BINARY):
        raise FileNotFoundError(f"Binary not found: {BINARY}")
    print(f"  exists ✓")
    print()

    submitted = {}
    for cell in GRID:
        try:
            job_id = stage_and_submit(cell, dry_run=args.dry_run)
            submitted[cell["id"]] = job_id
        except Exception as e:
            print(f"  ERROR {cell['id']}: {e}")
            submitted[cell["id"]] = f"ERROR: {e}"

    print()
    print("Submitted jobs:")
    for run_id, job_id in submitted.items():
        print(f"  {run_id:20s} → {job_id}")
    print()
    print("Monitor:")
    print("  watch squeue -u eaguerov")
    print()
    print("Check results:")
    print("  uv run python scripts/check_restart_run.py runs/grid_H*")


if __name__ == "__main__":
    main()
