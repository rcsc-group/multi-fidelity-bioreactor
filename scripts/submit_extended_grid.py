"""Submit extended MPI checkpoint-restart grid (H18–H25).

Covers axes that H12–H17 left open:
  - Non-rectangular geometry (n=2, fid=4)
  - Fill levels 0.3 / 0.6 / 0.7
  - Theta 5° (vs. 7° baseline)
  - Level-8 fidelity (highest available)
  - kLa agreement test: restart from completed l7 vs. original kLa_25

Usage:
    uv run python scripts/submit_extended_grid.py [--dry-run]
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

# Verified source checkpoints (identified from sweep scratch dirs)
CKPT = {
    "n2_fid4":   SCRATCH_BASE / "001e47bf" / "checkpoint.dump",   # n=2, fid=4
    "l7_f03_t7": SCRATCH_BASE / "6439dcba" / "checkpoint.dump",   # fid=7 fill=0.3 theta=7 omega=1.5708
    "l7_f07_t7": SCRATCH_BASE / "08832477" / "checkpoint.dump",   # fid=7 fill=0.7 theta=7 omega=2.618
    "l7_f06_t7": SCRATCH_BASE / "3438372e" / "checkpoint.dump",   # fid=7 fill=0.6 theta=7 omega=1.5708
    "l7_f05_t5": SCRATCH_BASE / "3c940365" / "checkpoint.dump",   # fid=7 fill=0.5 theta=5 omega=3.1416
    "l8_t4":     SCRATCH_BASE / "29dc4215" / "checkpoint.dump",   # fid=8 fill=0.5 theta=4 omega=1.5708
    "l8_t7":     SCRATCH_BASE / "63062937" / "checkpoint.dump",   # fid=8 fill=0.5 theta=7 omega=1.5708
    "l7_kla_ref": SCRATCH_BASE / "1f8c7d4b" / "checkpoint.dump",  # fid=7 fill=0.3 omega=1.8326 kLa_25=0.139
}

# t_checkpoint values (last t from logstats, period-aligned)
T_CKPT = {
    "n2_fid4":    13.3,
    "l7_f03_t7":  99.1,
    "l7_f07_t7":  198.8,
    "l7_f06_t7":  99.1,
    "l7_f05_t5":  96.3,
    "l8_t4":      95.4,
    "l8_t7":      99.1,
    "l7_kla_ref": 198.8,
}

ROCKING_BASE_N8 = {
    "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
    "fill_level": 0.5,
}

# n=2 geometry (elliptical bag, fill=0.7, theta=7)
ROCKING_N2 = {
    "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 2.0},
    "fill_level": 0.7,
}

GRID = [
    # ── geometry stress ──────────────────────────────────────────────────────
    {
        "id": "grid_H18", "hypothesis": "H18",
        "ckpt_key": "n2_fid4",
        "fidelity": 4, "ntasks": 2,
        "omega_b": 6.28, "n_mix_cycles": 5, "t_end": 10.0,
        "rocking": ROCKING_N2,
        "walltime": "00:15:00",
        "notes": "n=2 elliptical geometry (fid=4): curved embedded boundary, most fa≈0 cells at coarse level",
    },
    # ── fill level sweep ─────────────────────────────────────────────────────
    {
        "id": "grid_H19", "hypothesis": "H19",
        "ckpt_key": "l7_f03_t7",
        "fidelity": 7, "ntasks": 16,
        "omega_b": 1.5708, "n_mix_cycles": 5, "t_end": 12.0,
        "rocking": dict(ROCKING_BASE_N8, fill_level=0.3),
        "walltime": "00:30:00",
        "notes": "l7 fill=0.3 (lowest fill level) — free surface close to bottom wall",
    },
    {
        "id": "grid_H20", "hypothesis": "H20",
        "ckpt_key": "l7_f07_t7",
        "fidelity": 7, "ntasks": 16,
        "omega_b": 2.618, "n_mix_cycles": 5, "t_end": 12.0,
        "rocking": dict(ROCKING_BASE_N8, fill_level=0.7),
        "walltime": "00:30:00",
        "notes": "l7 fill=0.7 (highest fill level) — free surface close to top wall",
    },
    {
        "id": "grid_H21", "hypothesis": "H21",
        "ckpt_key": "l7_f06_t7",
        "fidelity": 7, "ntasks": 16,
        "omega_b": 1.5708, "n_mix_cycles": 5, "t_end": 12.0,
        "rocking": dict(ROCKING_BASE_N8, fill_level=0.6),
        "walltime": "00:30:00",
        "notes": "l7 fill=0.6",
    },
    # ── theta / omega stress ─────────────────────────────────────────────────
    {
        "id": "grid_H22", "hypothesis": "H22",
        "ckpt_key": "l7_f05_t5",
        "fidelity": 7, "ntasks": 16,
        "omega_b": 3.1416, "n_mix_cycles": 5, "t_end": 12.0,
        "rocking": dict(ROCKING_BASE_N8, theta_max=[5.0, 0.0, 0.0]),
        "walltime": "00:30:00",
        "notes": "l7 theta=5 omega=3.1416 — different angle and higher frequency",
    },
    # ── l8 fidelity (highest available) ──────────────────────────────────────
    {
        "id": "grid_H23", "hypothesis": "H23",
        "ckpt_key": "l8_t4",
        "fidelity": 8, "ntasks": 16,
        "omega_b": 1.5708, "n_mix_cycles": 5, "t_end": 12.0,
        "rocking": dict(ROCKING_BASE_N8, theta_max=[4.0, 0.0, 0.0]),
        "walltime": "01:30:00",
        "notes": "l8 theta=4 — highest production fidelity, 4096 cells",
    },
    {
        "id": "grid_H24", "hypothesis": "H24",
        "ckpt_key": "l8_t7",
        "fidelity": 8, "ntasks": 16,
        "omega_b": 1.5708, "n_mix_cycles": 5, "t_end": 12.0,
        "rocking": ROCKING_BASE_N8,
        "walltime": "01:30:00",
        "notes": "l8 theta=7 — highest fidelity, production angle",
    },
    # ── kLa agreement test ───────────────────────────────────────────────────
    # Restart from end of a COMPLETED l7 run; compare kLa to original kLa_25=0.139
    # Note: tracers reset on restart (by design). Agreement means quasi-steady-state
    # flow is preserved across checkpoint boundary → same mass-transfer rate.
    {
        "id": "grid_H25", "hypothesis": "H25",
        "ckpt_key": "l7_kla_ref",
        "fidelity": 7, "ntasks": 16,
        "omega_b": 1.8326, "n_mix_cycles": 80, "t_end": 100.0,
        "rocking": dict(ROCKING_BASE_N8, fill_level=0.3),
        "walltime": "01:00:00",
        "notes": "kLa agreement: restart l7 fill=0.3 omega=1.8326 from t=198.8; "
                 "reference kLa_25=0.139 (same params, fresh run 1f8c7d4b); "
                 "expect |restart_kLa - 0.139| / 0.139 < 30%",
    },
]


def build_params(cell: dict) -> dict:
    p = dict(cell["rocking"])
    p.update({
        "fidelity": cell["fidelity"],
        "omega_b": cell["omega_b"],
        "n_mix_cycles": cell["n_mix_cycles"],
        "t_checkpoint": T_CKPT[cell["ckpt_key"]],
        "t_end": cell["t_end"],
        "run_id": cell["id"],
        "_binary": BINARY,
        "_ntasks": cell["ntasks"],
        "_walltime": cell["walltime"],
        "_mem": "2G",
    })
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

    dump_src = CKPT[cell["ckpt_key"]]
    if not dump_src.exists():
        raise FileNotFoundError(f"Checkpoint not found: {dump_src}")
    dump_dst = scratch_dir / "checkpoint.dump"
    shutil.copy2(dump_src, dump_dst)

    export_str = f"NONE,PARAMS={scratch_dir}/params.json,DUMP={dump_dst}"
    cmd = [
        "sbatch", "--no-requeue",
        f"--job-name=ext-{cell['hypothesis']}",
        f"--ntasks={cell['ntasks']}",
        f"--time={cell['walltime']}",
        "--mem-per-cpu=2G",
        "--cpus-per-task=1",
        f"--export={export_str}",
        str(MPI_TEMPLATE.resolve()),
    ]

    if dry_run:
        print(f"  [DRY] {run_id}: ntasks={cell['ntasks']} fid={cell['fidelity']} "
              f"n={cell['rocking']['geometry']['n']} fill={cell['rocking']['fill_level']} "
              f"theta={cell['rocking']['theta_max'][0]} omega={cell['omega_b']}")
        return "DRY_RUN"

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    job_id = re.search(r"(\d+)", result.stdout).group(1)
    print(f"  {run_id}: job {job_id}  (n={cell['rocking']['geometry']['n']} "
          f"fill={cell['rocking']['fill_level']} theta={cell['rocking']['theta_max'][0]} "
          f"fid={cell['fidelity']} {cell['ntasks']}r)")
    return job_id


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    import os
    if not os.path.exists(BINARY):
        raise FileNotFoundError(f"Binary not found: {BINARY}")
    print(f"Binary: {BINARY}  ✓")
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
    print("Check results:")
    print("  uv run python scripts/check_restart_run.py runs/grid_H1[89] runs/grid_H2*")


if __name__ == "__main__":
    main()
