"""Diagnostic: test whether MPI binary hangs at t_mix with 16 ranks vs 4 ranks.

Submits two short runs (n_mix_cycles=5 → t_mix≈3 nondim, t_end=12).
Both use BioReactor-mpi-video at fidelity=7.
Expected runtime: ~5 min each.

Usage:
    uv run python scripts/mpi_diag_submit.py
"""
import json
import subprocess
import re
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
MPI_TEMPLATE = PROJECT_ROOT / "config" / "slurm_mpi_template.sh"
SCRATCH_BASE = Path("/oscar/scratch/eaguerov/mpi_runs")
RUNS_ROOT = PROJECT_ROOT / "runs"
LOGS_DIR = PROJECT_ROOT / "logs"
LOGS_DIR.mkdir(exist_ok=True)

BASE_PARAMS = {
    "omega_b": 3.93, "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_h": 0.0, "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
    "fill_level": 0.5, "fidelity": 7,
    "n_mix_cycles": 5,   # t_mix ≈ 3 nondim (quick test through t_mix event)
    "t_end": 12.0,       # well past t_mix, short total run
    "_walltime": "00:20:00",
    "_mem": "2G",
}


def stage_and_submit(run_id: str, ntasks: int) -> str:
    params = dict(BASE_PARAMS, run_id=run_id, _ntasks=ntasks)

    # Write to canonical runs dir
    run_dir = RUNS_ROOT / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    params["_canonical_run_dir"] = str(run_dir.resolve())
    (run_dir / "params.json").write_text(json.dumps(params, indent=2))

    # Stage params to scratch
    scratch_dir = SCRATCH_BASE / run_id
    scratch_dir.mkdir(parents=True, exist_ok=True)
    scratch_params = scratch_dir / "params.json"
    scratch_params.write_text(json.dumps(params, indent=2))

    cmd = [
        "sbatch",
        "--no-requeue",
        f"--job-name=mpi-diag-{ntasks}r",
        f"--ntasks={ntasks}",
        "--time=00:20:00",
        "--mem-per-cpu=2G",
        "--cpus-per-task=1",
        f"--export=NONE,PARAMS={scratch_params}",
        str(MPI_TEMPLATE.resolve()),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    match = re.search(r"(\d+)", result.stdout)
    job_id = match.group(1) if match else "unknown"
    print(f"  {ntasks}-rank job submitted: {job_id}  (run_id={run_id})")
    return job_id


print("Submitting MPI diagnostic jobs...")
j16 = stage_and_submit("mpi_diag_16r", ntasks=16)
j4  = stage_and_submit("mpi_diag_4r",  ntasks=4)

print(f"\nMonitor with:")
print(f"  squeue -j {j16},{j4}")
print(f"\nAfter completion check:")
print(f"  tail /oscar/scratch/eaguerov/mpi_runs/mpi_diag_16r/logstats.dat")
print(f"  tail /oscar/scratch/eaguerov/mpi_runs/mpi_diag_4r/logstats.dat")
print(f"\nExpected: both should reach t=12.0 if no hang; "
      f"16-rank hangs at t≈3 if 16-rank is the issue.")
