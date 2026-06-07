"""Submit verification run for the MPI movies_output fix.

Tests that BioReactor-mpi-video passes t_mix without hanging (16 ranks).
"""
import json
import subprocess
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
MPI_TEMPLATE = PROJECT_ROOT / "config" / "slurm_mpi_template.sh"
SCRATCH_BASE = Path("/oscar/scratch/eaguerov/mpi_runs")
RUNS_ROOT = PROJECT_ROOT / "runs"

PARAMS = {
    "omega_b": 3.93, "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_h": 0.0, "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
    "fill_level": 0.5, "fidelity": 7,
    "n_mix_cycles": 5,
    "t_end": 12.0,
    "run_id": "mpi_fix_16r",
    "_ntasks": 16,
    "_walltime": "00:20:00",
    "_mem": "2G",
}

run_dir = RUNS_ROOT / PARAMS["run_id"]
run_dir.mkdir(parents=True, exist_ok=True)
PARAMS["_canonical_run_dir"] = str(run_dir.resolve())
(run_dir / "params.json").write_text(json.dumps(PARAMS, indent=2))

scratch_dir = SCRATCH_BASE / PARAMS["run_id"]
scratch_dir.mkdir(parents=True, exist_ok=True)
scratch_params = scratch_dir / "params.json"
scratch_params.write_text(json.dumps(PARAMS, indent=2))

cmd = [
    "sbatch", "--no-requeue",
    "--job-name=mpi-fix-16r",
    "--ntasks=16",
    "--time=00:20:00",
    "--mem-per-cpu=2G",
    "--cpus-per-task=1",
    f"--export=NONE,PARAMS={scratch_params}",
    str(MPI_TEMPLATE.resolve()),
]
result = subprocess.run(cmd, capture_output=True, text=True, check=True)
match = re.search(r"(\d+)", result.stdout)
job_id = match.group(1) if match else "unknown"
print(f"Fix verification job: {job_id}  (run_id=mpi_fix_16r)")
print(f"Monitor: squeue -j {job_id}")
print(f"Check: tail /oscar/scratch/eaguerov/mpi_runs/mpi_fix_16r/logstats.dat")
print(f"Success if last t >= 12.0 (past t_mix≈3)")
