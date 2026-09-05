"""L9 smoke test: measure real wall-clock cost before committing to the
9-condition chained sweep (diary.md 2026-09-04). fidelity is a runtime
param (NN=1<<fidelity) -- the existing BioReactor-mpi-rampmatch binary
already supports L9, no rebuild needed.

Short run (t_end=2.0 nondim, ~3.3 cycles), 17.5rpm (first condition in
the low-to-high queue), 16 MPI ranks (matching the L8 sweep's rank
count as a baseline to scale from).

Usage:
    uv run python scripts/smoke_test_l9.py
"""
import math
import sys

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

PROJECT_ROOT = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"
TEMPLATE = f"{PROJECT_ROOT}/config/slurm_mpi_template.sh"
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-rampmatch"

omega_b = 17.5 * 2 * math.pi / 60.0
params = {
    "run_id": "smoke_l9_rpm17.5",
    "fidelity": 9,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
    "fill_level": 0.5,
    "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega_b,
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 2.0,
    "n_mix_cycles": 80,
    "_binary": BINARY,
}
job_id = submit_slurm(
    params, project_root=PROJECT_ROOT, walltime="02:00:00",
    template=TEMPLATE, mem="4G", cpus=1, ntasks=16,
)
print(f"smoke_l9_rpm17.5: job {job_id}")
