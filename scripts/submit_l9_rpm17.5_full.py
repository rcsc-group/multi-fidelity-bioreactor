"""L9 sweep, condition 1 of 9 (17.5rpm, first in the low-to-high warm-start
chain -- diary.md 2026-09-04). Full cold start to genuine convergence, at
32 ranks -- reduced from 64 (2026-09-04): the first 64-rank submission
(job 5823592) queued behind ~772 pending jobs cluster-wide (Exploratory
tier's lower scheduling priority + most nodes already at least partially
occupied), SLURM's own estimated start time was ~5 days out. Cancelled,
resubmitting at 32 to fit into partially-occupied nodes and clear the
queue faster, trading some per-job speed for a much shorter wait.

Target: ~40 total cycles (8.8 ramp + ~31 QSS) -- t_end=24.3 nondim.
At the smoke test's measured 29.4 min/cycle (16 ranks, unscaled), this
would take ~19.6h even with ZERO speedup from going to 32 ranks -- safely
inside the 48h Exploratory-tier walltime cap either way, so walltime=44h
is a comfortable margin, not a tight bound.

Usage:
    uv run python scripts/submit_l9_rpm17.5_full.py
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
    "run_id": "l9_sweep_rpm17.5",
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
    "t_end": 24.3,
    "n_mix_cycles": 80,
    "_binary": BINARY,
}
job_id = submit_slurm(
    params, project_root=PROJECT_ROOT, walltime="44:00:00",
    template=TEMPLATE, mem="4G", cpus=1, ntasks=32,
)
print(f"l9_sweep_rpm17.5: job {job_id} (32 ranks, t_end=24.3, walltime=44h)")
