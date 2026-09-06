"""Direct empirical check: does MPI rank count change results, for THIS
solver, not just inferred from the OpenMP-only tau-histogram race fix
(a648ca2) -- diary.md 2026-09-06. Reruns the exact same params as the
already-completed L9 smoke test (job 5806403, 16 ranks) at 8 ranks.

Usage:
    uv run python scripts/check_rank_invariance.py
"""
import sys

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

PROJECT_ROOT = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"
TEMPLATE = f"{PROJECT_ROOT}/config/slurm_mpi_template.sh"
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-rampmatch"

params = {
    "run_id": "rank_check_l9_rpm17.5_8ranks",
    "fidelity": 9,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
    "fill_level": 0.5,
    "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": 1.832595714594046,  # 17.5rpm, same as smoke_l9_rpm17.5
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 2.0,
    "n_mix_cycles": 80,
    "_binary": BINARY,
}
job_id = submit_slurm(
    params, project_root=PROJECT_ROOT, walltime="03:00:00",
    template=TEMPLATE, mem="4G", cpus=1, ntasks=8,
)
print(f"rank_check (8 ranks): job {job_id}")
