"""Direct empirical check: does MPI rank count change results, for THIS
solver, not just inferred from the OpenMP-only tau-histogram race fix
(a648ca2) -- diary.md 2026-09-06.

Redone at L6 (2026-09-06, corrected): the first attempt reused the L9
smoke test's params (fidelity=9) -- unnecessarily expensive, since rank
invariance is a structural/algorithmic property of the MPI code
(deterministic domain decomposition, no data races -- the only race
found was OpenMP-specific) and shouldn't depend on grid resolution at
all. Testing at L6 gives equally valid evidence for a fraction of the
cost. Reuses runs/fig13a_l6_rpm32.5 (ntasks=4) as the comparison
reference -- reruns the exact same params at 8 ranks.

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
    "run_id": "rank_check_l6_rpm32.5_8ranks",
    "fidelity": 6,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
    "fill_level": 0.5,
    "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": 3.4033920413889422,  # 32.5rpm, same as fig13a_l6_rpm32.5 (ntasks=4)
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 20.0,
    "n_mix_cycles": 80,
    "_binary": BINARY,
}
job_id = submit_slurm(
    params, project_root=PROJECT_ROOT, walltime="00:30:00",
    template=TEMPLATE, mem="2G", cpus=1, ntasks=8,
)
print(f"rank_check L6 (8 ranks): job {job_id}")
