"""Real cross-level warm-start pilot: L6 (converged) -> refine to L7 ->
run forward, compare settling against the existing L7 cold start
(kicktest_L7_th7 -- zero new compute for the control).

60 cycles, sized from REAL measured L7 throughput (job 6034591, 8 tasks:
120 cycles in 1240s real = 0.17 min/cycle) -- 60 cycles ~= 10.2 min real,
walltime set to 2h for large margin (this is cheap regardless). Submitted
under mbessa-condo, the user's explicit, narrowly-scoped exception for
this experiment only (2026-09-10) -- not reused beyond it.

Usage:
    uv run python scripts/submit_cross_level_pilot.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
src_dir = root / "runs" / "settling_baseline_L6"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
omega = 32.5 * 2 * 3.141592653589793 / 60.0

params = {
    "run_id": "crosslevel_L6toL7_pilot", "fidelity": 7,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 60 * 0.607329, "n_mix_cycles": 60,
    "t_checkpoint": t_dump,
    "omega_b_prev": omega, "theta_max_prev": [7.0, 0.0, 0.0],
    "phi_angular_prev": [0.0, 0.0, 0.0], "amplitude_h_prev": [0.0, 0.0, 0.0],
    "phi_horizontal_prev": [0.0, 0.0, 0.0], "omega_h_prev": 0.0,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-crosslevel",
}
job = submit_slurm(params, project_root=root, runs_root=root / "runs",
                   walltime="02:00:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   checkpoint=str((src_dir / "checkpoint.dump").resolve()),
                   cpus=1, ntasks=8, mem="4G",
                   account="mbessa-condo", qos="mbessa-condo")
print(f"cross-level pilot: job={job}  run_id=crosslevel_L6toL7_pilot")
