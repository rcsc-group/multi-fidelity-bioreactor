"""Regression guard for the 2026-09-14 cross-level fix: confirm an ORDINARY
(same-resolution, non-cross-level) restart is unaffected by it.

The fix (recompute cs/fs before refine()) sits inside
`#if CROSS_LEVEL_WARMSTART`, so an ordinary restart should be untouched --
but that is a code-structure argument, not a measurement. This restarts
kicktest_L7_th7 into itself with the cleaned PLAIN build (no cross-level
flag) and should reproduce the established control value: ~0% gap vs the
cold start (the pre-fix control measured -0.11%).

Usage:
    uv run python scripts/submit_trivial_restart_regression.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
src_dir = root / "runs" / "kicktest_L7_th7"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
omega = 32.5 * 2 * 3.141592653589793 / 60.0

params = {
    "run_id": "trivial_restart_regression", "fidelity": 7,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 20 * 0.607329, "n_mix_cycles": 20,
    "t_checkpoint": t_dump,
    "omega_b_prev": omega, "theta_max_prev": [7.0, 0.0, 0.0],
    "phi_angular_prev": [0.0, 0.0, 0.0], "amplitude_h_prev": [0.0, 0.0, 0.0],
    "phi_horizontal_prev": [0.0, 0.0, 0.0], "omega_h_prev": 0.0,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-plain-clean",
}
job = submit_slurm(params, project_root=root, runs_root=root / "runs",
                   walltime="02:00:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   checkpoint=str((src_dir / "checkpoint.dump").resolve()),
                   cpus=1, ntasks=8, mem="4G")
print(f"trivial-restart regression: job={job}  run_id=trivial_restart_regression")
