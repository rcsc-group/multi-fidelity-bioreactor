"""Tiny smoke test for cross-level (L6->L7) warm-starting before any real
compute is committed to it. 1-2 cycles only: does adapt_wavelet actually
refine to depth=7 with the right cell count, does the run survive it, does
anything blow up. Submitted under mbessa-condo -- explicitly, narrowly
authorized by the user for this experiment only (2026-09-10); never the
default (see feedback_mbessa_condo_scope.md).

Source: settling_baseline_L6 (32.5rpm/theta=7, already converged --
n_mix_cycles=80, so this is not a fresh-from-rest checkpoint).

Usage:
    uv run python scripts/smoke_cross_level_warmstart.py
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
    "run_id": "smoke_crosslevel_L6toL7", "fidelity": 7,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 2 * 0.607329, "n_mix_cycles": 2,
    "t_checkpoint": t_dump,
    "omega_b_prev": omega, "theta_max_prev": [7.0, 0.0, 0.0],
    "phi_angular_prev": [0.0, 0.0, 0.0], "amplitude_h_prev": [0.0, 0.0, 0.0],
    "phi_horizontal_prev": [0.0, 0.0, 0.0], "omega_h_prev": 0.0,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-crosslevel",
}
job = submit_slurm(params, project_root=root, runs_root=root / "runs",
                   walltime="00:30:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   checkpoint=str((src_dir / "checkpoint.dump").resolve()),
                   cpus=1, ntasks=8, mem="4G",
                   account="mbessa-condo", qos="mbessa-condo")
print(f"cross-level smoke test: job={job}  run_id=smoke_crosslevel_L6toL7")
