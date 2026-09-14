"""Real cross-level L6->L7 pilot re-test with fix attempt #8 (diary.md
2026-09-14): repeated (every-timestep, first 10 rocking cycles), not
one-shot, near-wall corrective smoothing, testing whether the checkerboard
is a weakly-damped mode that needs continuous suppression rather than a
one-time IC defect. Identical parameters to the original pilot
(scripts/submit_cross_level_pilot.py) and the csfix re-test, except the
binary points at the newly-built BioReactor-mpi-crosslevel-fix8.

Usage:
    uv run python scripts/submit_cross_level_fix8.py
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
    "run_id": "crosslevel_L6toL7_pilot_fix8", "fidelity": 7,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 60 * 0.607329, "n_mix_cycles": 60,
    "t_checkpoint": t_dump,
    "omega_b_prev": omega, "theta_max_prev": [7.0, 0.0, 0.0],
    "phi_angular_prev": [0.0, 0.0, 0.0], "amplitude_h_prev": [0.0, 0.0, 0.0],
    "phi_horizontal_prev": [0.0, 0.0, 0.0], "omega_h_prev": 0.0,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-crosslevel-fix8",
}
job = submit_slurm(params, project_root=root, runs_root=root / "runs",
                   walltime="02:00:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   checkpoint=str((src_dir / "checkpoint.dump").resolve()),
                   cpus=1, ntasks=8, mem="4G")
print(f"cross-level pilot (fix #8, repeated suppression): job={job}  run_id=crosslevel_L6toL7_pilot_fix8")
