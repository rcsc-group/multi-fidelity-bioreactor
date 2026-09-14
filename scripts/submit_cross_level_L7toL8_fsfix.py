"""L7->L8 cross-level pilot, re-run with the fs-before-refine ordering fix
(diary.md 2026-09-14). Generality check: the fix was verified to close the
L6->L7 gap from +29.7%/+219.8% (tau/ediss) to +0.04%/+0.01%; this pair
previously showed ~49% and has not been retested.

Same parameters as scripts/submit_cross_level_L7toL8.py, except the binary
(cleaned source with the fix, no smoothing hack) and the run_id. Runs on
the normal batch queue -- mbessa-condo's scoped authorization ended
2026-09-11 (feedback_mbessa_condo_scope.md).

Usage:
    uv run python scripts/submit_cross_level_L7toL8_fsfix.py
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
    "run_id": "crosslevel_L7toL8_fsfix", "fidelity": 8,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 60 * 0.607329, "n_mix_cycles": 60,
    "t_checkpoint": t_dump,
    "omega_b_prev": omega, "theta_max_prev": [7.0, 0.0, 0.0],
    "phi_angular_prev": [0.0, 0.0, 0.0], "amplitude_h_prev": [0.0, 0.0, 0.0],
    "phi_horizontal_prev": [0.0, 0.0, 0.0], "omega_h_prev": 0.0,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-crosslevel-clean",
}
job = submit_slurm(params, project_root=root, runs_root=root / "runs",
                   walltime="04:00:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   checkpoint=str((src_dir / "checkpoint.dump").resolve()),
                   cpus=1, ntasks=8, mem="4G")
print(f"L7->L8 cross-level pilot (fs fix): job={job}  run_id=crosslevel_L7toL8_fsfix")
