"""Cheap L6 smoke test of the restart phase-offset fix before spending the
full L7 verification sweep on it.

The first attempt at this fix rescaled the global t from inside event init
and died with an FPE in dtnext() on every theta-changing restart (the
theta=0 case survived only because its factor is exactly 1). This runs the
same shape of case -- a theta-changing restart, so t_phase_offset != 0 --
for 5 cycles at L6, just to confirm it starts, steps and exits cleanly.

Usage:
    uv run python scripts/smoke_phase_fix.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
src_dir = root / "runs" / "dtheta_src_L6_th4.1"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
omega = 32.5 * 2 * 3.141592653589793 / 60.0

params = {
    "run_id": "smoke_phasefix", "fidelity": 6,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [4.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 3.0, "n_mix_cycles": 5,
    "t_checkpoint": t_dump,
    "omega_b_prev": omega, "theta_max_prev": [4.1, 0.0, 0.0],
    "phi_angular_prev": [0.0, 0.0, 0.0], "amplitude_h_prev": [0.0, 0.0, 0.0],
    "phi_horizontal_prev": [0.0, 0.0, 0.0], "omega_h_prev": 0.0,
    "_binary": (sys.argv[1] if len(sys.argv) > 1
                else "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix"),
}
job = submit_slurm(params, project_root=root, runs_root=root / "runs",
                   walltime="00:30:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   checkpoint=str((src_dir / "checkpoint.dump").resolve()),
                   cpus=1, ntasks=8, mem="4G", exclude="node2336")
print(f"smoke: job={job}  run_id=smoke_phasefix")
