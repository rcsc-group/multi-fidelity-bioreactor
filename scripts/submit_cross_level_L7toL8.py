"""L7->L8 cross-level pilot: does the ~32%/coherent-structure discrepancy
found at L6->L7 persist at one finer resolution, or vanish like the
ω-checkpoint escape did going L7->L8 (diary.md 2026-09-08)? That precedent
is the direct reason to run this rather than assume either answer.

Source: kicktest_L7_th7 (139-cycle L7 cold start, most thoroughly
verified reference this session). Refine to L8, run 60 cycles (matching
the L6->L7 pilot's scope). Control: settling_baseline_L8 (already 80
cycles at the same condition -- zero new compute).

Sized from real measured L8 throughput (settling_baseline_L8's own job
5979997, 8 tasks: 80 cycles in 2671s real = 0.56 min/cycle) -- 60 cycles
~= 33.6 min real. mbessa-condo, continuing the same user-authorized
cross-level experiment.

Usage:
    uv run python scripts/submit_cross_level_L7toL8.py
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
    "run_id": "crosslevel_L7toL8_pilot", "fidelity": 8,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 60 * 0.607329, "n_mix_cycles": 60, "frames_per_period": 5,
    "t_checkpoint": t_dump,
    "omega_b_prev": omega, "theta_max_prev": [7.0, 0.0, 0.0],
    "phi_angular_prev": [0.0, 0.0, 0.0], "amplitude_h_prev": [0.0, 0.0, 0.0],
    "phi_horizontal_prev": [0.0, 0.0, 0.0], "omega_h_prev": 0.0,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-crosslevel-video",
}
job = submit_slurm(params, project_root=root, runs_root=root / "runs",
                   walltime="02:00:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   checkpoint=str((src_dir / "checkpoint.dump").resolve()),
                   cpus=1, ntasks=8, mem="4G",
                   account="mbessa-condo", qos="mbessa-condo")
print(f"L7->L8 cross-level pilot: job={job}  run_id=crosslevel_L7toL8_pilot")
