"""Delta-theta_max monotonicity study, DECREASING-direction mirror
(diary.md 2026-09-07): only 2 new cold-start sources needed for fine
resolution near the target (theta=2.1, 2.5) -- every other source
(3,4,5,6,6.5,6.9,7) already exists from the increasing-direction sweep and
its parent grid study, and gets reused directly by
submit_dtheta_transitions_decreasing.py.

Usage:
    uv run python scripts/submit_dtheta_sources_decreasing.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm
from scripts.settling_study_grid import T_END, GEOMETRY, FILL_LEVEL, PROJECT_ROOT, TEMPLATE, MAINLINE_BINARY, omega_b_of

RPM = 32.5
FIDELITY = 6
NTASKS = 8
NEW_THETA_SOURCES = [2.1, 2.5]


def source_run_id(theta: float) -> str:
    return f"dtheta_src_L6_th{theta:g}"


for theta in NEW_THETA_SOURCES:
    params = {
        "run_id": source_run_id(theta),
        "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
        "n_harmonics": 1, "theta_max": [theta, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": omega_b_of(RPM), "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": T_END, "n_mix_cycles": 80, "_binary": MAINLINE_BINARY,
    }
    job_id = submit_slurm(params, project_root=Path(PROJECT_ROOT), walltime="01:00:00",
                           template=Path(TEMPLATE), ntasks=NTASKS, mem="4G")
    print(f"L6 theta={theta:g} -> run={params['run_id']}  job={job_id}")
