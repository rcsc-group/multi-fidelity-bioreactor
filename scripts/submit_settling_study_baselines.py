"""Settling-time study, step 1 of 3: baseline cold starts (32.5rpm, 7deg)
at L6/L7/L8, mainline binary. These are the checkpoint sources for all
8 out-points at each level. See diary.md 2026-09-06.

Usage:
    uv run python scripts/submit_settling_study_baselines.py
"""
import sys

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm
from scripts.settling_study_grid import (
    LEVELS, BASELINE, T_END, PROJECT_ROOT, TEMPLATE, MAINLINE_BINARY,
    GEOMETRY, FILL_LEVEL, omega_b_of, baseline_run_id,
)

for fidelity, ntasks in LEVELS.items():
    params = {
        "run_id": baseline_run_id(fidelity),
        "fidelity": fidelity,
        "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL,
        "n_harmonics": 1,
        "theta_max": [BASELINE["theta"], 0.0, 0.0],
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": omega_b_of(BASELINE["rpm"]),
        "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0],
        "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": T_END,
        "n_mix_cycles": 80,
        "_binary": MAINLINE_BINARY,
    }
    job_id = submit_slurm(
        params, project_root=PROJECT_ROOT, walltime="04:00:00",
        template=TEMPLATE, mem="4G", cpus=1, ntasks=ntasks,
    )
    print(f"L{fidelity} baseline: job {job_id} (run_id={params['run_id']}, ntasks={ntasks})")
