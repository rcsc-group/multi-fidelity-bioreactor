"""Settling-time study, step 2 of 3: independent reference cold starts at
each of the 8 non-baseline (omega_b, theta_max) points, at L6/L7/L8 --
the "ground truth" each checkpointed transition gets compared against.
24 runs. See diary.md 2026-09-06.

Usage:
    uv run python scripts/submit_settling_study_references.py
"""
import sys

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm
from scripts.settling_study_grid import (
    LEVELS, OUT_POINTS, T_END, PROJECT_ROOT, TEMPLATE, MAINLINE_BINARY,
    GEOMETRY, FILL_LEVEL, omega_b_of, reference_run_id,
)

job_ids = []
for fidelity, ntasks in LEVELS.items():
    for p in OUT_POINTS:
        params = {
            "run_id": reference_run_id(fidelity, p),
            "fidelity": fidelity,
            "geometry": GEOMETRY,
            "fill_level": FILL_LEVEL,
            "n_harmonics": 1,
            "theta_max": [p["theta"], 0.0, 0.0],
            "phi_angular": [0.0, 0.0, 0.0],
            "omega_b": omega_b_of(p["rpm"]),
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
        job_ids.append((params["run_id"], job_id))
        print(f"L{fidelity} reference rpm={p['rpm']:g} theta={p['theta']:g}: job {job_id}")

print(f"\n{len(job_ids)} reference jobs submitted")
