"""Delta-theta_max monotonicity study, phase 1 (diary.md 2026-09-07): cold-
start reference runs at intermediate theta_max values, needed as checkpoint
SOURCES for phase 2's warm-start transitions into a single fixed TARGET
(32.5rpm, theta=7deg -- reuses the existing settling_baseline_L6 as both the
target's own converged reference and, for the trivial Delta_theta=0 case,
its own checkpoint source).

Design (per user request, diary.md 2026-09-07): fix ONE target condition,
vary the warm-start SOURCE across a range of theta_max values approaching
the target from below, and check two things:
  (1) settling time decreases monotonically as Delta_theta_max -> 0
  (2) settling time is exactly 0 at Delta_theta_max = 0 (source == target,
      "already converged")
This directly tests whether the cold-start-beats-warm-start anomaly found
in the main settling study (theta-only edges, L6/L7/L8) is at least
*locally* sane -- a small enough perturbation should settle fast, even if
larger ones behave anomalously.

theta_max=2.0, 4.0, 7.0 already have converged L6 references
(settling_ref_L6_rpm32.5_th2, _th4, settling_baseline_L6) from the earlier
grid study -- only the intermediate values need new cold starts here.

Usage:
    uv run python scripts/submit_dtheta_sources.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm
from scripts.settling_study_grid import (
    T_END, GEOMETRY, FILL_LEVEL, PROJECT_ROOT, TEMPLATE, MAINLINE_BINARY, omega_b_of,
)

RPM = 32.5
FIDELITY = 6
NTASKS = 8

NEW_THETA_SOURCES = [1.0, 3.0, 5.0, 6.0, 6.5, 6.9]


def source_run_id(theta: float) -> str:
    return f"dtheta_src_L6_th{theta:g}"


for theta in NEW_THETA_SOURCES:
    params = {
        "run_id": source_run_id(theta),
        "fidelity": FIDELITY,
        "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL,
        "n_harmonics": 1,
        "theta_max": [theta, 0.0, 0.0],
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": omega_b_of(RPM),
        "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0],
        "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": T_END,
        "n_mix_cycles": 80,
        "_binary": MAINLINE_BINARY,
    }
    job_id = submit_slurm(
        params, project_root=Path(PROJECT_ROOT), walltime="01:00:00",
        template=Path(TEMPLATE), ntasks=NTASKS, mem="4G",
    )
    print(f"L6 theta={theta:g} -> run={params['run_id']}  job={job_id}")
