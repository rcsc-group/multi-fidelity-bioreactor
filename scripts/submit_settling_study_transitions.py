"""Settling-time study, step 3 of 3: checkpoint transitions from each
level's baseline to each of the 8 out-points. 24 runs. Requires the
baselines (submit_settling_study_baselines.py) to have completed first --
reads each baseline's actual final t from its own shear_stress.dat rather
than assuming a value. See diary.md 2026-09-06.

Usage:
    uv run python scripts/submit_settling_study_transitions.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import (
    LEVELS, BASELINE, OUT_POINTS, CYCLES_PER_RUN, PROJECT_ROOT, MAINLINE_BINARY,
    GEOMETRY, FILL_LEVEL, omega_b_of, baseline_run_id,
)

chain.validate_params = lambda params: None

RUNS_DIR = Path(PROJECT_ROOT) / "runs"

for fidelity, ntasks in LEVELS.items():
    baseline_dir = RUNS_DIR / baseline_run_id(fidelity)
    results_path = baseline_dir / "results.json"
    shear_path = baseline_dir / "shear_stress.dat"
    if not shear_path.exists():
        print(f"L{fidelity}: baseline not done yet ({baseline_dir}), skipping")
        continue

    d = np.loadtxt(shear_path, skiprows=1)
    t_dump = float(d[-1, 1])  # actual final t reached, read from disk, not assumed
    print(f"L{fidelity} baseline actual final t = {t_dump:.4f}")

    for p in OUT_POINTS:
        cfg = {
            "motion": {"omega_b": omega_b_of(p["rpm"]), "theta_max": [p["theta"], 0.0, 0.0]},
            "fidelity": fidelity,
            "geometry": GEOMETRY,
            "fill_level": FILL_LEVEL,
            "n_mix_cycles": CYCLES_PER_RUN,        # unused when initial_checkpoint is set
            "n_transition_cycles": CYCLES_PER_RUN,
            "t_buffer": 0.0,
            "sweep": {"parameter": "omega_b", "values": [omega_b_of(p["rpm"])]},
            "initial_checkpoint": {
                "t_dump": t_dump,
                "omega_b": omega_b_of(BASELINE["rpm"]),
                "theta_max": [BASELINE["theta"], 0.0, 0.0],
                "checkpoint_path": str(baseline_dir / "checkpoint.dump"),
            },
            "mpi": True,
            "ntasks": ntasks,
            "mem_per_cpu": "4G",
            "walltime": "04:00:00",
            "binary": MAINLINE_BINARY,
            "submit": True,
        }
        job_run_ids = chain.submit_chain(cfg)
        print(f"L{fidelity} transition -> rpm={p['rpm']:g} theta={p['theta']:g}: {job_run_ids}")
