"""Test A at L7 (diary.md 2026-09-07): does the persistent-offset effect
shrink toward zero at small-enough Delta_theta, same as it did at L6
(Test A: 4.1->4.1 gave 1.1%, 4.1->4.05 gave 4.9%), or does it stay large
even at very small steps (much harder to explain as anything but a bug)?

Same identical-checkpoint design as the L6 Test A: ALL restarts share the
EXACT SAME source checkpoint (dtheta_src_L7_th6.9, already computed for
the L6-vs-L7 resolution check), varying only the requested target:
  (a) 6.9 -> 6.9   (trivial, Delta_theta=0)
  (b) 6.9 -> 6.92  (Delta_theta=0.02)
  (c) 6.9 -> 6.95  (Delta_theta=0.05)

Usage:
    uv run python scripts/test_A_L7.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import PROJECT_ROOT, GEOMETRY, FILL_LEVEL, MAINLINE_BINARY, omega_b_of

chain.validate_params = lambda params: None

RUNS_DIR = Path(PROJECT_ROOT) / "runs"
src_dir = RUNS_DIR / "dtheta_src_L7_th6.9"
d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
t_dump = float(d[-1, 1])

CASES = [("6.9", 6.9), ("6.92", 6.92), ("6.95", 6.95)]

for tag, target_theta in CASES:
    cfg = {
        "motion": {"omega_b": omega_b_of(32.5), "theta_max": [target_theta, 0.0, 0.0]},
        "fidelity": 7,
        "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL,
        "n_mix_cycles": 60,
        "n_transition_cycles": 60,
        "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(32.5)]},
        "initial_checkpoint": {
            "t_dump": t_dump,
            "omega_b": omega_b_of(32.5),
            "theta_max": [6.9, 0.0, 0.0],
            "checkpoint_path": str(src_dir / "checkpoint.dump"),
        },
        "mpi": True,
        "ntasks": 8,
        "mem_per_cpu": "4G",
        "walltime": "06:00:00",
        "binary": MAINLINE_BINARY,
        "submit": True,
    }
    job_run_ids = chain.submit_chain(cfg)
    print(f"theta_source=6.9 -> target={target_theta:g}: run={job_run_ids[0][0]}  {job_run_ids}")
