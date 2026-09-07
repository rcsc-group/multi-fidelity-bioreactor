"""Extended-window re-run of the L7 Test A fine-Delta_theta cases
(diary.md 2026-09-07), p-fix binary, 120 cycles instead of 60 -- checking
whether the small-Delta_theta cases (0.02, 0.05deg) that appeared to
"converge to zero" in a 60-cycle window actually show the same delayed
dip-then-rise-to-a-plateau drift found in the main L7 sweep once given
enough cycles to see it.

Same identical-checkpoint design as the original Test A: both restarts
share the EXACT SAME source checkpoint (dtheta_src_L7_th6.9).

Usage:
    uv run python scripts/test_A_L7_extended.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, omega_b_of, PROJECT_ROOT

chain.validate_params = lambda params: None

FIXED_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-testc"
RUNS_DIR = Path(PROJECT_ROOT) / "runs"
src_dir = RUNS_DIR / "dtheta_src_L7_th6.9"
d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
t_dump = float(d[-1, 1])

CASES = [("6.92", 6.92), ("6.95", 6.95)]

for tag, target_theta in CASES:
    cfg = {
        "motion": {"omega_b": omega_b_of(32.5), "theta_max": [target_theta, 0.0, 0.0]},
        "fidelity": 7,
        "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL,
        "n_mix_cycles": 120,
        "n_transition_cycles": 120,
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
        "binary": FIXED_BINARY,
        "submit": True,
    }
    job_run_ids = chain.submit_chain(cfg)
    print(f"theta_source=6.9 -> target={target_theta:g} (p-FIXED, 120 cycles): run={job_run_ids[0][0]}  {job_run_ids}")
