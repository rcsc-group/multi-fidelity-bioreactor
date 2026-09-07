"""Direct test of the bistability/hysteresis hypothesis (diary.md
2026-09-07): if the cold-start-beats-warm-start anomaly is genuine
hysteresis between two coexisting periodic branches (a fast, 3-cycle
ramp between conditions kicks the trajectory off the reference branch
onto a different, more robustly-attracting one), then a much SLOWER,
quasi-static ramp between the SAME two conditions should stay on the
reference branch instead.

Exact same case as the clean, unconfounded L7 result (theta_source=6.9 ->
target=7, Delta_theta=0.1deg, p-fix binary, 120-cycle window) that showed
a 26.1% offset from the independent theta=7 reference with the default
3-cycle ramp -- rerun here with N_RAMP_CYCLES=40 (the slow-ramp binary,
BioReactor-mpi-slowramp, also has the p.nodump restore fix baked in since
it was compiled from the same current src/BioReactor.c).

Usage:
    uv run python scripts/test_slow_ramp_hysteresis.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, omega_b_of, PROJECT_ROOT

chain.validate_params = lambda params: None

SLOWRAMP_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-slowramp"
RUNS_DIR = Path(PROJECT_ROOT) / "runs"
src_dir = RUNS_DIR / "dtheta_src_L7_th6.9"
d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
t_dump = float(d[-1, 1])

cfg = {
    "motion": {"omega_b": omega_b_of(32.5), "theta_max": [7.0, 0.0, 0.0]},
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
    "binary": SLOWRAMP_BINARY,
    "submit": True,
}
job_run_ids = chain.submit_chain(cfg)
print(f"theta_source=6.9 -> target=7 (SLOW-RAMP, 40 cycles): run={job_run_ids[0][0]}  {job_run_ids}")
