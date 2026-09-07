"""Falsifiable test of the dAk/dt fix (diary.md 2026-09-07): rerun the
smallest anomalous case -- theta_source=4.1 -> target=4 (Delta_theta_max=
0.1deg), which under the OLD binary plateaued at tau_mean=1.72e-4, a 16.6%
offset from the independent cold-start reference (1.475e-4) -- using the
FIXED binary (BioReactor-mpi-fixed-ramp, dAk/dt and dphk/dt correctly
included in Th_d/Th_2d). If the fix is the root cause, this should now
converge to within a few percent of 1.475e-4 instead of plateauing at
+16.6%.

Reuses the existing dtheta_src_L6_th4.1 checkpoint as source (checkpoint
field data is unaffected by a C-code fix in the ramp event; only the
POST-restart dynamics change).

Usage:
    uv run python scripts/test_ramp_fix.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import PROJECT_ROOT, GEOMETRY, FILL_LEVEL, omega_b_of

chain.validate_params = lambda params: None

RUNS_DIR = Path(PROJECT_ROOT) / "runs"
FIXED_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-fixed-ramp"

theta_source = 4.1
src_dir = RUNS_DIR / "dtheta_src_L6_th4.1"
d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
t_dump = float(d[-1, 1])

cfg = {
    "motion": {"omega_b": omega_b_of(32.5), "theta_max": [4.0, 0.0, 0.0]},
    "fidelity": 6,
    "geometry": GEOMETRY,
    "fill_level": FILL_LEVEL,
    "n_mix_cycles": 60,
    "n_transition_cycles": 60,
    "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(32.5)]},
    "initial_checkpoint": {
        "t_dump": t_dump,
        "omega_b": omega_b_of(32.5),
        "theta_max": [theta_source, 0.0, 0.0],
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
print(f"theta_source=4.1 -> target=4 (FIXED binary): run={job_run_ids[0][0]}  {job_run_ids}")
