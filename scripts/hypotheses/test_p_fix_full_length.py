"""Full-length verification of the p.nodump restore fix (diary.md
2026-09-07): rerun the exact case that showed a persistent 4.9% plateau
offset (theta_source=4.1 -> target=4.05, Delta_theta=0.05, L6 Test A),
60-cycle window, using the newly p-fixed binary, and check whether the
offset actually collapses now that pressure is genuinely restored instead
of silently reset to zero on every restart.

Usage:
    uv run python scripts/hypotheses/test_p_fix_full_length.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import scripts.chain as chain
from scripts.experiments.settling_study_grid import PROJECT_ROOT, GEOMETRY, FILL_LEVEL, omega_b_of
from scripts import paths  # noqa: E402  (machine-specific folders)

chain.validate_params = lambda params: None

RUNS_DIR = Path(PROJECT_ROOT) / "runs"
FIXED_BINARY = f"{paths.SCRATCH}/BioReactor-mpi-testc"  # now has BOTH the dAk/dt ramp fix and the p.nodump restore fix

src_dir = RUNS_DIR / "dtheta_src_L6_th4.1"
d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
t_dump = float(d[-1, 1])

cfg = {
    "motion": {"omega_b": omega_b_of(32.5), "theta_max": [4.05, 0.0, 0.0]},
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
        "theta_max": [4.1, 0.0, 0.0],
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
print(f"theta_source=4.1 -> target=4.05 (p-FIXED binary): run={job_run_ids[0][0]}  {job_run_ids}")
