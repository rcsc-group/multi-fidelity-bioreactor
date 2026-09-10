"""Resubmit the L9 rpm chain's final segment (37.5rpm), killed by TIMEOUT
after 4 cycles short of its 25-cycle target due to the chain.py walltime-
stamping bug (diary.md 2026-09-10) -- fixed there, this resubmits the
segment itself, correctly this time.

Restarts from a281a16f (35rpm, segment 6, validated clean by
scripts/validate_run.py) exactly as the original chain intended, this
time via submit_chain's single-segment path (submit_slurm called
directly, not through the buggy self-submission annotation).

Usage:
    uv run python scripts/resubmit_rpm_chain_seg7.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

src_dir = Path(PROJECT_ROOT) / "runs" / "a281a16f"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])

cfg = {
    "motion": {"omega_b": omega_b_of(37.5), "theta_max": [7.0, 0.0, 0.0]},
    "fidelity": 9, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
    "n_mix_cycles": 25, "n_transition_cycles": 25, "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(37.5)]},
    "initial_checkpoint": {
        "t_dump": t_dump, "omega_b": omega_b_of(35.0),
        "theta_max": [7.0, 0.0, 0.0],
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True, "ntasks": 32, "mem_per_cpu": "4G", "walltime": "10:00:00",
    "binary": "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed",
    "submit": True, "exclude": "node2336",
}
print(f"rpm chain seg7 retry (37.5rpm, from a281a16f): {chain.submit_chain(cfg)}")
