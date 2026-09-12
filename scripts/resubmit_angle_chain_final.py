"""Resubmit the angle-sweep chain's final segment (theta=2deg), killed by
the pre-fix walltime bug (same as the rpm chain's cbbd0063 -- this chain
was submitted before that fix landed). Restarts from 52e007f3 (theta=3,
validated clean), via submit_chain's direct (non-self-submitting) path.

Usage:
    uv run python scripts/resubmit_angle_chain_final.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

src_dir = Path(PROJECT_ROOT) / "runs" / "52e007f3"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])

cfg = {
    "motion": {"omega_b": omega_b_of(32.5), "theta_max": [2.0, 0.0, 0.0]},
    "fidelity": 9, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
    "n_mix_cycles": 25, "n_transition_cycles": 25, "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(32.5)]},
    "initial_checkpoint": {
        "t_dump": t_dump, "omega_b": omega_b_of(32.5),
        "theta_max": [3.0, 0.0, 0.0],
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True, "ntasks": 32, "mem_per_cpu": "4G", "walltime": "10:00:00",
    "binary": "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed",
    "submit": True, "exclude": "node2336",
}
print(f"angle chain final segment retry (theta=2, from 52e007f3): {chain.submit_chain(cfg)}")
