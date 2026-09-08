"""A SMALL genuine omega change: 30 -> 32.5 rpm (su = 0.923), theta=7, L7.

Companion to test_omega_discriminator.py. Large omega changes escape
(17.5->32.5 su=0.538 -> 1.258; 37.5->32.5 su=1.154 -> 1.260) with no phase
error to explain it. If small omega changes also escape, the escape is a
step function in "omega changed at all", pointing at the code path. If this
one converges, the escape is magnitude-dependent.

Usage:
    uv run python scripts/test_omega_small.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

src_dir = Path(PROJECT_ROOT) / "runs" / "omega_src_L7_rpm30_th7"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
cfg = {
    "motion": {"omega_b": omega_b_of(32.5), "theta_max": [7.0, 0.0, 0.0]},
    "fidelity": 7, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
    "n_mix_cycles": 120, "n_transition_cycles": 120, "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(32.5)]},
    "initial_checkpoint": {
        "t_dump": t_dump, "omega_b": omega_b_of(30.0),
        "theta_max": [7.0, 0.0, 0.0],
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
    "binary": "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix2",
    "submit": True, "exclude": "node2336",
}
print(f"rpm 30 -> 32.5 (su=0.923): {chain.submit_chain(cfg)}")
