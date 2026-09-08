"""Is the second attractor an L7 artifact, or real?

Everything in the omega investigation is L7. The branch the escaping
restarts land on is a thin wall film (uy_liq_max 4.3x, posY_max 1.65x at
constant mass and slightly LESS interfacial area) -- exactly the sort of
feature that may not survive refinement. This decides whether the su
threshold on omega checkpointing is a real physical limit or a discretisation
artifact.

Repeats the worst escaping case at L8: omega 17.5 -> 32.5, theta=7, which at
L7 lands at 1.2308 against its cold start. Reference is settling_baseline_L8
(a 25-cycle L8 cold start at 32.5rpm/theta=7); cold starts settle by ~cycle
7, so its late cycles are converged.

  ~1.00 at L8 -> the escaped branch does not survive refinement, the L7
      threshold is an artifact, and omega checkpointing is unconstrained
  ~1.23 at L8 -> the branch is real and large omega steps are genuinely
      unsafe; the staircase is the fix

Usage:
    uv run python scripts/test_omega_L8.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

src_dir = Path(PROJECT_ROOT) / "runs" / "settling_ref_L8_rpm17.5_th7"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
cfg = {
    "motion": {"omega_b": omega_b_of(32.5), "theta_max": [7.0, 0.0, 0.0]},
    "fidelity": 8, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
    "n_mix_cycles": 80, "n_transition_cycles": 80, "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(32.5)]},
    "initial_checkpoint": {
        "t_dump": t_dump, "omega_b": omega_b_of(17.5),
        "theta_max": [7.0, 0.0, 0.0],
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True, "ntasks": 16, "mem_per_cpu": "4G", "walltime": "24:00:00",
    "binary": "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix2",
    "submit": True, "exclude": "node2336",
}
print(f"L8 omega 17.5 -> 32.5: {chain.submit_chain(cfg)}")
