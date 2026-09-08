"""Does a SECOND chained segment stay on phase?

The phase fix shifts t_dump_checkpoint by t_phase_offset so that a restarted
segment writes its own checkpoint at a zero-crossing of ITS forcing, rather
than at an integer multiple of T_per_st that its offset has moved off-phase.
Without that, the fix would repair segment 1 and hand segment 2 exactly the
same bug one link down. That shift has never been exercised: every
verification run so far is a single segment.

Chain: settling_ref_L7_rpm32.5_th2 (theta=2, the largest phase offset in
the sweep at +116.6deg) -> segment 0 at theta=7 -> segment 1 at theta=7.
Segment 1 changes nothing, so it is a trivial restart from segment 0's own
dump; it can only go wrong if segment 0's dump is off-phase.

  both segments ~1.00 -> the shift works and chains are safe
  segment 0 ~1.00, segment 1 ~1.25 -> t_dump_checkpoint shift is wrong and
      the bug survives one link down

Usage:
    uv run python scripts/test_chain_two_segments.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

RPM, TARGET, FIDELITY, CYCLES = 32.5, 7.0, 7, 90
src_dir = Path(PROJECT_ROOT) / "runs" / "settling_ref_L7_rpm32.5_th2"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])

cfg = {
    "motion": {"omega_b": omega_b_of(RPM), "theta_max": [TARGET, 0.0, 0.0]},
    "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
    "n_mix_cycles": CYCLES, "n_transition_cycles": CYCLES, "t_buffer": 0.0,
    # two identical sweep points -> two chained segments
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(RPM), omega_b_of(RPM)]},
    "initial_checkpoint": {
        "t_dump": t_dump, "omega_b": omega_b_of(RPM),
        "theta_max": [2.0, 0.0, 0.0],
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
    "binary": "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix",
    "submit": True, "exclude": "node2336",
}
print(f"two-segment chain: {chain.submit_chain(cfg)}")
