"""Can a large omega change be checkpointed by stepping it?

The omega escape is not a restart defect: the state conversion is correct,
g turned out to be irrelevant (leaving it, zeroing it and scaling it by su^2
all give 1.227-1.231 on the worst case), and the omega code path is clean at
small steps. What decides the outcome is the SIZE of the jump -- clean at
|su-1| <= 0.126, escaping at >= 0.154 -- i.e. the warm start lands outside
the reference basin when the step is too big.

If that reading is right, walking 17.5 -> 32.5 rpm in steps that each stay
under the threshold should arrive on the reference branch, where the single
jump lands at 1.2308. Steps here are ~+15% in rpm, i.e. su ~ 0.87 per
segment, the largest value measured clean:

    17.5 -> 20 -> 23 -> 26.5 -> 30 -> 32.5

Only the final segment's value is compared against the 32.5rpm cold start;
the intermediates just have to stay on their own branches.

Usage:
    uv run python scripts/test_omega_staircase.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

STEPS = [20.0, 23.0, 26.5, 30.0, 32.5]
src_dir = Path(PROJECT_ROOT) / "runs" / "settling_ref_L7_rpm17.5_th7"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])

prev = 17.5
for a, b in zip([17.5] + STEPS[:-1], STEPS):
    print(f"  step {a:g} -> {b:g} rpm: su = {a/b:.3f}")

cfg = {
    "motion": {"omega_b": omega_b_of(32.5), "theta_max": [7.0, 0.0, 0.0]},
    "fidelity": 7, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
    "n_mix_cycles": 60, "n_transition_cycles": 60, "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(r) for r in STEPS]},
    "initial_checkpoint": {
        "t_dump": t_dump, "omega_b": omega_b_of(17.5),
        "theta_max": [7.0, 0.0, 0.0],
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
    "binary": "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix2",
    "submit": True, "exclude": "node2336",
}
print(f"\nstaircase: {chain.submit_chain(cfg)}")
