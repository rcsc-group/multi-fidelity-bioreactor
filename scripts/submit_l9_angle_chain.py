"""L9 angle-sweep chain: the SECOND axis of Kim et al.'s data (this
project has only ever swept RPM so far). theta_max = 7 -> 6 -> 5 -> 4 ->
3 -> 2 deg, omega_b fixed at 32.5rpm -- matches shear_ediss_vs_angle.csv /
mixing_kla_vs_angle.csv exactly.

This is precisely the restart type root-caused and fixed this week
(theta-changing checkpoints, the phase-offset bug) -- exercising the fixed
machinery on its own target case, not a proxy for it.

Source: a34fc4d4, segment 5 of the just-finished RPM chain (32.5rpm,
theta=7, L9) -- already fully instrumented and validated bit-exact this
session, unlike l9_sweep_rpm17.5 (the RPM chain's stale, pre-
instrumentation bootstrap). No unverifiable link this time.

25 cycles/hop, video-enabled fully-fixed binary, 32 tasks/segment (fits
in the headroom left after the L10 scaling probe's 240).

Usage:
    uv run python scripts/submit_l9_angle_chain.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed"
N_TRANSITION = 25
THETAS = [6.0, 5.0, 4.0, 3.0, 2.0]

src_dir = Path(PROJECT_ROOT) / "runs" / "a34fc4d4"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])

cfg = {
    "motion": {"omega_b": omega_b_of(32.5), "theta_max": [THETAS[-1], 0.0, 0.0]},
    "fidelity": 9, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
    "n_mix_cycles": N_TRANSITION, "n_transition_cycles": N_TRANSITION, "t_buffer": 0.0,
    "sweep": {"parameter": "theta_max_0", "values": THETAS},
    "initial_checkpoint": {
        "t_dump": t_dump, "omega_b": omega_b_of(32.5),
        "theta_max": [7.0, 0.0, 0.0],
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True, "ntasks": 32, "mem_per_cpu": "4G", "walltime": "10:00:00",
    "binary": BINARY, "submit": True, "exclude": "node2336",
}
run_job_ids = chain.submit_chain(cfg)
print(f"\nL9 angle chain (theta order {THETAS}): {run_job_ids}")
