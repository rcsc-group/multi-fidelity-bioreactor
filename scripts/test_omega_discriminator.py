"""Is the omega_b escape a step function (any omega change) or magnitude-
dependent (a genuinely large perturbation leaving the basin)?

Both omega restarts tested land on the ~1.25 branch, but there is no phase
error: T_per_st is independent of omega_b (0.607329 at 17.5/32.5/37.5) and
the solver's computed offset is ~0. Magnitude does not explain it either --
theta 2->7 starts 53% off (ratio 0.473 at cycle 0) and converges cleanly,
while omega 37.5->32.5 starts 33% off and escapes. So this is a separate
defect from the phase bug.

Same 2x2 logic that cracked the theta case:

  W-B  omega code path live, state on target: restart from the 32.5rpm dump
       at 32.5rpm, but declare omega_b_prev = 32.4rpm. su = 0.997, a 0.3%
       velocity rescale of an already-correct state -- and 0.3% is a
       perturbation the kicknoise run measured decaying to nothing. If this
       escapes, the omega code path itself is broken and magnitude is
       irrelevant, exactly as the ramp/phase step function was.

  W-src  a nearby-rpm source (30rpm, theta=7, L7, 25 cycles) so that a
       genuinely SMALL omega change (30 -> 32.5, su = 0.923) can be tested
       against the large ones. Submitted here; the restart follows once it
       finishes.

Usage:
    uv run python scripts/test_omega_discriminator.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.simulate import submit_slurm
from scripts.settling_study_grid import (GEOMETRY, FILL_LEVEL, PROJECT_ROOT,
                                         TEMPLATE, omega_b_of)

chain.validate_params = lambda params: None

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix2"
ROOT = Path(PROJECT_ROOT)
runs_dir = ROOT / "runs"
T_PER_ND = 0.607329

# --- W-B: omega code path live, state already on target -------------------
src_dir = runs_dir / "settling_baseline_L7"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
cfg = {
    "motion": {"omega_b": omega_b_of(32.5), "theta_max": [7.0, 0.0, 0.0]},
    "fidelity": 7, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
    "n_mix_cycles": 120, "n_transition_cycles": 120, "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(32.5)]},
    "initial_checkpoint": {
        "t_dump": t_dump, "omega_b": omega_b_of(32.4),   # the only lie
        "theta_max": [7.0, 0.0, 0.0],
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
    "binary": BINARY, "submit": True, "exclude": "node2336",
}
print(f"W-B (omega path live, su=0.997): {chain.submit_chain(cfg)}")

# --- W-src: a 30rpm source for a small genuine omega change ---------------
src = {
    "run_id": "omega_src_L7_rpm30_th7", "fidelity": 7,
    "geometry": GEOMETRY, "fill_level": FILL_LEVEL, "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega_b_of(30.0), "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 25 * T_PER_ND, "n_mix_cycles": 25, "_binary": BINARY,
}
job = submit_slurm(src, project_root=ROOT, walltime="02:00:00",
                   template=Path(TEMPLATE), cpus=1, ntasks=8, mem="4G",
                   exclude="node2336")
print(f"W-src (30rpm theta=7 cold start, 25 cycles): job={job}")
