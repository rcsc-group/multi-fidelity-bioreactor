"""Which variable actually causes the branch jump: the ramp code path, or
the restored state?

analyze_branch_separation.py (L7, 32.5rpm, target theta=7deg, last 20 of
120 cycles, tau_mean per-cycle peak / cold-start reference):

    Delta_theta = 0      ->  0.996
    Delta_theta = 0.1    ->  1.255
    Delta_theta = 3      ->  1.257
    Delta_theta = 5      ->  1.253
    Delta_theta = 0.1, 40-cycle ramp -> 1.250

Identical destination (+/-0.3%, inside the 1.5-1.8% cycle-to-cycle spread)
across a 50x range of Delta_theta and a 13x range of ramp duration, with
Delta_theta=0 exactly on reference. That is a switch, not a response -- but
in the existing runs the switch is confounded: Delta_theta=0 is also the
only case restarted from the TARGET condition's own checkpoint, so
"restored state is already on the target attractor" and "the ramp branch is
a no-op" flip together.

These two runs break the confound. The ramp is algebraically a no-op iff
theta_max_prev == theta_max (it only ever feeds Ak/phk/omega_h/Ah/phh), so:

  B  ramp ON, state on-target:  restart from the theta=7 checkpoint, but
     declare theta_max_prev = 6.9999999. Delta_theta = 1e-7 -- physically
     nil, U_bio rescale differs by 1e-8 -- yet every ramp expression is
     now live.
  C  ramp OFF, state off-target: restart from the theta=6.9 checkpoint,
     but declare theta_max_prev = 7.0. No ramp, and a state that is a
     genuine ~0.3% away from the target attractor.

Predictions:
  B -> 1.25 and C -> 1.00  : the ramp code path is the cause, and the
      magnitude of the perturbation it is asked to apply is irrelevant.
  B -> 1.00 and C -> 1.25  : the ramp is innocent; the restored state is,
      which would contradict the kick test (a 1e-4 velocity perturbation
      on the same condition decayed over 120 cycles).
  both 1.25 / both 1.00    : neither variable alone; look elsewhere.

Usage:
    uv run python scripts/test_ramp_path_vs_state.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

FIXED_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-testc"
RPM, TARGET_THETA, FIDELITY, CYCLES = 32.5, 7.0, 7, 120
runs_dir = Path(PROJECT_ROOT) / "runs"

CASES = [
    ("B_ramp_on_state_on_target",  "settling_baseline_L7", 6.9999999),
    ("C_ramp_off_state_off_target", "dtheta_src_L7_th6.9",  7.0),
]

for label, source_run, theta_prev in CASES:
    src_dir = runs_dir / source_run
    t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
    cfg = {
        "motion": {"omega_b": omega_b_of(RPM), "theta_max": [TARGET_THETA, 0.0, 0.0]},
        "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
        "n_mix_cycles": CYCLES, "n_transition_cycles": CYCLES, "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(RPM)]},
        "initial_checkpoint": {
            "t_dump": t_dump, "omega_b": omega_b_of(RPM),
            "theta_max": [theta_prev, 0.0, 0.0],
            "checkpoint_path": str(src_dir / "checkpoint.dump"),
        },
        "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
        "binary": FIXED_BINARY, "submit": True, "exclude": "node2336",
    }
    job_run_ids = chain.submit_chain(cfg)
    print(f"{label}: source={source_run} theta_prev={theta_prev!r} "
          f"-> run={job_run_ids[0][0]}  {job_run_ids}")
