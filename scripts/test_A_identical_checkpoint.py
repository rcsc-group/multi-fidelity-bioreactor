"""Test A -- identical-checkpoint isolation (diary.md 2026-09-07). Both
runs restart from the EXACT SAME source checkpoint (dtheta_src_L6_th4.1),
so restored fields at t=t_checkpoint are byte-identical; only the
requested TARGET theta_max differs:
  (a) theta_source=4.1 -> target=4.1  (trivial, Delta_theta=0 -- NEVER
      actually run before now; earlier trivial controls were 4->4 and
      2->2, never FROM this specific source)
  (b) theta_source=4.1 -> target=4.05 (Delta_theta=0.05, half the step
      that already showed a 15% persistent offset at Delta_theta=0.1)
Isolates whether the persistent-offset anomaly originates strictly in
what happens AFTER restore (ramp/forcing), with zero confound from which
checkpoint file was read. Uses the FIXED binary (dAk/dt-corrected
Th_d/Th_2d) since that's now the physically-correct formula, even though
it already failed to explain the theta_source=4.1->4.0 case.

Usage:
    uv run python scripts/test_A_identical_checkpoint.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import PROJECT_ROOT, GEOMETRY, FILL_LEVEL, omega_b_of

chain.validate_params = lambda params: None

RUNS_DIR = Path(PROJECT_ROOT) / "runs"
FIXED_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-fixed-ramp"
src_dir = RUNS_DIR / "dtheta_src_L6_th4.1"
d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
t_dump = float(d[-1, 1])

CASES = [("4.1", 4.1), ("4.05", 4.05)]

for tag, target_theta in CASES:
    cfg = {
        "motion": {"omega_b": omega_b_of(32.5), "theta_max": [target_theta, 0.0, 0.0]},
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
    print(f"theta_source=4.1 -> target={target_theta:g}: run={job_run_ids[0][0]}  {job_run_ids}")
