"""Does an omega_b-changing restart converge to the cold start, and faster?

Two gaps left open by the theta-only verification.

CORRECTNESS. omega_b restarts should be immune to the phase bug by
derivation: T_per_st = T_per/T_bio is independent of omega_b, because
U_bio ~ omega_b and so T_bio ~ 1/omega_b. That is why the original code
survived the omega_b sweeps it was written for, and why theta sweeps broke
it. Derivation, not experiment -- so measure it.

These are also by far the hardest test of the su velocity rescale, which the
A/B found inert at up to 14%: 17.5 -> 32.5 rpm is su = 0.538, a 46% rescale
of u (and 71% of p, which goes as su^2).

SPEED. On tau the theta-changing warm starts saved only single-digit cycles,
partly because they pay the same 3-cycle forcing ramp as a cold start. An
omega_b change alters the flow scale rather than its geometry, so it is a
genuinely different question.

Target is 32.5 rpm / theta=7 / L7 both ways, so the reference is the same
kicktest_L7_th7 cold start used throughout.

Usage:
    uv run python scripts/test_omega_restart.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

TARGET_RPM, THETA, FIDELITY, CYCLES = 32.5, 7.0, 7, 120
runs_dir = Path(PROJECT_ROOT) / "runs"

for src_rpm in (17.5, 37.5):
    src_dir = runs_dir / f"settling_ref_L7_rpm{src_rpm:g}_th7"
    t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
    cfg = {
        "motion": {"omega_b": omega_b_of(TARGET_RPM), "theta_max": [THETA, 0.0, 0.0]},
        "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
        "n_mix_cycles": CYCLES, "n_transition_cycles": CYCLES, "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(TARGET_RPM)]},
        "initial_checkpoint": {
            "t_dump": t_dump, "omega_b": omega_b_of(src_rpm),
            "theta_max": [THETA, 0.0, 0.0],
            "checkpoint_path": str(src_dir / "checkpoint.dump"),
        },
        "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
        "binary": "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix2",
        "submit": True, "exclude": "node2336",
    }
    su = src_rpm / TARGET_RPM
    print(f"rpm {src_rpm:g} -> {TARGET_RPM:g} (su={su:.3f}): {chain.submit_chain(cfg)}")
