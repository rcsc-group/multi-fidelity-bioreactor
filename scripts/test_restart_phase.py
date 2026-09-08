"""Falsification test for the restart phase-error hypothesis.

Hypothesis (diary.md 2026-09-07 (4)). The nondimensional time unit is
T_bio = L_bio/U_bio and U_bio depends on tan(theta_max)
(src/BioReactor.c:339-341), so T_per_st = T_per/T_bio differs between
theta_max values at the same omega_b. A checkpoint is dumped at an integer
multiple of the SOURCE run's T_per_st (src/BioReactor.c:383-384), i.e. at a
zero-crossing of that run. Restart it under a different theta_max and the
same t is no longer an integer multiple of the TARGET's T_per_st, so the
forcing resumes at a nonzero phase while the restored velocity and interface
fields correspond to phase 0.

Measured phase errors line up perfectly with the measured escape:
  source theta 7.0 -> -1.9deg   -> tau/tau_ref 0.996  (clean)
  source theta 6.9 -> -30.0deg  -> 1.255
  source theta 4.0 -> -53.7deg  -> 1.257
  source theta 2.0 -> +116.6deg -> 1.253

Correlation is not causation, and the phase error is confounded with
theta_max changing at all. These two runs break that, using phi_angular to
move the forcing phase with no theta change anywhere:

  P1 INJECT -- the known-clean configuration (restart from the theta=7 dump,
     target 7, prev 7, phase error -1.9deg, measured 0.996) with
     phi_angular = phi_angular_prev set to put a -30deg error in by hand.
     Hypothesis predicts it ESCAPES to ~1.25 with zero theta change.

  P2 CANCEL -- the ordinary Delta_theta=0.1 restart (theta 6.9 -> 7, correct
     U_bio rescale, measured 1.255) with phi_angular = phi_angular_prev set
     to exactly cancel its -30deg error. Hypothesis predicts it comes back
     CLEAN at ~1.00, which is also the fix.

Both set phi_angular == phi_angular_prev, so dphk = 0 and the phase is a
constant offset, not something the ramp touches.

Usage:
    uv run python scripts/test_restart_phase.py
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

RPM, TARGET, FIDELITY, CYCLES = 32.5, 7.0, 7, 120
L_BIO, H_BIO = 0.25, 2 * 0.03575
ROOT = Path(PROJECT_ROOT)
runs_dir = ROOT / "runs"


def T_per_st(theta):
    """Nondimensional rocking period, exactly as src/BioReactor.c:339-344."""
    omega_b = omega_b_of(RPM)
    T_per = 2 * math.pi / omega_b
    V_bio = L_BIO / 4 * (H_BIO + 0.5 * L_BIO * math.tan(math.radians(theta)))
    U_bio = V_bio / (H_BIO * 0.5) / T_per
    return T_per / (L_BIO / U_bio)


def phase_error(t_dump, theta_target):
    """Forcing phase at t_dump, in radians, wrapped to (-pi, pi]."""
    n = t_dump / T_per_st(theta_target)
    return 2 * math.pi * (n - round(n))


# Exact restore-time t values, read from the runs' own post-restore
# diagnostics rather than from the coarser shear_stress.dat sampling.
def exact_t(run):
    for line in (runs_dir / run / "restart_diagnostic_post_restore.txt").read_text().splitlines():
        if line.startswith("t "):
            return float(line.split()[1])
    raise KeyError(run)


t_clean = exact_t("8fd81f04")     # restart from settling_baseline_L7 (theta 7)
t_esc = exact_t("46196ea2")       # restart from dtheta_src_L7_th6.9 (theta 6.9)
print(f"t at restore: clean {t_clean:.12f}, escaping {t_esc:.12f}")
print(f"native phase error: clean {math.degrees(phase_error(t_clean, TARGET)):+.2f}deg, "
      f"escaping {math.degrees(phase_error(t_esc, TARGET)):+.2f}deg\n")

# P1: add a -30deg error to the clean case.
phi_inject = -math.radians(30.0)
# P2: cancel the escaping case's error. Want w*t + phi == 0 (mod 2pi).
phi_cancel = -phase_error(t_esc, TARGET)
print(f"P1 phi = {phi_inject:+.6f} rad ({math.degrees(phi_inject):+.2f}deg)")
print(f"P2 phi = {phi_cancel:+.6f} rad ({math.degrees(phi_cancel):+.2f}deg)\n")

CASES = [
    ("P1_inject_phase", "settling_baseline_L7", 7.0, phi_inject),
    ("P2_cancel_phase", "dtheta_src_L7_th6.9",  6.9, phi_cancel),
]

for label, source_run, theta_prev, phi in CASES:
    src_dir = runs_dir / source_run
    t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
    cfg = {
        "motion": {"omega_b": omega_b_of(RPM), "theta_max": [TARGET, 0.0, 0.0],
                   "phi_angular": [phi, 0.0, 0.0]},
        "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
        "n_mix_cycles": CYCLES, "n_transition_cycles": CYCLES, "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(RPM)]},
        "initial_checkpoint": {
            "t_dump": t_dump, "omega_b": omega_b_of(RPM),
            "theta_max": [theta_prev, 0.0, 0.0],
            "phi_angular": [phi, 0.0, 0.0],
            "checkpoint_path": str(src_dir / "checkpoint.dump"),
        },
        "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
        "binary": "/oscar/scratch/eaguerov/BioReactor-mpi-testc",
        "submit": True, "exclude": "node2336",
    }
    print(f"{label}: {chain.submit_chain(cfg)}")
