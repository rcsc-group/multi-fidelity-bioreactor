"""Is the restart resuming at the wrong forcing phase?

The nondimensional time unit is T_bio = L_bio/U_bio, and U_bio depends on
theta_max through V_bio = L_bio/4*(H_bio + 0.5*L_bio*tan(theta_max)). So the
nondimensional rocking period T_per/T_bio is DIFFERENT for different
theta_max at the same omega_b.

A checkpoint is dumped at a zero-crossing of its own run, i.e. at a t that
is an integer number of ITS OWN nondimensional periods. Restart that dump
under a different theta_max and the same t is no longer an integer number
of the TARGET's nondimensional periods -- so the forcing resumes at some
nonzero phase while the restored velocity and interface fields correspond
to phase 0.

This prints, for every source checkpoint used in the sweep:
  - t_dump / T_nd(source theta): should be an integer (dump is at a
    zero-crossing of the source run)
  - t_dump / T_nd(target theta): the fractional part is the phase error the
    restart begins with
and checks it against the measured escape.

Usage:
    uv run python scripts/check_restart_phase.py
"""
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
RPM, TARGET = 32.5, 7.0


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(math.radians(theta)))
    return T_per / (L_bio / (V_bio / (H_bio * 0.5) / T_per))


T_TARGET = t_per_nd(RPM, TARGET)

# source run, its theta, and the restart run that used it (and its measured
# tau_mean ratio from analyze_branch_separation.py)
SOURCES = [
    ("settling_baseline_L7",        7.0, "restart dtheta=0 / B", 0.996),
    ("dtheta_src_L7_th6.9",         6.9, "restart dtheta=0.1 / C", 1.255),
    ("settling_ref_L7_rpm32.5_th4", 4.0, "restart dtheta=3", 1.257),
    ("settling_ref_L7_rpm32.5_th2", 2.0, "restart dtheta=5", 1.253),
]

print(f"target theta={TARGET:g}, T_nd(target) = {T_TARGET:.6f}\n")
print(f"{'source':30s} {'theta':>6s} {'t_dump':>10s} {'n_src':>9s} "
      f"{'n_target':>9s} {'phase err':>10s} {'measured':>9s}")
for run, theta, label, measured in SOURCES:
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t_dump = float(d[-1, 1])
    n_src = t_dump / t_per_nd(RPM, theta)
    n_tgt = t_dump / T_TARGET
    frac = n_tgt - round(n_tgt)
    print(f"{run:30s} {theta:6.2f} {t_dump:10.5f} {n_src:9.4f} "
          f"{n_tgt:9.4f} {frac * 360:9.1f}d {measured:9.3f}")

print("\nn_src is the dump's own period count -- integer means the dump really")
print("is at a zero-crossing of the run that wrote it. phase err is what the")
print("restart begins with, in degrees of the target's rocking cycle.")
