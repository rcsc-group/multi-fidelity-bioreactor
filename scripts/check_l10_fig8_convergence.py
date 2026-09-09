"""Is the L10 tau signal actually periodic/settled, checked against the
real physics rather than against n_mix_cycles (a mixing-time-study
budget, unrelated to how fast the FLOW itself settles)?

Stitches together the two chained segments that carry this L10 flow's
history as far back as it's preserved on disk (l10_kim_seg2, t=10.34-
20.64, then l10_kim_fig8_signed, t=20.66-22.46 -- whatever came before
t=10.34, the run's true t=0 origin, isn't in runs/ under any name found).
That's ~20 continuous cycles, not the 2 a single file suggested.

Usage:
    uv run python scripts/check_l10_fig8_convergence.py
"""
import math
from pathlib import Path

import numpy as np

RUNS = Path(__file__).parent.parent / "runs"


def t_per_nd(rpm=32.5, theta=7.0):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    return T / (L / U)


T_ND = t_per_nd()

d1 = np.loadtxt(RUNS / "l10_kim_seg2" / "shear_stress.dat", skiprows=1)
d2 = np.loadtxt(RUNS / "l10_kim_fig8_signed" / "shear_stress.dat", skiprows=1)
# [BUG, CAUGHT] first version of this script used column 10
# (tau_mean_signed) for d2 but column 5 (tau_mean, unsigned) for d1 --
# different quantities (tau_mean_signed cancels sign, so it's roughly
# HALF of tau_mean by construction), producing a fake ~2x "drop" at the
# segment boundary that looked like a real transient. Both files actually
# carry the same "tau_mean" column at index 5 -- use that, consistently.
t = np.concatenate([d1[:, 1], d2[:, 1]])
tau = np.concatenate([d1[:, 5], d2[:, 5]])
order = np.argsort(t)
t, tau = t[order], tau[order]

cyc = t / T_ND  # ABSOLUTE cycle index (t=0 = the flow's true origin, not this file's start)
n0, n1 = int(cyc[0]), int(cyc[-1])
print(f"combined coverage: absolute cycle {n0} to {n1} "
      f"({n1 - n0} cycles, t={t[0]:.2f} to {t[-1]:.2f})\n")
print(f"{'cycle':>6s} {'max |tau_mean|':>16s}")
for k in range(n0, n1 + 1):
    m = (cyc >= k) & (cyc < k + 1)
    if m.sum():
        print(f"{k:6d} {tau[m].max():16.5e}")
