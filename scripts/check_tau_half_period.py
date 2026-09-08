"""Is tau_mean half-period symmetric?

The alignment scan returned shifts clustered near 170-175deg rather than 0,
which would be a large residual phase error if the tau_mean waveform had the
rocking period. It should not: tau_mean is a magnitude-like quantity (a
volume mean of |tau|), so it completes a full cycle on each half-stroke and
its period is T_per/2. If so, 0deg and 180deg are degenerate minima and the
scan's choice between them is meaningless -- the shifts must be read modulo
180deg, leaving 5-10deg, which is about one output sample (dt_out ~ 0.02 vs
T_per_st 0.607 -> 11.9deg) and therefore an artifact of anchoring analysis
phase on the run's first sample.

Usage:
    uv run python scripts/check_tau_half_period.py
"""
import math
from pathlib import Path

import numpy as np

RUNS = Path(__file__).parent.parent / "runs"
NPH = 256
PHASE = np.linspace(0, 1, NPH, endpoint=False)


def t_per_nd(rpm, theta):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T / (L / (V / (H * 0.5) / T))


T_ND = t_per_nd(32.5, 7.0)
d = np.loadtxt(RUNS / "kicktest_L7_th7" / "shear_stress.dat", skiprows=1)
t, y = d[:, 1], d[:, 5]
c = (t - t[0]) / T_ND
w = [np.interp(PHASE, c[m] - k, y[m]) for k in range(int(c[-1]))
     if (m := ((c >= k) & (c < k + 1))).sum() >= 8]
conv = np.mean(w[-20:], axis=0)
scale = np.sqrt(np.mean(conv ** 2))

half = np.interp(np.mod(PHASE + 0.5, 1.0), PHASE, conv, period=1.0)
print(f"converged cold-start tau_mean waveform:")
print(f"  RMS( w(phi+180deg) - w(phi) ) / RMS(w) = "
      f"{100*np.sqrt(np.mean((half-conv)**2))/scale:.2f}%")
print(f"  cycle-to-cycle noise floor              = "
      f"{100*np.mean([np.sqrt(np.mean((x-conv)**2))/scale for x in w[-20:]]):.2f}%")
print(f"\none output sample as a phase = "
      f"{360*np.median(np.diff(t))/T_ND:.1f}deg")
