"""Does the escaped branch saturate, or keep growing?

Run d6e8428e: case C (ramp off, state off target -- i.e. a restart carrying
a -23.9deg forcing phase error) taken out to 400 cycles instead of 120. The
120-cycle runs were all still climbing at the end, so "a second attractor"
was never established. Now that the phase error is identified as the cause,
this says what the spurious response actually does over long times: a
saturating offset is a different (wrong) periodic state, monotone growth is
a numerical instability being fed continuously.

Usage:
    uv run python scripts/analyze_longrun_escape.py
"""
import math
from pathlib import Path

import numpy as np

R = Path(__file__).parent.parent / "runs"


def t_per_nd(rpm, theta):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T / (L / (V / (H * 0.5) / T))


T_ND = t_per_nd(32.5, 7.0)


def peaks(run, col):
    d = np.loadtxt(R / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, col]
    c = (t - t[0]) / T_ND
    return np.array([y[(c >= k) & (c < k + 1)].max() for k in range(int(c[-1]))])


ref = float(np.mean(peaks("kicktest_L7_th7", 5)[-20:]))
ref100 = float(np.mean(peaks("kicktest_L7_th7", 4)[-20:]))
p, q = peaks("d6e8428e", 5), peaks("d6e8428e", 4)
print(f"escaped run taken to {len(p)} cycles, 40-cycle block means vs cold start\n")
print(f"{'block':>13s} {'tau_mean/ref':>13s} {'tau_100/ref':>12s}")
for a in range(0, len(p) - 39, 40):
    print(f"{a:4d}-{a + 39:<8d} {np.mean(p[a:a + 40]) / ref:13.4f}"
          f" {np.mean(q[a:a + 40]) / ref100:12.4f}")
