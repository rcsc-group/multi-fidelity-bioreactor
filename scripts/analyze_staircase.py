"""Does stepping omega_b in sub-threshold hops land on the reference branch?

17.5 -> 20 -> 23 -> 26.5 -> 30 -> 32.5 rpm, each hop su ~ 0.87-0.92, all
below the escaping threshold measured at |su-1| >= 0.154. Only the final
segment (at 32.5rpm/theta=7/L7) is compared against the cold-start
reference; intermediate segments just have to stay converged on their own
branch, which is unverifiable without their own cold-start references, so
they're reported but not judged.

Usage:
    uv run python scripts/analyze_staircase.py
"""
import math
from pathlib import Path

import numpy as np

RUNS = Path(__file__).parent.parent / "runs"


def t_per_nd(rpm, theta):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T / (L / (V / (H * 0.5) / T))


SEGS = [("17.5->20", "b2f44191", 20.0), ("20->23", "e9d01d94", 23.0),
        ("23->26.5", "f48b804a", 26.5), ("26.5->30", "677f5bc5", 30.0),
        ("30->32.5", "5c877b33", 32.5)]

ref_run = "kicktest_L7_th7"
T_ND_ref = t_per_nd(32.5, 7.0)
d = np.loadtxt(RUNS / ref_run / "shear_stress.dat", skiprows=1)
t, y = d[:, 1], d[:, 5]
c = (t - t[0]) / T_ND_ref
n = int(c[-1])
peaks = np.array([y[(c >= k) & (c < k + 1)].mean() for k in range(n)])
ref = float(np.mean(peaks[-20:]))

print(f"cold-start reference (32.5rpm/theta=7): {ref:.5e}\n")
for label, run, rpm in SEGS:
    T_ND = t_per_nd(rpm, 7.0)
    dd = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    tt, yy = dd[:, 1], dd[:, 5]
    cc = (tt - tt[0]) / T_ND
    nn = int(cc[-1])
    p = np.array([yy[(cc >= k) & (cc < k + 1)].mean() for k in range(nn)])
    late = float(np.mean(p[-20:])) if nn >= 20 else float(np.mean(p[-min(nn,5):]))
    tag = f"  ratio/32.5rpm-ref = {late/ref:.4f}" if rpm == 32.5 else ""
    print(f"{label:12s} ({rpm:g}rpm, {nn} cycles): late-window mean = {late:.5e}{tag}")
