"""Why do omega_b-changing restarts escape when theta ones no longer do?

T_per_st = T_per/T_bio really is independent of omega_b (U_bio ~ omega_b so
T_bio ~ 1/omega_b), so these carry no phase error and the v2 fix computes an
offset of ~0 for them. Yet both land on the same ~1.25 branch. Something
else is throwing them out of the reference basin.

Prints the phase offset the solver would compute, the trajectory shape (a
slow drift up means a growing mode as with the phase bug; starting high and
staying means a large initial kick), and the su rescale each case applies.

Usage:
    uv run python scripts/analyze_omega_escape.py
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


T_ND = t_per_nd(32.5, 7.0)
print(f"T_per_st at 17.5/32.5/37.5 rpm, theta=7: "
      + ", ".join(f"{t_per_nd(r, 7.0):.6f}" for r in (17.5, 32.5, 37.5))
      + "  (independent of rpm, as derived)\n")


def cyc(run, col=5, fn=np.mean):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, col]
    c = (t - t[0]) / T_ND
    return np.array([fn(y[(c >= k) & (c < k + 1)]) for k in range(int(c[-1]))])


ref = float(np.mean(cyc("kicktest_L7_th7")[-20:]))
CASES = [("rpm 17.5->32.5", "b231f1fc", 17.5), ("rpm 37.5->32.5", "58e2202b", 37.5),
         ("theta 2->7 (ok)", "2b90f675", 32.5)]

for label, run, src_rpm in CASES:
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t0 = float(d[0, 1])
    off = T_ND * round(t0 / T_ND) - t0
    p = cyc(run) / ref
    print(f"{label:18s} su={src_rpm/32.5:.3f}  solver phase offset "
          f"{360*off/T_ND:+.1f}deg")
    print("   cycle  " + "".join(f"{k:>7d}" for k in [0, 1, 2, 3, 5, 10, 20, 40, 80, 110]))
    print("   ratio  " + "".join(f"{p[k]:7.3f}" if k < len(p) else "      -"
                                 for k in [0, 1, 2, 3, 5, 10, 20, 40, 80, 110]))
