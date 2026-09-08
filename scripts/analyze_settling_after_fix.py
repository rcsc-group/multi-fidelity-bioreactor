"""Post-fix: do warm and cold starts converge to the same mean, and does
warm-starting converge FASTER?

The second question is the settling-time study's original one, and every
earlier measurement of it was contaminated by the restart forcing-phase
bug, so it has to be asked again from scratch.

Settling cycle = the first cycle from which the per-cycle-peak tau_mean
stays within tol of the target for the REST of the run (not merely touches
it once). Target is the cold start's own converged value, so cold and warm
starts are judged against the same number.

Usage:
    uv run python scripts/analyze_settling_after_fix.py
"""
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"


def t_per_nd(rpm, theta):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T / (L / (V / (H * 0.5) / T))


T_ND = t_per_nd(32.5, 7.0)


def peaks(run):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, 5]
    c = (t - t[0]) / T_ND
    return np.array([y[(c >= k) & (c < k + 1)].max() for k in range(int(c[-1]))])


def settle(p, target, tol):
    rel = np.abs(p - target) / target
    for s in range(len(rel)):
        if np.all(rel[s:] < tol):
            return s
    return None


cold = peaks("kicktest_L7_th7")
TARGET = float(np.mean(cold[-20:]))
NOISE = float(np.std(cold[-20:])) / TARGET

CASES = [("cold start", "kicktest_L7_th7"), ("warm, dtheta=0", "ebd941af"),
         ("warm, dtheta=0.1", "89b3359a"), ("warm, dtheta=3", "2d527106"),
         ("warm, dtheta=5", "2b90f675"), ("chain seg0", "ef6f956d"),
         ("chain seg1", "e3fa8c4e")]

print(f"target = cold start's converged value; its own cycle-to-cycle "
      f"noise is {100 * NOISE:.2f}%\n")
print(f"{'case':20s} {'ncyc':>5s} {'mean/target':>12s} "
      + "".join(f"{'settle@' + t:>12s}" for t in ("5%", "2%", "1%")))
means = []
for label, run in CASES:
    p = peaks(run)
    m = float(np.mean(p[-20:])) / TARGET
    if label != "cold start":
        means.append(m)
    row = f"{label:20s} {len(p):5d} {m:12.4f}"
    for tol in (0.05, 0.02, 0.01):
        s = settle(p, TARGET, tol)
        row += f"{('never' if s is None else s):>12}"
    print(row)

a = np.array(means)
print(f"\nfour warm starts: mean {a.mean():.4f}, spread {a.max() - a.min():.4f}, "
      f"max deviation from cold start {np.abs(a - 1).max() * 100:.2f}%")
print(f"cold start's own 20-cycle sampling error on the target: "
      f"{100 * NOISE / math.sqrt(20):.2f}%")
