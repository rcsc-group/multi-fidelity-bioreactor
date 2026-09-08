"""Does the warm/cold agreement survive the choice of tau statistic?

scripts/postprocess.py (docstring at the tau KPI block) states the project's
Kim et al. Fig.13a analog explicitly:

    tau_liq_mean = max over time of the spatially-averaged tau
                 = max_t <tau(t)>_space, NOT a mean over time

so the headline KPI really is a per-period MAX of the spatial mean, which is
what the settling analysis has been using. postprocess also computes a
median-over-QSS-window variant (tau_*_qss). This checks the warm/cold
comparison against three statistics of the same per-cycle data:

    peak   max_t <tau>_space over the cycle   (Kim analog, the KPI)
    mean   time-average of <tau>_space        (a genuine time average)
    rms    root-mean-square over the cycle

If the agreement is an artifact of picking the max, the three will disagree.

Usage:
    uv run python scripts/analyze_kpi_robustness.py
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
STATS = {"peak": np.max, "mean": np.mean,
         "rms": lambda v: float(np.sqrt(np.mean(np.asarray(v) ** 2)))}


def per_cycle(run, fn):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, 5]
    c = (t - t[0]) / T_ND
    return np.array([fn(y[(c >= k) & (c < k + 1)]) for k in range(int(c[-1]))])


CASES = [("warm, dtheta=0", "ebd941af"), ("warm, dtheta=0.1", "89b3359a"),
         ("warm, dtheta=3", "2d527106"), ("warm, dtheta=5", "2b90f675"),
         ("chain seg0", "ef6f956d"), ("chain seg1", "e3fa8c4e"),
         ("warm, rpm 17.5->32.5", "b231f1fc"), ("warm, rpm 37.5->32.5", "58e2202b"),
         ("W-B omega path, su=0.997", "6ca01df0"),
         ("warm, rpm 30->32.5", "c3d1b1c4"),
         ("17.5->32.5 g mode1", "0f42fbb2"), ("17.5->32.5 g mode2", "54c714d3")]

print("late-window (last 20 cycles) value / cold start, by tau statistic\n")
print(f"{'case':18s}" + "".join(f"{n:>10s}" for n in STATS))
refs = {n: float(np.mean(per_cycle("kicktest_L7_th7", f)[-20:]))
        for n, f in STATS.items()}
rows = {n: [] for n in STATS}
for label, run in CASES:
    line = f"{label:18s}"
    for n, f in STATS.items():
        v = float(np.mean(per_cycle(run, f)[-20:])) / refs[n]
        rows[n].append(v)
        line += f"{v:10.4f}"
    print(line)
print()
for n in STATS:
    a = np.array(rows[n])
    sd = float(np.std(per_cycle("kicktest_L7_th7", STATS[n])[-20:]) / refs[n])
    print(f"  {n:5s}: max deviation {100*np.abs(a-1).max():.2f}%   "
          f"cold-start cycle-to-cycle noise {100*sd:.2f}%")
