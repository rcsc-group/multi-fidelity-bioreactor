"""L9, 22.5rpm: does the single-hop warm restart (su=0.778, from the
17.5rpm source) converge to the same value as an independent L9 cold
start, and does it save time?

Usage:
    uv run python scripts/analyze_l9_22p5_final.py
"""
import math
from pathlib import Path

import numpy as np

RUNS = Path(__file__).parent.parent / "runs"


def t_per_nd(rpm, theta=7.0):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    return T / (L / U)


T_ND = t_per_nd(22.5)


def cycle_series(run, col=5):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, col]
    c = (t - t[0]) / T_ND
    n = int(c[-1])
    peak = np.array([y[(c >= k) & (c < k + 1)].max() for k in range(n)])
    mean = np.array([y[(c >= k) & (c < k + 1)].mean() for k in range(n)])
    return peak, mean


cp, cm = cycle_series("l9_cold_rpm22.5")
wp, wm = cycle_series("l9_warm_rpm22.5")

ref_peak = float(np.mean(cp[-15:]))
ref_mean = float(np.mean(cm[-15:]))
warm_peak = float(np.mean(wp[-15:]))
warm_mean = float(np.mean(wm[-15:]))

print(f"cold start: {len(cp)} cycles, warm restart: {len(wp)} cycles\n")
print(f"tau_mean (peak-of-cycle KPI):  cold={ref_peak:.5e}  warm={warm_peak:.5e}"
      f"  ratio={warm_peak/ref_peak:.4f}")
print(f"tau_mean (time-avg, tighter):  cold={ref_mean:.5e}  warm={warm_mean:.5e}"
      f"  ratio={warm_mean/ref_mean:.4f}")

print(f"\ncycle-by-cycle (time-avg tau_mean / cold-start converged value):")
for k in list(range(0, 10)) + list(range(10, min(len(cp), len(wp)), 5)):
    if k < len(cp) and k < len(wp):
        print(f"  cyc {k:3d}: cold={cp[k]/ref_peak:.3f}  warm={wp[k]/ref_peak:.3f}")
