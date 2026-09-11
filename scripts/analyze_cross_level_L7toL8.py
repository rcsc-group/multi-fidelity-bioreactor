"""Does the L6->L7 cross-level discrepancy persist at L7->L8, or vanish
like the omega-checkpoint escape did going L7->L8?

Reference: settling_baseline_L8 (32.5rpm/theta=7, 80-cycle cold start --
zero new compute).

Usage:
    uv run python scripts/analyze_cross_level_L7toL8.py
"""
import math
from pathlib import Path

import numpy as np

RUNS = Path(__file__).parent.parent / "runs"
NPH = 64
PHASE = np.linspace(0, 1, NPH, endpoint=False)


def t_per_nd(rpm=32.5, theta=7.0):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    return T / (L / U)


T_ND = t_per_nd()


def waveforms(run):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, 5]
    c = (t - t[0]) / T_ND
    out = []
    for k in range(int(c[-1])):
        m = (c >= k) & (c < k + 1)
        out.append(np.interp(PHASE, c[m] - k, y[m]) if m.sum() >= 8 else None)
    return out


cold = waveforms("settling_baseline_L8")
warm = waveforms("crosslevel_L7toL8_pilot")
conv = np.mean([w for w in cold[-20:] if w is not None], axis=0)
scale = np.sqrt(np.mean(conv ** 2))


def err(ws):
    return np.array([np.sqrt(np.mean((w - conv) ** 2)) / scale if w is not None
                     else np.nan for w in ws])


ec, ew = err(cold), err(warm)
print(f"{'cycle':>6s} {'L8 cold start':>14s} {'L7->L8 warm':>14s}")
for k in range(0, min(len(ec), len(ew)), 5):
    print(f"{k:6d} {100*ec[k]:13.2f}% {100*ew[k]:13.2f}%")

print(f"\nL8 cold start, last 10 cycles:  {100*np.nanmean(ec[-10:]):.2f}%")
print(f"L7->L8 warm,   last 10 cycles:  {100*np.nanmean(ew[-10:]):.2f}%")
print(f"\n(for comparison, L6->L7 pilot's last-10-cycle residual was ~32%)")
