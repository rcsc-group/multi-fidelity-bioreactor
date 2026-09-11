"""Does L6->L7 spatial-interpolation warm-starting settle faster than an
ordinary L7 cold start?

Reference: kicktest_L7_th7 (32.5rpm/theta=7, 139-cycle cold start, already
validated and used as the reference throughout this session's phase-bug
investigation -- zero new compute).

Usage:
    uv run python scripts/analyze_cross_level_pilot.py
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
NPH = 64
PHASE = np.linspace(0, 1, NPH, endpoint=False)


def waveforms(run):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, 5]
    c = (t - t[0]) / T_ND
    out = []
    for k in range(int(c[-1])):
        m = (c >= k) & (c < k + 1)
        out.append(np.interp(PHASE, c[m] - k, y[m]) if m.sum() >= 8 else None)
    return out


cold = waveforms("kicktest_L7_th7")
warm = waveforms("crosslevel_L6toL7_pilot")
conv = np.mean([w for w in cold[-20:] if w is not None], axis=0)
scale = np.sqrt(np.mean(conv ** 2))


def err(ws):
    return np.array([np.sqrt(np.mean((w - conv) ** 2)) / scale if w is not None
                     else np.nan for w in ws])


ec, ew = err(cold), err(warm)


def settle(e, tol):
    return next((k for k in range(len(e)) if np.all(e[k:] < tol)), None)


print(f"waveform RMS error vs the converged cold-start cycle\n")
print(f"{'cycle':>6s} {'cold start':>12s} {'L6->L7 warm':>12s}")
for k in range(0, 15):
    cs = f"{100*ec[k]:.2f}%" if k < len(ec) else "-"
    ws = f"{100*ew[k]:.2f}%" if k < len(ew) else "-"
    print(f"{k:6d} {cs:>12s} {ws:>12s}")

print()
for tol in (0.05, 0.02, 0.01):
    print(f"settles @{int(100*tol)}%:  cold={settle(ec,tol)}  warm(L6->L7)={settle(ew,tol)}")
