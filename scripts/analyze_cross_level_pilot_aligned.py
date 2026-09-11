"""Re-check analyze_cross_level_pilot.py's ~32% plateau with phase
alignment -- ruling out the exact artifact caught earlier this session
(anchoring intra-cycle phase on a run's first sample gives a constant,
fake offset that looks like non-convergence).

Usage:
    uv run python scripts/analyze_cross_level_pilot_aligned.py
"""
import math
from pathlib import Path

import numpy as np

RUNS = Path(__file__).parent.parent / "runs"
NPH = 256
PHASE = np.linspace(0, 1, NPH, endpoint=False)


def t_per_nd(rpm=32.5, theta=7.0):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    return T / (L / U)


T_ND = t_per_nd()


def waveforms(run, shift=0.0):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, 5]
    c = (t - t[0]) / T_ND
    out = []
    for k in range(int(c[-1])):
        m = (c >= k) & (c < k + 1)
        if m.sum() < 8:
            out.append(None)
            continue
        ph, v = c[m] - k, y[m]
        q = np.mod(PHASE + shift, 1.0)
        out.append(np.interp(q, ph, v, period=1.0))
    return out


def converged(ws, n=20):
    return np.mean([w for w in ws[-n:] if w is not None], axis=0)


cold = waveforms("kicktest_L7_th7")
conv = converged(cold)
scale = np.sqrt(np.mean(conv ** 2))

warm_base = waveforms("crosslevel_L6toL7_pilot")
cwarm = converged(warm_base)

best, bshift = None, 0.0
for s in np.linspace(-0.5, 0.5, 2001):
    q = np.mod(PHASE + s, 1.0)
    cand = np.interp(q, PHASE, cwarm, period=1.0)
    r = np.sqrt(np.mean((cand - conv) ** 2))
    if best is None or r < best:
        best, bshift = r, s

print(f"best-fit phase shift for the warm run's converged tail: "
      f"{bshift*360:.1f}deg (mod 180 since tau_mean is a |.|-mean: "
      f"{(bshift*360) % 180:.1f}deg)")

warm = waveforms("crosslevel_L6toL7_pilot", bshift)
err = np.array([np.sqrt(np.mean((w - conv) ** 2)) / scale if w is not None
               else np.nan for w in warm])
print(f"\nresidual after best-fit alignment, last 10 cycles: "
      f"{100*np.nanmean(err[-10:]):.2f}% (was ~32% unaligned)")
print("per-cycle (aligned):")
for k in range(0, len(err), 5):
    print(f"  cycle {k:3d}: {100*err[k]:.2f}%")
