"""Warm vs cold settling, measured on the tau_mean WAVEFORM.

The per-cycle-peak metric used until now is invalid for this question: the
peak tracks the smooth-step ramp amplitude almost exactly (cycle 0/1/2 peaks
0.303/0.697/0.986 against alpha 0.259/0.741/1.000), so it reports settling
as soon as N_RAMP_CYCLES elapses regardless of the real transient, and its
tightest tolerances are dominated by ~0.4% inter-run offsets.

This compares each cycle's full tau_mean waveform, resampled onto a common
intra-cycle phase grid, against the cold start's converged waveform. Warm
and cold starts are phase-comparable because the phase fix puts the restart
instant at forcing phase 0, which is also where the cold start begins.

Usage:
    uv run python scripts/analyze_settling_waveform.py
"""
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
NPH = 64
PHASE = np.linspace(0, 1, NPH, endpoint=False)


def t_per_nd(rpm, theta):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T / (L / (V / (H * 0.5) / T))


T_ND = t_per_nd(32.5, 7.0)


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
conv = np.mean([x for x in cold[-20:] if x is not None], axis=0)
scale = np.sqrt(np.mean(conv ** 2))


def err(w):
    return np.array([np.sqrt(np.mean((x - conv) ** 2)) / scale if x is not None
                     else np.nan for x in w])


def settle(e, tol):
    return next((k for k in range(len(e))
                 if np.all(e[k:] < tol)), None)


CASES = [("cold start", "kicktest_L7_th7"), ("warm, dtheta=0", "ebd941af"),
         ("warm, dtheta=0.1", "89b3359a"), ("warm, dtheta=3", "2d527106"),
         ("warm, dtheta=5", "2b90f675"), ("chain seg1", "e3fa8c4e")]

floor = float(np.nanmean(err(cold)[-20:]))
print(f"waveform RMS error vs the cold start's converged cycle; "
      f"noise floor {100*floor:.2f}%\n")
print(f"{'case':18s}" + "".join(f"{'cyc' + str(k):>8s}" for k in range(8))
      + f"{'@5%':>6s}{'@2%':>6s}{'@1%':>6s}")
for label, run in CASES:
    e = err(waveforms(run))
    row = f"{label:18s}" + "".join(
        f"{100*e[k]:7.1f}%" if k < len(e) else f"{'':>8s}" for k in range(8))
    for tol in (0.05, 0.02, 0.01):
        s = settle(e, tol)
        row += f"{('--' if s is None else s):>6}"
    print(row)
