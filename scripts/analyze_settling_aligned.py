"""Warm vs cold settling, with analysis-side phase misalignment removed.

analyze_settling_waveform.py used each run's first sample as intra-cycle
phase zero. That is not exactly the forcing zero-crossing, so every cycle of
a run picks up the SAME constant phase error -- which shows up as a flat,
non-decaying residual (5.2% for the bit-exact dtheta=0 restart, which by
construction has no transient at all). That is an artifact of the analysis,
not of the runs.

Here each run's converged waveform is cross-correlated against the cold
start's to find the best sub-sample phase shift; that shift is applied to
all of the run's cycles before measuring the residual. The shift is reported
too, because after the phase fix it should be small -- it is a direct
measurement of residual phase error, independent of the solver.

Usage:
    uv run python scripts/analyze_settling_aligned.py
"""
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
NPH = 256
PHASE = np.linspace(0, 1, NPH, endpoint=False)


def t_per_nd(rpm, theta):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T / (L / (V / (H * 0.5) / T))


T_ND = t_per_nd(32.5, 7.0)


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
        # periodic resample at PHASE + shift
        q = np.mod(PHASE + shift, 1.0)
        out.append(np.interp(q, ph, v, period=1.0))
    return out


def converged(w, n=20):
    return np.mean([x for x in w[-n:] if x is not None], axis=0)


cold = waveforms("kicktest_L7_th7")
conv = converged(cold)
scale = np.sqrt(np.mean(conv ** 2))

CASES = [("cold start", "kicktest_L7_th7"), ("warm, dtheta=0", "ebd941af"),
         ("warm, dtheta=0.1", "89b3359a"), ("warm, dtheta=3", "2d527106"),
         ("warm, dtheta=5", "2b90f675"), ("chain seg1", "e3fa8c4e")]

print(f"residual after optimal phase alignment; cold start's own noise "
      f"floor sets the scale\n")
print(f"{'case':18s}{'phase shift':>12s}" +
      "".join(f"{'cyc' + str(k):>8s}" for k in range(8)) +
      f"{'@5%':>5s}{'@2%':>5s}{'@1%':>5s}")
for label, run in CASES:
    # scan shifts; pick the one minimising the converged-waveform residual
    grid = np.linspace(-0.5, 0.5, 2001)
    base = waveforms(run)
    cbase = converged(base)
    best, bshift = None, 0.0
    for s in grid:
        q = np.mod(PHASE + s, 1.0)
        cand = np.interp(q, PHASE, cbase, period=1.0)
        r = np.sqrt(np.mean((cand - conv) ** 2))
        if best is None or r < best:
            best, bshift = r, s
    w = waveforms(run, bshift)
    e = np.array([np.sqrt(np.mean((x - conv) ** 2)) / scale
                  if x is not None else np.nan for x in w])
    row = f"{label:18s}{bshift * 360:10.1f}d " + "".join(
        f"{100 * e[k]:7.1f}%" if k < len(e) else f"{'':>8s}" for k in range(8))
    for tol in (0.05, 0.02, 0.01):
        s = next((k for k in range(len(e)) if np.all(e[k:] < tol)), None)
        row += f"{('--' if s is None else s):>5}"
    print(row)
