"""Is "the cold start settles in 2 cycles" real, or an artifact of the metric?

N_RAMP_CYCLES=3, so the forcing amplitude is still ramping through cycle 2:
the smooth-step alpha = 3x^2-2x^3 over x = elapsed/(3*T_per_st) reaches 1
only at the END of cycle 2. The settling metric takes the per-cycle MAXIMUM
of tau_mean, so a cycle that finishes at full amplitude registers close to
full amplitude -- the metric cannot tell "ramp just completed" from
"converged", and would report settling at the ramp's end no matter how long
the physical transient actually is.

This prints, per cycle: alpha at the cycle's start and end, the per-cycle
peak, and -- as a metric that does NOT share that blind spot -- the RMS
distance between the cycle's full tau_mean waveform and the converged
waveform, resampled onto a common phase grid. A peak can match while the
waveform shape is still wrong; this catches that.

Usage:
    uv run python scripts/analyze_ramp_vs_settling.py
"""
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
N_RAMP = 3


def t_per_nd(rpm, theta):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T / (L / (V / (H * 0.5) / T))


T_ND = t_per_nd(32.5, 7.0)
NPH = 64
PHASE = np.linspace(0, 1, NPH, endpoint=False)


def alpha_at(cycle_frac):
    x = min(cycle_frac / N_RAMP, 1.0)
    return 3 * x * x - 2 * x * x * x


def load(run):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    return d[:, 1], d[:, 5]


def waveforms(run):
    """Per-cycle tau_mean resampled onto a common intra-cycle phase grid."""
    t, y = load(run)
    c = (t - t[0]) / T_ND
    out = []
    for k in range(int(c[-1])):
        m = (c >= k) & (c < k + 1)
        if m.sum() < 8:
            out.append(None)
            continue
        out.append(np.interp(PHASE, c[m] - k, y[m]))
    return out


w = waveforms("kicktest_L7_th7")
conv = np.mean([x for x in w[-20:] if x is not None], axis=0)
scale = np.sqrt(np.mean(conv ** 2))
peaks = np.array([x.max() if x is not None else np.nan for x in w])
pconv = float(np.nanmean(peaks[-20:]))

print(f"cold start, L7 32.5rpm theta=7.  ramp = {N_RAMP} cycles\n")
print(f"{'cycle':>5s} {'alpha start':>12s} {'alpha end':>10s} "
      f"{'peak/conv':>10s} {'waveform RMS err':>17s}")
for k in list(range(10)) + [12, 15, 20, 30, 40, 60, 80, 100]:
    if k >= len(w) or w[k] is None:
        continue
    err = np.sqrt(np.mean((w[k] - conv) ** 2)) / scale
    print(f"{k:5d} {alpha_at(k):12.3f} {alpha_at(k + 1):10.3f} "
          f"{peaks[k] / pconv:10.4f} {100 * err:16.2f}%")

for tol in (0.05, 0.02, 0.01):
    s = next((k for k in range(len(w))
              if all(w[j] is not None
                     and np.sqrt(np.mean((w[j] - conv) ** 2)) / scale < tol
                     for j in range(k, len(w)))), None)
    print(f"waveform settles @{int(tol*100)}%: cycle {s}")
