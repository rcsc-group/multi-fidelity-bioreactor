"""Does the escaped omega branch survive refinement to L8?

At L7, omega 17.5 -> 32.5 (su=0.538) lands at tau/tau_ref = 1.2308 -- the
escaped branch. This repeats the identical restart at L8 and compares
against settling_baseline_L8, an independent L8 cold start at the same
condition (32.5rpm/theta=7).

Usage:
    uv run python scripts/analyze_omega_L8.py
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


def late_stats(run, n_last=20):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, 5]
    c = (t - t[0]) / T_ND
    n = int(c[-1])
    peaks = np.array([y[(c >= k) & (c < k + 1)].max() for k in range(n)])
    means = np.array([y[(c >= k) & (c < k + 1)].mean() for k in range(n)])
    return n, peaks, means


n_ref, pk_ref, mn_ref = late_stats("settling_baseline_L8")
n_esc, pk_esc, mn_esc = late_stats("c075d658")

print(f"L8 cold-start reference: {n_ref} cycles")
print(f"L8 omega 17.5->32.5 restart: {n_esc} cycles\n")

pk_r = float(np.mean(pk_ref[-20:]))
mn_r = float(np.mean(mn_ref[-20:]))
n_last = min(20, n_esc)
pk_e = float(np.mean(pk_esc[-n_last:]))
mn_e = float(np.mean(mn_esc[-n_last:]))

print(f"peak(tau_mean)/ref:  {pk_e/pk_r:.4f}   (last {n_last} of {n_esc} cycles)")
print(f"mean(tau_mean)/ref:  {mn_e/mn_r:.4f}")
print(f"\nper-cycle peak ratio trajectory (escaped run / L8 ref late mean):")
for k in list(range(0, min(10, n_esc))) + list(range(10, n_esc, 5)):
    print(f"  cycle {k:3d}: {pk_esc[k]/pk_r:.4f}")

print(f"\nfor comparison, L7 result at the same condition: 1.2308 (peak-based)")
