"""When does the escaped branch's 3.8x shear maximum appear, and is the
velocity field itself different or only the stress diagnostic?

normf.dat carries independent liquid-phase diagnostics (vorticity and
velocity: avg / rms / volume / max), computed by a completely separate
reduction from the tau columns in shear_stress.dat.

  - if Omega_liq_max jumps ~4x alongside tau_100 while Omega_liq_rms rises
    only ~25% like tau_mean, the escaped branch carries a real localized
    vorticity spike on top of a globally more energetic flow
  - if the velocity norms are flat and only tau moves, the stress
    diagnostic is the thing that differs, not the flow

Usage:
    uv run python scripts/analyze_branch_onset.py
"""
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"

REF, ESC = "kicktest_L7_th7", "c2001da9"   # cold start vs dtheta=0.1 restart


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(math.radians(theta)))
    return T_per / (L_bio / (V_bio / (H_bio * 0.5) / T_per))


T_ND = t_per_nd(32.5, 7.0)


def cycle_peaks(path, col):
    d = np.loadtxt(path, skiprows=1)
    t, y = d[:, 1], d[:, col]
    t0 = t[0]
    cyc = (t - t0) / T_ND
    n = int(cyc[-1])
    return np.array([y[(cyc >= c) & (cyc < c + 1)].max() for c in range(n)])


# --- 1. onset: cycle-by-cycle tau_100 and tau_mean, escaped / reference
print("cycle-by-cycle ratio to the cold-start reference (same cycle index)\n")
print(f"{'cycle':>6s} {'tau_100':>9s} {'tau_mean':>9s}")
r100, rmean = cycle_peaks(RUNS / REF / "shear_stress.dat", 4), cycle_peaks(RUNS / REF / "shear_stress.dat", 5)
e100, emean = cycle_peaks(RUNS / ESC / "shear_stress.dat", 4), cycle_peaks(RUNS / ESC / "shear_stress.dat", 5)
n = min(len(r100), len(e100))
for c in list(range(0, 12)) + list(range(12, n, 8)):
    if c < n:
        print(f"{c:6d} {e100[c]/r100[c]:9.3f} {emean[c]/rmean[c]:9.3f}")

# --- 2. independent velocity/vorticity diagnostics
NORMF = {2: "Omega_liq_avg", 3: "Omega_liq_rms", 5: "Omega_liq_max",
         7: "ux_liq_rms", 9: "ux_liq_max", 11: "uy_liq_rms", 13: "uy_liq_max"}
print("\nnormf.dat, per-cycle peak, mean of last 20 cycles, escaped / reference")
for col, name in NORMF.items():
    try:
        a = cycle_peaks(RUNS / REF / "normf.dat", col)
        b = cycle_peaks(RUNS / ESC / "normf.dat", col)
    except (IndexError, OSError):
        continue
    print(f"  {name:16s} {np.mean(b[-20:]) / np.mean(a[-20:]):8.4f}")
