"""How far apart are the two "branches", in nondimensional terms?

All L7, 32.5rpm, target theta_max=7deg, so every run below shares the same
U_bio and the raw shear_stress.dat columns (already Basilisk-nondimensional)
are directly comparable without any rescaling.

Reference branch = kicktest_L7_th7: a single continuous cold start, 140
cycles, whose one-shot 1e-4 velocity kick at cycle 20 provably did nothing.

Usage:
    uv run python scripts/analyze_branch_separation.py
"""
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"

RPM, THETA = 32.5, 7.0


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(math.radians(theta)))
    U_bio = V_bio / (H_bio * 0.5) / T_per
    return T_per / (L_bio / U_bio)


T_ND = t_per_nd(RPM, THETA)


def per_cycle_peak(path):
    d = np.loadtxt(path, skiprows=1)
    t, tau = d[:, 1], d[:, 5]
    t0 = t[0]
    n = int((t[-1] - t0) / T_ND)
    out = np.full(n, np.nan)
    cyc = (t - t0) / T_ND
    for c in range(n):
        m = (cyc >= c) & (cyc < c + 1)
        if m.sum():
            out[c] = tau[m].max()
    return out


CASES = [
    ("cold start (reference)", "kicktest_L7_th7", None),
    ("restart dtheta=0",       "8fd81f04",        0.0),
    ("restart dtheta=0.1",     "c2001da9",        0.1),
    ("restart dtheta=0.1 slow ramp", "f1ec9a45",  0.1),
    ("restart dtheta=3",       "1ddbba85",        3.0),
    ("restart dtheta=5",       "3324f34b",        5.0),
    ("B ramp on, state on-target",  "ba2f3cc1",   1e-7),
    ("C ramp off, state off-target", "46196ea2",  0.0),
    ("P1 clean case + 30deg phase error", "35bde376", 0.0),
    ("P2 dtheta=0.1, phase corrected",    "57eeee1a", 0.1),
    ("FIX dtheta=0",   "2e69509a", 0.0),
    ("FIX dtheta=0.1", "3537850c", 0.1),
    ("FIX dtheta=3",   "8108abc0", 3.0),
    ("FIX dtheta=5",   "fa718f3a", 5.0),
    ("FIX dtheta=0.1, LEGACY su", "c1437263", 0.1),
    ("FIX dtheta=5,   LEGACY su", "277f6876", 5.0),
    ("v1 CHAIN seg0 (th2->7)", "5f931126", 5.0),
    ("v1 CHAIN seg1 (th7->7)", "475a8839", 0.0),
    ("v2 dtheta=0",   "ebd941af", 0.0),
    ("v2 dtheta=0.1", "89b3359a", 0.1),
    ("v2 dtheta=3",   "2d527106", 3.0),
    ("v2 dtheta=5",   "2b90f675", 5.0),
    ("v2 CHAIN seg0 (th2->7)", "ef6f956d", 5.0),
    ("v2 CHAIN seg1 (th7->7)", "e3fa8c4e", 0.0),
]

ref_peaks = per_cycle_peak(RUNS / "kicktest_L7_th7" / "shear_stress.dat")
ref = float(np.nanmean(ref_peaks[-20:]))
ref_sd = float(np.nanstd(ref_peaks[-20:]))
print(f"reference branch (cold start, last 20 of {len(ref_peaks)} cycles):")
print(f"  tau_nd = {ref:.5e}  +/- {ref_sd:.2e} (cycle-to-cycle sd, {100*ref_sd/ref:.2f}%)\n")

print(f"{'case':32s} {'ncyc':>5s} {'tau_nd':>11s} {'tau/tau_ref':>12s} {'sd%':>6s}")
for label, run, dth in CASES:
    p = per_cycle_peak(RUNS / run / "shear_stress.dat")
    m = float(np.nanmean(p[-20:]))
    s = float(np.nanstd(p[-20:]))
    print(f"{label:32s} {len(p):5d} {m:11.5e} {m/ref:12.4f} {100*s/m:6.2f}")
