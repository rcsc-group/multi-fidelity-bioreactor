"""Does a generic 0.3% perturbation of a pure cold start escape?

kicknoise_L7_th7: continuous cold start, no checkpoint anywhere, one-shot
per-cell u *= (1 + 3e-3*noise()) at cycle 40 (confirmed in the .err log:
"kick: mode 1, relative amplitude 0.003, at t=24.3"), 200 cycles. 3e-3
matches the ~0.3% state mismatch that does escape on the restart path.

  escapes -> the reference state is linearly stable with a finite basin,
      0.3% is outside it, and nothing about the escape is restart-specific
  decays  -> 0.3% generic is inside the basin, so the restart path is
      injecting something more than a state offset

Usage:
    uv run python scripts/analyze_kick_noise.py
"""
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
KICK_CYCLE = 40


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(math.radians(theta)))
    return T_per / (L_bio / (V_bio / (H_bio * 0.5) / T_per))


T_ND = t_per_nd(32.5, 7.0)


def cycle_peaks(run, col, fname="shear_stress.dat"):
    d = np.loadtxt(RUNS / run / fname, skiprows=1)
    t, y = d[:, 1], d[:, col]
    cyc = (t - t[0]) / T_ND
    return np.array([y[(cyc >= c) & (cyc < c + 1)].max() for c in range(int(cyc[-1]))])


ref = float(np.mean(cycle_peaks("kicktest_L7_th7", 5)[-20:]))
p_mean = cycle_peaks("kicknoise_L7_th7", 5)
p_max = cycle_peaks("kicknoise_L7_th7", 4)
ref100 = float(np.mean(cycle_peaks("kicktest_L7_th7", 4)[-20:]))

print(f"kicknoise: {len(p_mean)} cycles, kick at cycle {KICK_CYCLE}\n")
print(f"  pre-kick  (cycles 20-39) tau_mean/ref = {np.mean(p_mean[20:40]) / ref:.4f}")
print(f"{'cycles':>12s} {'tau_mean/ref':>13s} {'tau_100/ref':>12s}")
for a in range(40, len(p_mean) - 19, 20):
    print(f"{a:5d}-{a + 19:<6d} {np.mean(p_mean[a:a + 20]) / ref:13.4f}"
          f" {np.mean(p_max[a:a + 20]) / ref100:12.4f}")

# The escape's fingerprint, for comparison with the restart case
# (uy_liq_max 4.30, posY_max 1.65x, Omega_liq_rms 1.03, area 0.965).
print("\nlate-window fingerprint (last 20 cycles), kicknoise / cold-start reference")
for fname, cols in [("normf.dat", {3: "Omega_liq_rms", 5: "Omega_liq_max",
                                   11: "uy_liq_rms", 13: "uy_liq_max"}),
                    ("vol_frac_interf.dat", {3: "f_liq_interf", 4: "posY_max"})]:
    for col, name in cols.items():
        a = cycle_peaks("kicktest_L7_th7", col, fname)
        b = cycle_peaks("kicknoise_L7_th7", col, fname)
        print(f"  {name:16s} {np.mean(b[-20:]) / np.mean(a[-20:]):8.4f}")
