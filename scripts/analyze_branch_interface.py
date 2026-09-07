"""Is the escaped branch a genuinely bigger sloshing wave, or parasitic
interfacial currents -- and has it converged at all?

analyze_branch_onset.py: bulk vorticity is the same to 3% (Omega_liq_avg
1.028, Omega_liq_rms 1.031) while uy_liq_max is 4.3x and uy_liq_rms 1.27x.
Large vertical velocity with essentially unchanged vorticity is either
(a) a larger, near-irrotational standing wave, or (b) spurious currents in
interfacial cells. vol_frac_interf.dat separates them: a real wave moves
posY_max/posY_min and the interfacial area; parasitic currents do not.

It also prints the late-window trend, because the tau_100 ratio was still
climbing at cycle 116 (1.25 at cycle 36 -> 5.0 at cycle 116), which would
mean "the two branches' converged values" is the wrong frame -- one of them
may not be converged.

Usage:
    uv run python scripts/analyze_branch_interface.py
"""
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
REF, ESC = "kicktest_L7_th7", "c2001da9"


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(math.radians(theta)))
    return T_per / (L_bio / (V_bio / (H_bio * 0.5) / T_per))


T_ND = t_per_nd(32.5, 7.0)


def cycle_stat(path, col, fn=np.max):
    d = np.loadtxt(path, skiprows=1)
    t, y = d[:, 1], d[:, col]
    cyc = (t - t[0]) / T_ND
    return np.array([fn(y[(cyc >= c) & (cyc < c + 1)]) for c in range(int(cyc[-1]))])


print("vol_frac_interf.dat, mean of last 20 cycles, escaped / reference")
for col, name, fn in [(2, "f_liq_sum (mass)", np.mean), (3, "f_liq_interf (area)", np.max),
                      (4, "posY_max", np.max), (5, "posY_min", np.min)]:
    a = cycle_stat(RUNS / REF / "vol_frac_interf.dat", col, fn)
    b = cycle_stat(RUNS / ESC / "vol_frac_interf.dat", col, fn)
    print(f"  {name:22s} ref={np.mean(a[-20:]):11.5g}  esc={np.mean(b[-20:]):11.5g}"
          f"  ratio={np.mean(b[-20:]) / np.mean(a[-20:]):7.4f}")

# Wave height is the physically meaningful combination.
for tag, run in [("ref", REF), ("esc", ESC)]:
    hi = cycle_stat(RUNS / run / "vol_frac_interf.dat", 4, np.max)
    lo = cycle_stat(RUNS / run / "vol_frac_interf.dat", 5, np.min)
    print(f"  wave height ({tag}): {np.mean((hi - lo)[-20:]):.5g}")

print("\nis the escaped branch converged? per-20-cycle-block means, escaped run")
for col, name in [(4, "tau_100"), (5, "tau_mean")]:
    p = cycle_stat(RUNS / ESC / "shear_stress.dat", col)
    blocks = [f"{np.mean(p[i:i + 20]):.4g}" for i in range(0, len(p) - 19, 20)]
    print(f"  {name:9s} " + "  ".join(blocks))
for col, name in [(4, "tau_100"), (5, "tau_mean")]:
    p = cycle_stat(RUNS / REF / "shear_stress.dat", col)
    blocks = [f"{np.mean(p[i:i + 20]):.4g}" for i in range(0, len(p) - 19, 20)]
    print(f"  {name:9s} " + "  ".join(blocks) + "   (reference)")
