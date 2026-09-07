"""Is the branch offset in the FLOW, or only in the averaged diagnostics?

analyze_branch_separation.py showed every nonzero-Delta_theta restart lands
on tau_mean/tau_ref = 1.253 +/- 0.003, identical across Delta_theta = 0.1 to
5deg and across a 3- vs 40-cycle ramp -- a switch, not a response. Two
families of explanation:

  (a) the flow really differs   -> every column (maxima AND means AND the
      independent interface/norm diagnostics) shifts together
  (b) a diagnostic masking or normalization bug on the restart path -> the
      averaged columns shift but the pure maxima (tau_100) do not, or the
      shift differs column to column

Usage:
    uv run python scripts/analyze_branch_all_columns.py
"""
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"

COLS = {2: "tau_95", 3: "tau_98", 4: "tau_100", 5: "tau_mean",
        6: "tau_100_strict", 7: "tau_mean_strict", 9: "ediss_mean"}


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(math.radians(theta)))
    U_bio = V_bio / (H_bio * 0.5) / T_per
    return T_per / (L_bio / U_bio)


T_ND = t_per_nd(32.5, 7.0)


def late_mean(path, col, n_last=20):
    d = np.loadtxt(path, skiprows=1)
    t, y = d[:, 1], d[:, col]
    t0 = t[0]
    n = int((t[-1] - t0) / T_ND)
    cyc = (t - t0) / T_ND
    peaks = [y[(cyc >= c) & (cyc < c + 1)].max() for c in range(n)
             if ((cyc >= c) & (cyc < c + 1)).sum()]
    return float(np.mean(peaks[-n_last:]))


CASES = [("dtheta=0", "8fd81f04"), ("dtheta=0.1", "c2001da9"),
         ("dtheta=3", "1ddbba85"), ("dtheta=5", "3324f34b")]

print("per-cycle-peak, mean of last 20 cycles, as a ratio to the cold-start reference\n")
hdr = f"{'column':18s}" + "".join(f"{lbl:>12s}" for lbl, _ in CASES)
print(hdr)
for col, name in COLS.items():
    ref = late_mean(RUNS / "kicktest_L7_th7" / "shear_stress.dat", col)
    row = f"{name:18s}"
    for _, run in CASES:
        row += f"{late_mean(RUNS / run / 'shear_stress.dat', col) / ref:12.4f}"
    print(row)

# Independent diagnostics: interfacial area and the velocity norm.
for fname, col, name in [("vol_frac_interf.dat", None, "interfacial area"),
                         ("normf.dat", None, "normf")]:
    p = RUNS / "kicktest_L7_th7" / fname
    if not p.exists():
        continue
    print(f"\n--- {fname} header: {p.read_text().splitlines()[0]}")
