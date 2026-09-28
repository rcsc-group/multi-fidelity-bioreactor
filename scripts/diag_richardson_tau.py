"""Is Fig 13(a) tau_max in the asymptotic range across L6 -> L8 -> L10?

Level N has h = L / 2^N, so L6, L8, L10 are a constant-ratio triplet (r = 4).
If q_N = q_inf + C h_N^p, the convergence ratio
    R = (q10 - q8) / (q8 - q6) = r^-p
must lie in (0, 1); then p = -ln R / ln r and
    q_inf = q10 + (q10 - q8) / (r^p - 1),   GCI_10 = Fs |q10 - q8| / (r^p - 1).
Falsifier: R <= 0 (oscillatory) or R >= 1 (diverging) at a point means no
asymptotic range there, and the discretisation error of L10 cannot be bounded
from these three levels.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.plot_mf_fig13 import per_cycle_max  # noqa: E402

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
LEVELS = {6: ("fig13a_l6_mf_rpm{:g}", 0.6), 8: ("fig13a_rm_mf_rpm{:g}", 0.6),
          9: ("fig9_l9_rpm{:g}", 0.6), 10: ("l10_fig13a_xlevel_rpm{:g}", 0.0)}
r, FS = 4.0, 1.25

print(f"{'rpm':>5} {'L6':>13} {'L8':>13} {'L10':>13} {'R':>7} {'p':>5} "
      f"{'q_inf':>7} {'GCI10':>7} {'SE10':>7}")
for rpm in RPMS:
    q, se = {}, {}
    try:
        for N, (tmpl, s) in LEVELS.items():
            m = per_cycle_max(tmpl.format(rpm), "tau_100_signed", "tau", s)
            q[N], se[N] = m.mean(), m.std(ddof=1) / math.sqrt(len(m))
    except (FileNotFoundError, ValueError, KeyError):
        continue
    R = (q[10] - q[8]) / (q[8] - q[6])
    cell = " ".join(f"{q[N]:.3f}±{se[N]:.3f}" for N in (6, 8, 10))
    if 0 < R < 1:
        p = -math.log(R) / math.log(r)
        qi = q[10] + (q[10] - q[8]) / (r ** p - 1)
        gci = FS * abs(q[10] - q[8]) / (r ** p - 1)
        print(f"{rpm:5g} {cell} {R:7.3f} {p:5.2f} {qi:7.3f} {gci:7.3f} {se[10]:7.3f}")
    else:
        print(f"{rpm:5g} {cell} {R:7.3f}   -- no asymptotic range --       {se[10]:7.3f}")

# Second, finer triplet at constant ratio r = 2: L8 -> L9 -> L10.
print(f"\n{'rpm':>5} {'L9':>13} {'R_8,9,10':>9} {'p':>5}")
for rpm in RPMS:
    try:
        q = {N: per_cycle_max(LEVELS[N][0].format(rpm), "tau_100_signed", "tau",
                              LEVELS[N][1]) for N in (8, 9, 10)}
    except (FileNotFoundError, ValueError, KeyError):
        continue
    m = {N: v.mean() for N, v in q.items()}
    R2 = (m[10] - m[9]) / (m[9] - m[8])
    p2 = -math.log(R2) / math.log(2) if 0 < R2 < 1 else float("nan")
    print(f"{rpm:5g} {m[9]:.3f}±{q[9].std(ddof=1)/math.sqrt(len(q[9])):.3f} {R2:9.3f} {p2:5.2f}")
