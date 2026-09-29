"""Do mean / percentile statistics converge across L6, L8, L9, L10 at 7 deg?

Three levels fix (q_inf, C, p) exactly: zero degrees of freedom, no test.
Four levels give two triplets with different ratios, and the test is whether
they agree:
    A = (L6, L8, L10), r = 4    B = (L8, L9, L10), r = 2
Each gives R = (q3 - q2)/(q2 - q1). The asymptotic model q_N = q_inf + C h^p,
h = 2^-N, predicts R_A = 4^-p and R_B = 2^-p, i.e. p_A = p_B, with both R in (0,1).
Falsifier: either R outside (0,1), or |p_A - p_B| large. Then the model does not hold
from L6 onwards. If B alone is in (0,1), the asymptotic range may start at L8,
but then only a triplet with zero degrees of freedom is left.

Statistic: time average over the settled window, before release. The
cycle-to-cycle SE is reported so level jumps can be compared with it.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RPMS = [17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
LEVELS = {6: ("fig13a_l6_mf_rpm{:g}", 0.6), 8: ("fig13a_rm_mf_rpm{:g}", 0.6),
          9: ("fig9_l9_rpm{:g}", 0.6), 10: ("l10_fig13a_xlevel_rpm{:g}", 0.0)}
COLUMNS = [("tau_mean", "tau"), ("tau_95", "tau"), ("tau_98", "tau"),
           ("ediss_mean", "ediss")]


def cycle_means(run: str, column: str, kind: str, settled_from: float) -> np.ndarray:
    d = ROOT / "runs" / run
    p = json.load(open(d / "params.json"))
    L, H = p["geometry"]["a"], 2 * p["geometry"]["b"]
    Tp = 2 * math.pi / p["omega_b"]
    U = L / 4 * (H + 0.5 * L * math.tan(math.radians(p["theta_max"][0]))) / (H * 0.5) / Tp
    T_nd = Tp / (L / U)
    df = pd.read_csv(d / "shear_stress.dat", sep=r"\s+")
    t = df["t"].to_numpy()
    t_rel = (p.get("t_checkpoint") or 0) + p["n_mix_cycles"] * T_nd
    keep = t < t_rel if t_rel < t[-1] else np.ones_like(t, bool)
    t = t[keep]
    y = df[column].to_numpy()[keep] * (1000 * U**2 if kind == "tau" else 1000 * U**3 / L)
    c = (t - t[0]) / T_nd
    k0, k1 = int(math.ceil(settled_from * c[-1])), int(math.floor(c[-1]))
    return np.array([y[(c >= k) & (c < k + 1)].mean() for k in range(k0, k1)])


def triplet(q1, q2, q3, r):
    R = (q3 - q2) / (q2 - q1)
    return R, (-math.log(R) / math.log(r) if 0 < R < 1 else float("nan"))


for col, kind in COLUMNS:
    print(f"\n== {col}")
    print(f"{'rpm':>5} {'L6':>9} {'L8':>9} {'L9':>9} {'L10':>9} {'SE10':>8} "
          f"{'R_A':>6} {'p_A':>5} {'R_B':>6} {'p_B':>5} {'q_inf_B':>9} {'GCI_B':>8}")
    for rpm in RPMS:
        try:
            m = {N: cycle_means(t.format(rpm), col, kind, s) for N, (t, s) in LEVELS.items()}
        except (FileNotFoundError, KeyError, ValueError):
            continue
        q = {N: v.mean() for N, v in m.items()}
        se10 = m[10].std(ddof=1) / math.sqrt(len(m[10]))
        RA, pA = triplet(q[6], q[8], q[10], 4)
        RB, pB = triplet(q[8], q[9], q[10], 2)
        if 0 < RB < 1:
            qi = q[10] + (q[10] - q[9]) / (2**pB - 1)
            gci = 1.25 * abs(q[10] - q[9]) / (2**pB - 1)
            tail = f"{qi:9.4g} {gci:8.2g}"
        else:
            tail = f"{'-':>9} {'-':>8}"
        print(f"{rpm:5g} " + " ".join(f"{q[N]:9.4g}" for N in (6, 8, 9, 10))
              + f" {se10:8.2g} {RA:6.2f} {pA:5.2f} {RB:6.2f} {pB:5.2f} " + tail)
sys.exit(0)
