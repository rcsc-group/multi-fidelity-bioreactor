"""Numerical uncertainty of L10 statistics by least-squares grid fits.

A simplified version of Eca & Hoekstra (2014, J. Comput. Phys. 262:104),
the standard procedure for scattered or oscillatory multi-grid data where
three-grid Richardson extrapolation is undefined (R < 0 or R > 1).
With levels N, h_i = 2^(10 - N), so h_i = 1 at L10.
Candidate error models, each fitted by least squares:
    free p:   q = q0 + a h^p      (kept only if 0.5 <= p <= 2)
    p = 1:    q = q0 + a h
    p = 2:    q = q0 + a h^2
    1 & 2:    q = q0 + a h + b h^2
The model with the smallest fit standard deviation U_s is selected.
eps = |q10 - q0| is the error estimate, and Dq = (max q - min q)/(n - 1) is the data range.
    p in [0.5, 2.1) and U_s < Dq:   U = 1.25 eps + U_s + |q10 - fit(1)|
    U_s < Dq, other p:              U = 3 eps + U_s + |q10 - fit(1)|
    U_s >= Dq (fit is poor):        U = 3 U_s/Dq (eps + U_s + |q10 - fit(1)|)
U is an approximately 95% bound on |q10 - q_exact|.
Two sets of levels are reported: {6, 8, 9, 10} and {8, 9, 10}. If U
changes a lot when L6 is dropped, the estimate depends on a pre-asymptotic grid.
"""
from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.diag_richardson_means import LEVELS, RPMS, cycle_means  # noqa: E402

warnings.filterwarnings("ignore")
COLUMNS = [("tau_mean", "tau"), ("tau_98", "tau"), ("ediss_mean", "ediss")]


def fits(h, q):
    out = []
    if len(h) >= 4:
        try:
            (q0, a, p), _ = curve_fit(lambda h, q0, a, p: q0 + a * h**p, h, q,
                                      p0=(q[-1], q[0] - q[-1], 1.0), maxfev=20000)
            if 0.5 <= p <= 2:
                out.append(("p=%.2f" % p, p, lambda x, q0=q0, a=a, p=p: q0 + a * x**p, 3))
        except RuntimeError:
            pass
    for name, cols, p in (("p=1", [h], 1.0), ("p=2", [h**2], 2.0),
                          ("p=1&2", [h, h**2], 1.5)):
        A = np.column_stack([np.ones_like(h), *cols])
        if A.shape[1] >= len(h):
            continue
        c = np.linalg.lstsq(A, q, rcond=None)[0]
        k = len(cols)
        out.append((name, p, lambda x, c=c, k=k: c[0] + sum(c[j + 1] * x**(j + 1 if k == 2 else p)
                                                          for j in range(k)), A.shape[1]))
    best = None
    for name, p, f, npar in out:
        us = math.sqrt(np.sum((q - f(h))**2) / max(len(h) - npar, 1))
        if best is None or us < best[3]:
            best = (name, p, f, us)
    return best


def uncertainty(levels, q):
    h = np.array([2.0 ** (10 - N) for N in levels])
    q = np.array(q)
    name, p, f, us = fits(h, q)
    q0 = f(0.0)
    eps, dq, res = abs(q[-1] - q0), (q.max() - q.min()) / (len(q) - 1), abs(q[-1] - f(1.0))
    if us < dq:
        U = (1.25 if 0.5 <= p < 2.1 else 3.0) * eps + us + res
    else:
        U = 3 * us / dq * (eps + us + res)
    return name, q0, U


for col, kind in COLUMNS:
    print(f"\n== {col}")
    print(f"{'rpm':>5} {'L10':>9} {'SE10':>8} | {'fit':>7} {'q0':>9} {'U/L10':>6} "
          f"| {'fit':>7} {'q0':>9} {'U/L10':>6}   (levels 6,8,9,10 | 8,9,10)")
    for rpm in RPMS:
        try:
            m = {N: cycle_means(t.format(rpm), col, kind, s) for N, (t, s) in LEVELS.items()}
        except (FileNotFoundError, KeyError, ValueError):
            continue
        q = {N: v.mean() for N, v in m.items()}
        se = m[10].std(ddof=1) / math.sqrt(len(m[10]))
        row = f"{rpm:5g} {q[10]:9.4g} {se:8.2g}"
        for lv in ((6, 8, 9, 10), (8, 9, 10)):
            name, q0, U = uncertainty(lv, [q[N] for N in lv])
            row += f" | {name:>7} {q0:9.4g} {U / q[10]:6.1%}"
        print(row)
