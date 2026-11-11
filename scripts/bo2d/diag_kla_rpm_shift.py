"""Hypothesis (2026-11-11 08:40): on theta = 7 the L8 -> L9 kLa mismatch is a regime change the coarse grid places at a
different rpm. Test: ln kLa9(r) = a + b ln kLa8(s * r), kLa8 interpolated (log-linear in rpm) from the dense theta = 7 L8
rows (truth snapshot, median per rpm). Scan s in [0.8, 1.25]; report LOO RMSE (ln) over the 10 fig9 L9 sites vs s = 1
(no shift) and vs the mean. Same for tau_95 as a control (expect best s ~ 1).
  uv run python scripts/bo2d/diag_kla_rpm_shift.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

S = Path(__file__).parent / "state"


def series(q, lev):
    d = pd.read_csv(S / f"truth_snapshot_{q}.csv")
    d = d[(d.level == lev) & np.isclose(d.theta, 7) & np.isfinite(d.y) & (d.y > 0)]
    g = d.groupby(d.rpm.round(3)).y.median()
    return g.index.to_numpy(float), np.log(g.to_numpy(float))


def loo(f, y):
    A = np.column_stack([np.ones_like(f), f]); e = []
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None); e.append(y[i] - A[i] @ c)
    return float(np.sqrt(np.mean(np.square(e))))


for q in ("kLa_1T_10", "tau_95"):
    r8, l8 = series(q, 8); r9, l9 = series(q, 9)
    res = []
    for s in np.round(np.arange(0.80, 1.251, 0.01), 2):
        x = s * r9
        ok = (x >= r8.min()) & (x <= r8.max())
        if ok.sum() < 7:
            continue
        res.append((loo(np.interp(x[ok], r8, l8), l9[ok]), s, int(ok.sum()), float(np.std(l9[ok]))))
    best = min(res); one = [r for r in res if r[1] == 1.0][0]
    print(f"{q:10s} n9={len(r9)}  no shift: LOO {one[0]:.3f} (n {one[2]}, sd {one[3]:.3f})  best s={best[1]:.2f}: "
          f"LOO {best[0]:.3f} (n {best[2]}, sd {best[3]:.3f})")
    print("   s:LOO " + " ".join(f"{s:.2f}:{e:.2f}" for e, s, _, _ in res[::3]))
