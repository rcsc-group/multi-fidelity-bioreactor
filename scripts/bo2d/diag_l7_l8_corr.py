"""How well does L7 track L8 (and L9) at the same (rpm, theta)? ln-scale Pearson and Spearman correlation and the
mean ln ratio, for kLa_1T_10 and tau_95 (2026-11-07). The chain only helps if L7 ranks points like L8/L9."""
import sys
from pathlib import Path
import numpy as np
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import loop as L  # noqa: E402
C = L.C
for q in (L.QK, L.QT):
    d7 = L.rows(q, ["bo_l7_"], {7}); db = C.load_db(q)
    for lev in (8, 9):
        d = db[db.level == lev]
        m = d7.assign(k=lambda z: z.rpm.round(1).astype(str) + "_" + z.theta.round(1).astype(str)).merge(
            d.assign(k=lambda z: z.rpm.round(1).astype(str) + "_" + z.theta.round(1).astype(str)), on="k", suffixes=("7", "x"))
        if len(m) < 3:
            print(q, f"L7 vs L{lev}: n={len(m)}"); continue
        m = m[np.isfinite(np.log(m.y7)) & np.isfinite(np.log(m.yx))]; a, b = np.log(m.y7), np.log(m.yx)
        print(f"{q:10s} L7 vs L{lev}: n={len(m):2d} pearson {stats.pearsonr(a, b)[0]:.2f} spearman {stats.spearmanr(a, b)[0]:.2f} "
              f"mean ratio L7/L{lev} {np.exp((a - b).mean()):.2f}")
