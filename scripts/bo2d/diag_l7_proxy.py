"""Does L7 kLa_1T_25 (defined at all 36 L7 seeds) track L8/L9 kLa_1T_10 as well as L7 kLa_1T_10 does (2026-11-08)?
L7 kLa_1T_10 is NaN at 9/36 seeds (all rpm <= 21.4): the 10 % crossing comes within half a period of release."""
import sys
from pathlib import Path
import numpy as np
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import loop as L  # noqa: E402
C = L.C
db = C.load_db("kLa_1T_10")
key = lambda z: z.rpm.round(1).astype(str) + "_" + z.theta.round(1).astype(str)
for q7 in ("kLa_1T_10", "kLa_1T_25"):
    d7 = L.rows(q7, ["bo_l7_"], {7})
    for lev in (8, 9):
        d = db[db.level == lev]
        m = d7.assign(k=key).merge(d.assign(k=key), on="k", suffixes=("7", "x"))
        m = m[np.isfinite(np.log(m.y7)) & np.isfinite(np.log(m.yx))]; a, b = np.log(m.y7), np.log(m.yx)
        print(f"L7 {q7} vs L{lev} kLa_1T_10: n={len(m):2d} pearson {stats.pearsonr(a, b)[0]:.2f} "
              f"spearman {stats.spearmanr(a, b)[0]:.2f}")
