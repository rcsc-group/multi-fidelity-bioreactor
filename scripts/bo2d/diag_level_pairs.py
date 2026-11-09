"""Do cheap levels carry information about L9 in the REAL data? Pair every L9 kLa_1T_10 / tau_95 row with the
L8 and L7 rows at the same (rpm, theta) (|d rpm| < 0.05, |d theta| < 0.02) and report Spearman, Pearson (ln),
and the sd of ln(L9/Lk). Usage: uv run python scripts/bo2d/diag_level_pairs.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import loop as L  # noqa: E402

C = L.C


def table(q):
    db = C.load_db(q)[["level", "rpm", "theta", "run", "y"]]
    extra = L.rows(q, ["bo1_", "bo4_", "bo_l7_"], {7, 8, 9})[["level", "rpm", "theta", "run", "y"]]
    d = pd.concat([db, extra]).drop_duplicates("run")
    return d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]


QS = sys.argv[1:] or ["kLa_1T_10", "tau_95"]
for q in QS:
    d = table(q)
    d9 = d[d.level == 9]
    for lk in (8, 7):
        dk = d[d.level == lk]
        pr = []
        for _, r in d9.iterrows():
            m = dk[(abs(dk.rpm - r.rpm) < 0.05) & (abs(dk.theta - r.theta) < 0.02)]
            if len(m):
                pr.append((r.rpm, r.theta, r.y, float(m.y.mean())))
        if len(pr) < 3:
            print(q, f"L{lk}: {len(pr)} pairs"); continue
        p = pd.DataFrame(pr, columns=["rpm", "theta", "y9", "yk"])
        a, b = np.log(p.y9), np.log(p.yk)
        print(f"{q} L{lk}-L9: n={len(p)} Spearman {p.y9.corr(p.yk, method='spearman'):+.2f} "
              f"Pearson(ln) {a.corr(b):+.2f} sd ln(L9/L{lk}) {np.std(a - b):.3f} sd ln L9 {np.std(a):.3f}")
        if q == "kLa_1T_10":
            print(p.sort_values(["theta", "rpm"]).round(3).to_string(index=False))
