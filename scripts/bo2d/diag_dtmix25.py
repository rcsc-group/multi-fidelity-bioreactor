"""Fill the dtmix cache for every database / campaign / BO run and check the recomputed dtmix_0.50 against the
database values (kla_mix_dataset.csv, results.json). Usage: uv run python scripts/bo2d/diag_dtmix25.py"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import loop as L  # noqa: E402

C = L.C
d25 = C.load_db("dtmix_0.25")
d50 = C.load_db("dtmix_0.50")
bo = L.rows("dtmix_0.25", ["bo1_", "bo4_", "bo_l7_"], {7, 8, 9})
print("db rows with dtmix_0.25:", len(d25), "by level", d25.level.value_counts().to_dict())
print("bo rows with dtmix_0.25:", len(bo), "by level", bo.level.value_counts().to_dict())
m = d25.merge(d50, on="run", suffixes=("25", "50db"))
recomp = np.array([C.mix_value(C.ROOT / "runs" / r, "dtmix_0.50") or np.nan for r in m.run])
ok = np.isfinite(recomp)
rel = recomp[ok] / m.y50db.to_numpy(float)[ok] - 1
print(f"dtmix_0.50 recomputed vs database: n={ok.sum()} median rel diff {np.median(rel):+.3f}, "
      f"|diff|>10%: {(abs(rel) > 0.10).sum()}")
for r, a, b in zip(m.run[ok], recomp[ok], m.y50db[ok]):
    if abs(a / b - 1) > 0.10:
        print(f"   {r:32s} recomputed {a:7.1f}  db {b:7.1f}")
print("ratio dtmix_0.50/dtmix_0.25 median:", round(float(np.nanmedian(m.y50db / m.y25)), 2))
