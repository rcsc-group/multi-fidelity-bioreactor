"""kLa L8->L9: are the transfer misfits concentrated at high rpm, where the 2026-11-06 finding says L9 kLa_1T_10 is
measured while the wave is still decaying (release at 12 warm cycles, settles after ~27)? Paired real sites, LOO
residual of y9 = a + b y8 + c theta, printed per site sorted by rpm. Usage: uv run python scripts/bo2d/diag_kla_residual_rpm.py"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import loop as LP  # noqa: E402
q = "kLa_1T_10"
d = pd.read_csv(Path("scripts/bo2d/state") / f"truth_snapshot_{q}.csv")
d = d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]
d["ly"] = np.log(d.y.astype(float)); d["key"] = d.rpm.round(2).astype(str) + "_" + d.theta.round(2).astype(str)
p = d[d.level == 8].groupby("key").agg(y8=("ly", "mean"), th=("theta", "first"), rpm=("rpm", "first")).join(
    d[d.level == 9].groupby("key").agg(y9=("ly", "mean"), run9=("run", "first")), how="inner").sort_values("rpm")
A = np.column_stack([np.ones(len(p)), p.y8, p.th]); y = p.y9.to_numpy(); e = []
for i in range(len(y)):
    m = np.arange(len(y)) != i; c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None); e.append(y[i] - A[i] @ c)
p["loo_res"] = np.round(e, 3)
print(p[["rpm", "th", "y8", "y9", "loo_res", "run9"]].round(3).to_string())
lo, hi = p.rpm < 25, p.rpm >= 25
print(f"rms LOO residual: rpm<25 {np.sqrt(np.mean(np.square(p.loo_res[lo]))):.3f} (n={lo.sum()}), "
      f"rpm>=25 {np.sqrt(np.mean(np.square(p.loo_res[hi]))):.3f} (n={hi.sum()})")
