"""Falsification (2026-11-04): is the dedup LF optimum (lengthscales 0.02-0.05, noise 2e-4) stable across the
LF seed? Refits the LF only on the dedup kLa data for seeds 1..5 and prints lengthscales, noise and crit at (15, 2).
Compute node only. Usage: uv run python scripts/campaign2d/diag_kla_seed.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import campaign as c  # noqa: E402
from scripts.test_t3_learning_curve import CachedLooKRR  # noqa: E402

kla = pd.read_csv(c.DB_KLA)
d = kla[["level", "rpm", "theta", "run", "kLa_1T_10"]].rename(columns={"kLa_1T_10": "y"}).dropna()
d = d[d.level.isin([8, 9, 10]) & (d.rpm >= 15) & (d.rpm <= 37.5) & (d.theta >= 2) & (d.theta <= 7) & (d.y > 0)]
for label, dd in (("dup", d), ("dedup", c.dedup_runs(d))):
  lo = dd[dd.level == 8]
  X = lo[["rpm", "theta"]].to_numpy(float)
  Y = np.log(lo.y.to_numpy(float)).reshape(-1, 1)
  for seed in (42, 1, 2, 3):
    print(label, "rows", len(lo), end=" | ")
    lf = CachedLooKRR(design_space=c.DOMAIN, params_optimize=True, noise_data=True, optimizer_restart=10, seed=seed)
    lf.train(X, Y, portion_test=0.2)
    ls = np.round(np.asarray(lf.kernel.param, float), 4).tolist()
    nz = round(float(np.asarray(lf.noise_std).ravel()[0]), 5)
    lo_pred = lf.predict(np.array([[15.0, 2.0]])) if hasattr(lf, "predict") else None
    print(f"seed {seed}: lengthscale params {ls}, noise_std {nz}, prediction at (15,2) {np.round(lo_pred, 3).tolist() if lo_pred is not None else 'n/a'}")
