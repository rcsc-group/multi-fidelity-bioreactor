"""Two checks (2026-11-10, user questions):
(a) How long does ONE mfbml fit take? Time fit_single (L9, n=8) and fit_chain3 (L7/L8 at 44 sites, L9 at 8) and one
    prediction on the 101 x 101 grid.
(b) Is the truth fair to the chain, judged by the TRANSFER FUNCTION (what KRR-LR-GPR uses), not by correlation?
    mfbml's L9 layer = basis [1, f8, theta] times coefficients + GP discrepancy. Fit y9 = a + b y8 + c theta by least
    squares: (i) at paired real sites (observed), (ii) on the truth grid. Report b, c, and the discrepancy fraction
    R = var(residual) / var(y9). If R_truth >> R_obs the truth breaks the transfer function that the data have.
Usage: uv run python scripts/bo2d/diag_fit_time_transfer.py"""
import os
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import synth_bench as S  # noqa: E402
from scripts.bo2d import bo  # noqa: E402
from scripts.bo2d import loop as LP  # noqa: E402


def lsq(y9, y8, th):
    A = np.column_stack([np.ones_like(y8), y8, th]); c, *_ = np.linalg.lstsq(A, y9, rcond=None)
    r = y9 - A @ c
    return c, r.var() / y9.var()


G = bo.grid()
for obj, q in (("kla", "kLa_1T_10"), ("mix", "dtmix_0.25")):
    os.environ["BO_OBJ"] = obj; S.setup("real")
    for name, i in ((q, 0), ("tau_95", 1)):
        d = pd.read_csv(Path(__file__).parent / "state" / f"truth_snapshot_{name}.csv")
        d = d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]; d = LP.as_objective(d, name)
        d["ly"] = np.log(d.y.astype(float)); d["key"] = d.rpm.round(2).astype(str) + "_" + d.theta.round(2).astype(str)
        p = d[d.level == 8].groupby("key").agg(y8=("ly", "mean"), th=("theta", "first")).join(
            d[d.level == 9].groupby("key").ly.mean().rename("y9"), how="inner")
        co, Ro = lsq(p.y9.to_numpy(), p.y8.to_numpy(), p.th.to_numpy())
        ct, Rt = lsq(S.F[9][i](G), S.F[8][i](G), G[:, 1])
        print(f"{obj} {name}: observed n={len(p)} b={co[1]:.2f} c={co[2]:+.3f} R={Ro:.2f} | truth grid b={ct[1]:.2f} "
              f"c={ct[2]:+.3f} R={Rt:.2f}", flush=True)
# (a) timing
os.environ["BO_OBJ"] = "mix"; S.setup("real"); C = S.C
rng = np.random.default_rng(0); lo, hi = C.DOMAIN[:, 0], C.DOMAIN[:, 1]
X = lo + (hi - lo) * rng.random((44, 2)); K, T = {}, {}
S.add(K, T, 9, X[:8], rng); S.add(K, T, 8, X, rng); S.add(K, T, 7, X, rng)
t0 = time.time(); m = bo.fit_single(T[9], "tau_95"); t1 = time.time(); bo.single_moments(m, G); t2 = time.time()
ms = bo.fit_chain3(T[7], T[8], T[9], "tau_95"); t3 = time.time(); bo.chain_moments(ms, G); t4 = time.time()
print(f"timing: fit_single {t1-t0:.1f} s, predict grid {t2-t1:.1f} s | fit_chain3 {t3-t2:.1f} s, predict grid {t4-t3:.1f} s")
print(f"hyper-grid nodes per model: {len(m.w)}")
