"""How many runs per level does a good 2D surrogate need? (user, 2026-11-11 ~09:40: "200 for 2D, is that not overkill?
how many L9 are expected to be needed for a good surrogate? then backtrack the core-h budget".)

Fixed before any result:
A. L7 density. 200 scr_l7 runs; a fixed random 40 are the test set (seed 0). Train a GP (FullBayes, ordinary, ln y) on
   N7 in (10, 20, 40, 80, 160) random sites of the other 160 (3 draws each). Scores: (i) RMSE ln on the 40 L7 test
   sites (map error at L7 itself); (ii) "L7 map -> L9": LOO calibrated error of L9 ~ a + b * map over the 33 L9 sites
   (the number that matters for the chain), as a ratio to the mean-only error.
B. L9 count. 33 L9 sites (bo1 / bo4 / fig9; NOT space-filling: they sit near the tau cap and on theta = 7 - stated
   as a bias). L9-only GP on N9 in (6, 10, 15, 20, 25) random sites (5 draws), test = the rest; RMSE ln.
One QoI per process (memory). Writes state/learning_curve_<q>.json.
  uv run python scripts/bo2d/learning_curve.py tau_95
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.bo2d import loop as LP  # noqa: E402
from scripts.mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

C = LP.C
q = sys.argv[1]
# Same patches as diag_l7_screen_map (copied: importing that module RUNS its full grid map -> OOM, 2026-11-11 09:38)
C.warm_run_gate = lambda d: ("PASS", {})
_row = C.campaign_row


def _safe_row(*a, **k):
    try:
        return _row(*a, **k)
    except (ValueError, IndexError, KeyError, FileNotFoundError):
        return None


C.campaign_row = _safe_row


class _m:  # noqa: N801
    @staticmethod
    def clean(d):
        y = d.y.astype(float)
        return d[np.isfinite(y) & (y > 0)]


def fit(X, y):
    return FullBayesKRRLRGP(C.DOMAIN, basis="ordinary", grid=C.hf_grid(C.noise_floor(q), float(np.std(y)))).train(X, y)


def xy(d):
    return d[["rpm", "theta"]].to_numpy(float), np.log(d.y.to_numpy(float))


X7, y7 = xy(_m.clean(LP.as_objective(LP.rows(q, ["scr_l7_"], {7}), q)))
X9, y9 = xy(_m.clean(LP.as_objective(LP.rows(q, ["bo1_", "bo4_", "fig9_l9_rpm"], {9}), q)))
rng = np.random.default_rng(0)
perm = rng.permutation(len(y7)); te, pool = perm[:40], perm[40:]
out = {"q": q, "n7": len(y7), "n9": len(y9), "A": [], "B": []}


def calib(f, y):
    A = np.column_stack([np.ones_like(f), f]); e, e0 = [], []
    for i in range(len(y)):
        k = np.arange(len(y)) != i
        c, *_ = np.linalg.lstsq(A[k], y[k], rcond=None); e.append(y[i] - A[i] @ c); e0.append(y[i] - y[k].mean())
    return float(np.sqrt(np.mean(np.square(e))) / np.sqrt(np.mean(np.square(e0))))


for n in (10, 20, 40, 80, 160):
    for s in range(1 if n == 160 else 3):
        tr = np.random.default_rng(100 + s).choice(pool, min(n, len(pool)), replace=False)   # kLa: only 156 finite L7
        m = fit(X7[tr], y7[tr])
        e7 = float(np.sqrt(np.mean((np.ravel(m.predict_mean(X7[te])) - y7[te]) ** 2)))
        r9 = calib(np.ravel(m.predict_mean(X9)), y9)
        out["A"].append(dict(n=n, seed=s, rmse_l7=e7, sd_l7=float(np.std(y7[te])), l9_calib_ratio=r9))
        print(f"A {q} N7 {n:3d} draw {s}: L7 test RMSE ln {e7:.3f} (sd {np.std(y7[te]):.3f})  L7map->L9 calib/mean {r9:.2f}", flush=True)
for n in (6, 10, 15, 20, 25):
    for s in range(5):
        tr = np.random.default_rng(200 + s).choice(len(y9), n, replace=False); ho = np.setdiff1d(np.arange(len(y9)), tr)
        m = fit(X9[tr], y9[tr])
        e = float(np.sqrt(np.mean((np.ravel(m.predict_mean(X9[ho])) - y9[ho]) ** 2)))
        out["B"].append(dict(n=n, seed=s, rmse_l9=e, sd_l9=float(np.std(y9[ho]))))
        print(f"B {q} N9 {n:2d} draw {s}: L9 held-out RMSE ln {e:.3f} (sd {np.std(y9[ho]):.3f})", flush=True)
(ROOT / f"scripts/bo2d/state/learning_curve_{q}.json").write_text(json.dumps(out, indent=1))
