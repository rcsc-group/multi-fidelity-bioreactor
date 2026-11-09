"""Leave-one-out check of L9 GP variants on the real L9 kLa_1T_10 data (ln units).
Question: is the flat L9-only GP the data or the model? Variants:
  V0 ordinary (constant mean), noise ratio >= floor/ystd (floor 0.095)  [locked L9-only model]
  V1 ordinary, floor 0.04 (L9 replicate sd, diary 2026-10-23: 28.7/26.5/27.9/28.9)
  V2 linear trend in (rpm, theta) + GP, floor 0.095
  V3 linear trend, floor 0.04
Data sets: arm = 4 seed + L9-only arm's own runs; nofig9 = all L9 except theta-7 fig9; all = every L9 row.
Score: LOO RMSE of ln kLa, mean log predictive density (mixture of Student t), mean predictive sd, and
the spread (sd over a 41x41 grid) of the posterior mean: 'flatness'.
Usage: uv run python scripts/bo2d/diag_l9_loo.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special, stats

sys.path.insert(0, str(Path(__file__).parent))
import loop as L  # noqa: E402
from scripts.bo2d import bo  # noqa: E402
import loop as L  # noqa: E402

C = L.C
FB = bo.FullBayesKRRLRGP


class Lin(FB):
    def _basis(self, X):
        Xn = self.normalize_input(np.atleast_2d(X))
        return np.column_stack([np.ones(len(Xn)), Xn])


def make(var, y):
    floor = 0.095 if var in ("V0", "V2") else 0.04
    cls = FB if var in ("V0", "V1") else Lin
    return cls(C.DOMAIN, basis="ordinary", grid=C.hf_grid(floor, float(np.std(y))))


def logpdf(m, x, y):
    mu, sc, nu, w = m.predict_dist(x[None])
    return float(special.logsumexp(np.log(w) + stats.t.logpdf(y, nu, mu[:, 0], sc[:, 0])))


def main():
    q = "kLa_1T_10"
    db = C.load_db(q); db = db[db.level == 9][["rpm", "theta", "run", "y"]]
    bo1 = L.rows(q, ["bo1_"], {9})[["rpm", "theta", "run", "y"]]
    allr = pd.concat([db, bo1]).drop_duplicates("run").reset_index(drop=True)
    seed = allr[allr.run.isin(L.SEED9)]
    sets = {"arm": pd.concat([seed, bo1]).drop_duplicates("run"),
            "nofig9": allr[~allr.run.str.startswith("fig9")],
            "all": allr}
    g = np.array([[r, t] for r in np.linspace(*C.DOMAIN[0], 41) for t in np.linspace(*C.DOMAIN[1], 41)])
    print(f"{'set':7s} {'var':3s} {'n':>3s} {'RMSE':>6s} {'logpd':>7s} {'sd_pred':>7s} {'flat':>6s} {'ystd':>6s}")
    for name, d in sets.items():
        X = d[["rpm", "theta"]].to_numpy(float); y = np.log(d.y.to_numpy(float))
        for var in ("V0", "V1", "V2", "V3"):
            err, lp, sds = [], [], []
            for i in range(len(y)):
                k = np.arange(len(y)) != i
                m = make(var, y[k]).train(X[k], y[k])
                mu, sc, nu, w = m.predict_dist(X[i:i + 1])
                mean = float(w @ mu[:, 0])
                var_t = float(w @ (sc[:, 0] ** 2 * nu / (nu - 2) + mu[:, 0] ** 2)) - mean ** 2
                err.append(mean - y[i]); lp.append(logpdf(m, X[i], y[i])); sds.append(np.sqrt(var_t))
            mf = make(var, y).train(X, y)
            flat = float(np.std(mf.predict_mean(g)))
            print(f"{name:7s} {var:3s} {len(y):3d} {np.sqrt(np.mean(np.square(err))):6.3f} {np.mean(lp):7.3f} "
                  f"{np.mean(sds):7.3f} {flat:6.3f} {np.std(y):6.3f}")


if __name__ == "__main__":
    main()
