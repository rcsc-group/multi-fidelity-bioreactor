"""Does L8 data help predict L9? The predictive test (replaces the Spearman test, which only asked whether the two
levels rank points in the same order; the MF model HF = rho*LF + beta (+ theta) + GP discrepancy does not need that).
Leave one L9 point out; predict it from the other L9 points with
  L9-only : constant mean + GP (the locked L9-only model)
  MF      : L8 GP (all L8 rows) as the basis [1, f8] + GP discrepancy   (locked chain L9 layer without theta)
  MF+th   : basis [1, f8, theta] + GP                                   (locked chain L9 layer, CHAIN_THETA)
Reported per QoI: LOO RMSE (ln), mean log predictive density, and the LOO predictions for the figure (csv).
Usage: uv run python scripts/bo2d/diag_mf_loo.py [qoi ...]"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special, stats

sys.path.insert(0, str(Path(__file__).parent))
import loop as L  # noqa: E402
from scripts.bo2d import bo  # noqa: E402

C = L.C
OUT = Path(__file__).parent / "state" / "diag_mf_loo.csv"


def data(q):
    db = C.load_db(q)[["level", "rpm", "theta", "run", "y"]]
    extra = L.rows(q, ["bo1_", "bo4_"], {8, 9})[["level", "rpm", "theta", "run", "y"]]
    d = pd.concat([db, extra]).drop_duplicates("run")
    return d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]


def lpd(m, x, y):
    mu, sc, nu, w = m.predict_dist(x[None])
    return float(special.logsumexp(np.log(w) + stats.t.logpdf(y, nu, mu[:, 0], sc[:, 0])))


def main():
    rows = []
    for q in sys.argv[1:] or ["kLa_1T_10", "tau_95", "edr_95", "dtmix_0.50"]:
        d = data(q); floor = C.noise_floor(q)
        d8, d9 = d[d.level == 8], d[d.level == 9]
        X8, y8 = d8[["rpm", "theta"]].to_numpy(float), np.log(d8.y.to_numpy(float))
        X9, y9 = d9[["rpm", "theta"]].to_numpy(float), np.log(d9.y.to_numpy(float))
        m8 = bo.FullBayesKRRLRGP(C.DOMAIN, basis="ordinary", grid=C.hf_grid(floor, float(np.std(y8)))).train(X8, y8)
        lf = bo._Mean(m8)
        res = {k: ([], []) for k in ("L9-only", "MF", "MF+th")}
        for i in range(len(y9)):
            k = np.arange(len(y9)) != i
            g = C.hf_grid(floor, float(np.std(y9[k])))
            ms = {"L9-only": bo.FullBayesKRRLRGP(C.DOMAIN, basis="ordinary", grid=g),
                  "MF": bo.FullBayesKRRLRGP(C.DOMAIN, lf_model=lf, grid=g),
                  "MF+th": bo._ThetaBasis(C.DOMAIN, lf_model=lf, grid=g)}
            for name, m in ms.items():
                m.train(X9[k], y9[k])
                pm = float(m.predict_mean(X9[i:i + 1])[0])
                res[name][0].append(pm - y9[i]); res[name][1].append(lpd(m, X9[i], y9[i]))
                rows.append(dict(qoi=q, model=name, rpm=X9[i, 0], theta=X9[i, 1], y9=y9[i], pred=pm,
                                 f8=float(lf.predict(X9[i:i + 1])[0])))
        print(f"{q}: n9={len(y9)} n8={len(y8)} sd(ln L9)={np.std(y9):.3f}")
        for name, (e, lp) in res.items():
            print(f"   {name:8s} RMSE {np.sqrt(np.mean(np.square(e))):.3f}  logpd {np.mean(lp):+.3f}")
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
