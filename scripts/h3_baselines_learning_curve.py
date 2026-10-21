"""Baselines and learning curve for the hydrodynamic chain (pre-registered, diary 2026-10-21).

(1) 13-run leave-one-out and the theta-5.5 test set: chain (L9 at the point given), L9 x constant ratio, L9 x LOG-LINEAR ratio
    (ln(L10/L9) = a + b rpm + c theta by least squares), L10-only GP, L8 x ratio.
(2) Learning curve: random training subsets of the 13 L10 runs (n = 3, 4, 6, 8, 10, 12; 15 draws), scored on the held-out runs.
(3) Variant C (no L9 at the held-out point): chain vs 'L9 GP interpolated x log-linear ratio'.
RULE: the chain has a REAL ADVANTAGE if its RMSE is >= 20 % lower than the best baseline for >= 2 of 3 quantities at n <= 6 AND on the full LOO.
Usage: uv run python scripts/h3_baselines_learning_curve.py
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from scripts.h3_band_check import Offset  # noqa: E402
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from scripts.test_h1_mfbml import BOUNDS, LF_THETAS, RPMS, hf_points, hf_run, lf_run, load_points, train_lf2d  # noqa: E402
from scripts.test_h2_mfbml import l9_points, l9_run  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

QS = ("tau_95_t", "tau_mean_t", "edr_mean_t")
SIZES = (3, 4, 6, 8, 10, 12)
DRAWS = 15


def design(X):
    return np.column_stack([np.ones(len(X)), X[:, 0], X[:, 1]])


def loglin_ratio(Xtr, y10, y9, Xte, y9te):
    c, *_ = np.linalg.lstsq(design(Xtr), np.log(y10) - np.log(y9), rcond=None)
    return np.exp(np.log(y9te) + design(Xte) @ c)


def rel_rmse(pred, obs):
    return float(np.sqrt(np.mean((np.asarray(pred) / np.asarray(obs) - 1) ** 2)))


def main():
    warnings.filterwarnings("ignore")
    lf_pts, _ = load_points([(r, t) for r in RPMS for t in LF_THETAS], lf_run)
    l9_pts, _ = load_points(l9_points(), l9_run)
    hf_pts, _ = load_points(hf_points(), hf_run)
    hk = sorted(hf_pts, key=lambda k: (k[1], k[0]))
    X = np.array(hk, float)
    out = {}
    test_rpm = [r for r in (17.5, 20, 22.5, 27.5, 30, 32.5) if (ROOT / "runs" / f"x8_l10_rpm{r:g}_th5.5" / "results.json").exists()
               and (ROOT / "runs" / f"x8_l9_rpm{r:g}_th5.5" / "results.json").exists()]
    Xt = np.array([[r, 5.5] for r in test_rpm])
    rng = np.random.default_rng(20261021)
    for q in QS:
        k8 = sorted(lf_pts, key=lambda k: (k[1], k[0]))
        lf = train_lf2d(np.array(k8), np.log([lf_pts[k][q] for k in k8]))
        k9 = sorted(l9_pts, key=lambda k: (k[1], k[0]))
        m1 = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(np.array(k9), np.log([l9_pts[k][q] for k in k9]))
        ad = Offset(m1)
        y10 = np.array([hf_pts[k][q] for k in hk])
        y9 = np.array([l9_pts[k][q] for k in hk])
        y8 = np.array([lf_pts[k][q] for k in hk])
        res = {}

        def chain_pred(tr, te):
            m2 = FullBayesKRRLRGP(BOUNDS, lf_model=ad, basis="linear").train(X[tr], np.log(y10[tr]))
            return np.exp(np.ravel(m2.predict_mean(X[te])))

        def baselines(tr, te):
            r9 = np.exp(np.mean(np.log(y10[tr] / y9[tr])))
            r8 = np.exp(np.mean(np.log(y10[tr] / y8[tr])))
            h = FullBayesKRRLRGP(BOUNDS, lf_model=None, basis="ordinary").train(X[tr], np.log(y10[tr]))
            b = dict(l9_const=y9[te] * r9, l8_const=y8[te] * r8, l10_only=np.exp(np.ravel(h.predict_mean(X[te]))))
            if tr.sum() >= 4:
                b["l9_loglin"] = loglin_ratio(X[tr], y10[tr], y9[tr], X[te], y9[te])
            return b

        # (1) leave-one-out
        acc = {k: [] for k in ("chain", "l9_const", "l9_loglin", "l10_only", "l8_const")}
        for i in range(len(hk)):
            tr = np.arange(len(hk)) != i
            te = ~tr
            acc["chain"].append(chain_pred(tr, te)[0] / y10[i] - 1)
            for k, v in baselines(tr, te).items():
                acc[k].append(v[0] / y10[i] - 1)
        res["loo"] = {k: float(np.sqrt(np.mean(np.square(v)))) for k, v in acc.items()}
        # test set
        if len(test_rpm):
            obs = np.array([qois(f"x8_l10_rpm{r:g}_th5.5")[q] for r in test_rpm])
            p9 = np.array([qois(f"x8_l9_rpm{r:g}_th5.5")[q] for r in test_rpm])
            tr = np.ones(len(hk), bool)
            m2 = FullBayesKRRLRGP(BOUNDS, lf_model=ad, basis="linear").train(X, np.log(y10))
            d9 = {**l9_pts, **{(r, 5.5): qois(f"x8_l9_rpm{r:g}_th5.5")[q] for r in test_rpm}}
            kk = sorted(d9, key=lambda k: (k[1], k[0]))
            m1t = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(np.array(kk), np.log([d9[k] for k in kk]))
            m2t = FullBayesKRRLRGP(BOUNDS, lf_model=Offset(m1t), basis="linear").train(X, np.log(y10))
            res["test"] = dict(n=len(test_rpm), chain=rel_rmse(np.exp(np.ravel(m2t.predict_mean(Xt))), obs),
                               l9_const=rel_rmse(p9 * np.exp(np.mean(np.log(y10 / y9))), obs),
                               l9_loglin=rel_rmse(loglin_ratio(X, y10, y9, Xt, p9), obs))
        # (2) learning curve
        lc = {}
        for n in SIZES:
            a = {k: [] for k in ("chain", "l9_const", "l9_loglin", "l10_only", "l8_const")}
            for _ in range(DRAWS):
                idx = rng.choice(len(hk), n, replace=False)
                tr = np.zeros(len(hk), bool)
                tr[idx] = True
                te = ~tr
                a["chain"].append(rel_rmse(chain_pred(tr, te), y10[te]))
                for k, v in baselines(tr, te).items():
                    a[k].append(rel_rmse(v, y10[te]))
            lc[n] = {k: float(np.mean(v)) for k, v in a.items() if len(v)}
        res["learning_curve"] = lc
        # (3) variant C: no L9 at the held-out point
        c_chain, c_base = [], []
        for i, k in enumerate(hk):
            tr = np.arange(len(hk)) != i
            k9c = [p for p in k9 if p != k]
            m1c = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(np.array(k9c), np.log([l9_pts[p][q] for p in k9c]))
            m2c = FullBayesKRRLRGP(BOUNDS, lf_model=Offset(m1c), basis="linear").train(X[tr], np.log(y10[tr]))
            c_chain.append(float(np.exp(m2c.predict_mean(X[i:i + 1])[0])) / y10[i] - 1)
            g9 = FullBayesKRRLRGP(BOUNDS, lf_model=None, basis="ordinary").train(np.array(k9c), np.log([l9_pts[p][q] for p in k9c]))
            p9i = float(np.exp(g9.predict_mean(X[i:i + 1])[0]))
            c_base.append(float(loglin_ratio(X[tr], y10[tr], y9[tr], X[i:i + 1], np.array([p9i]))[0]) / y10[i] - 1)
        res["variant_C"] = dict(chain=float(np.sqrt(np.mean(np.square(c_chain)))), l9gp_loglin=float(np.sqrt(np.mean(np.square(c_base)))))
        out[q] = res
        print(f"\n{q}\n  LOO RMSE: " + ", ".join(f"{k} {v:.1%}" for k, v in res["loo"].items()))
        if "test" in res:
            print(f"  test set (n={res['test']['n']}): chain {res['test']['chain']:.1%}  L9 const {res['test']['l9_const']:.1%}  L9 loglin {res['test']['l9_loglin']:.1%}")
        print("  variant C (no L9 at the point): chain %.1f%%  L9-GP x loglin %.1f%%" % (100 * res["variant_C"]["chain"], 100 * res["variant_C"]["l9gp_loglin"]))
        for n, d in lc.items():
            print(f"  n={n:2d}: " + "  ".join(f"{k} {v:.1%}" for k, v in d.items()))
    # rule
    adv = 0
    for q in QS:
        r = out[q]
        base_loo = min(r["loo"][k] for k in ("l9_const", "l9_loglin", "l10_only", "l8_const"))
        small = [r["learning_curve"][n] for n in (3, 4, 6)]
        ok_small = all(d["chain"] <= 0.8 * min(v for k, v in d.items() if k != "chain") for d in small)
        adv += (r["loo"]["chain"] <= 0.8 * base_loo) and ok_small
    out["verdict"] = "REAL ADVANTAGE" if adv >= 2 else "NO ACCURACY ADVANTAGE (intervals only)"
    print("\nverdict:", out["verdict"], f"({adv} of 3 quantities meet the 20 % rule)")
    (ROOT / "experiments/multifidelity/h3_baselines_learning_curve.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
