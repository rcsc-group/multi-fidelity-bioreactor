"""H2 evaluation (pre-registered, diary.md): recursive mfbml, L8 -> L9 -> L10, on (rpm, theta).

Rule (quoted): "recursive KRR-LR-GPR. Model 1 = FullBayesKRRLRGP L8 -> L9 (LF KRR on the 40 L8 points;
trained on all available L9 points: theta 7 fig9_l9 at 10 rpm + 6 new x8_l9 warm starts at the H1
(rpm, theta)). Model 2 = FullBayesKRRLRGP L9 -> L10 with lf_model = model 1 predictive mean, trained and
scored by the SAME LOO over the 13 HF points and the SAME PASS rule as H1 (tau_95 AND EDR: rel RMSE <=
10%, below HF-only GP and below L8 x ratio, coverage >= 11/13). Model 1 keeps the L9 point at the
held-out x (deployment assumption)."
PENDING while fewer than 13 HF runs have runs/<id>/results.json (same logic as H1).
Decisive baselines: HF-only GP, L8 x global ratio. Non-decisive (reported only): L9 x global ratio and
the plain L9 value as predictor.

Usage: uv run python scripts/test_h2_mfbml.py
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from scripts.test_h1_mfbml import (  # noqa: E402,F401
    BOUNDS, CRITERIA_QOIS, LF_THETAS, N_REQUIRED, QOIS, RPMS, h1_verdict, hf_points, hf_run, lf_run,
    load_points, ratio_baseline, rel_rmse_ratio, train_lf2d,
)

OUT = ROOT / "experiments" / "multifidelity" / "test_h2_mfbml.json"
X8_RPMS = (17.5, 25, 32.5)
X8_THETAS = (4, 9)


def l9_run(rpm, theta):
    if theta == 7:
        return f"fig9_l9_rpm{rpm:g}"
    return f"x8_l9_rpm{rpm:g}_th{theta:g}"


def l9_points():
    pts = [(float(r), 7.0) for r in RPMS]
    pts += [(float(r), float(t)) for r in X8_RPMS for t in X8_THETAS]
    return sorted(pts, key=lambda p: (p[1], p[0]))


h2_verdict = h1_verdict   # same PASS rule as H1 (pre-registered)


class MeanAdapter:
    """Exposes a trained FullBayesKRRLRGP as an LF model: .predict(X) -> predictive mean, shape (n,)."""

    def __init__(self, model):
        self.model = model

    def predict(self, X):
        return np.asarray(self.model.predict_mean(np.atleast_2d(np.asarray(X, float))), float).ravel()


def _fit_predict(X, y, xt, lf_model, basis, log=False):
    """log=True: train on log(y); return exp(posterior mean), exp(2.5/97.5% quantiles) (H3)."""
    from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP
    if log:
        y = np.log(y)
    m = FullBayesKRRLRGP(BOUNDS, lf_model=lf_model, basis=basis).train(X, y)
    qq = m.quantile(xt, [0.025, 0.975])[0]
    mean = float(m.predict_mean(xt)[0])
    if log:
        return float(np.exp(mean)), float(np.exp(qq[0])), float(np.exp(qq[1]))
    return mean, float(qq[0]), float(qq[1])


def loo_h2(q, lf_model, l9_pts, lf_pts, hf_pts, log=False):
    """LOO over hf_pts (each must have an L8 and L9 counterpart). lf_model = trained L8 KRR.
    Model 1 (L8 -> L9) is trained once on ALL l9_pts (it never sees L10).
    log=True (H3): model 1, model 2 and HF-only work on log(QoI); the caller must pass an lf_model
    trained on log L8. Predictions/intervals are returned on the raw scale; ratio and L9 columns stay raw."""
    from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP
    k9 = sorted(l9_pts, key=lambda k: (k[1], k[0]))
    model1 = FullBayesKRRLRGP(BOUNDS, lf_model=lf_model, basis="linear").train(
        np.array(k9), (np.log if log else np.asarray)([l9_pts[k][q] for k in k9]))
    adapter = MeanAdapter(model1)
    hk = sorted(hf_pts, key=lambda k: (k[1], k[0]))
    X = np.array(hk)
    yh = np.array([hf_pts[k][q] for k in hk])
    y8 = np.array([lf_pts[k][q] for k in hk])
    y9 = np.array([l9_pts[k][q] for k in hk])
    rows = []
    for i, k in enumerate(hk):
        tr = np.arange(len(hk)) != i
        xt = X[i:i + 1]
        r = dict(rpm=k[0], theta=k[1], obs=float(yh[i]), l8=float(y8[i]), l9=float(y9[i]))
        for name, lfm, basis in (("mf", adapter, "linear"), ("hf", None, "ordinary")):
            r[f"{name}_pred"], r[f"{name}_lo"], r[f"{name}_hi"] = _fit_predict(
                X[tr], yh[tr], xt, lfm, basis, log)
        r["ratio_pred"] = ratio_baseline(yh[tr], y8[tr], y8[i])        # decisive baseline (as H1)
        r["l9ratio_pred"] = ratio_baseline(yh[tr], y9[tr], y9[i])      # non-decisive
        r["l9_pred"] = float(y9[i])                                    # non-decisive
        rows.append(r)
    obs = np.array([r["obs"] for r in rows])
    hit = lambda n: int(sum(r[f"{n}_lo"] <= r["obs"] <= r[f"{n}_hi"] for r in rows))  # noqa: E731
    rm = lambda n: rel_rmse_ratio([r[f"{n}_pred"] for r in rows], obs)  # noqa: E731
    s = dict(n=len(rows), mf_rmse=rm("mf"), hf_rmse=rm("hf"), ratio_rmse=rm("ratio"),
             mf_cover=hit("mf"), hf_cover=hit("hf"),
             l9_ratio_rmse=rm("l9ratio"), l9_rmse=rm("l9"), n_l9_model1=len(k9))
    return s, rows


def main():
    warnings.filterwarnings("ignore")
    lf_pts, lf_missing = load_points([(r, t) for r in RPMS for t in LF_THETAS], lf_run)
    l9_pts, l9_missing = load_points(l9_points(), l9_run)
    hf_pts, hf_missing = load_points(hf_points(), hf_run)
    n = len(hf_pts)
    print(f"LF (L8) points {len(lf_pts)}/40 (missing: {lf_missing})")
    print(f"L9 points {len(l9_pts)}/{len(l9_points())} (missing: {l9_missing})")
    print(f"HF (L10) points {n}/{N_REQUIRED} (missing: {hf_missing})")
    hf_use = {k: v for k, v in hf_pts.items() if k in lf_pts and k in l9_pts}
    if len(hf_use) < n:
        print("WARNING: HF points without L8 and L9 counterpart dropped:", sorted(set(hf_pts) - set(hf_use)))
    res, stats = {}, {}
    if len(hf_use) >= 3 and len(lf_pts) >= 8 and len(l9_pts) >= 3:
        for q in QOIS:
            lf = train_lf2d(np.array(sorted(lf_pts, key=lambda k: (k[1], k[0]))),
                            [lf_pts[k][q] for k in sorted(lf_pts, key=lambda k: (k[1], k[0]))])
            s, rows = loo_h2(q, lf, l9_pts, lf_pts, hf_use)
            stats[q] = s
            res[q] = dict(stats=s, points=rows)
            print(f"\n{q}: n={s['n']}  rel RMSE  MF {s['mf_rmse']:.1%}  HF-only {s['hf_rmse']:.1%}  "
                  f"L8xratio {s['ratio_rmse']:.1%}  coverage MF {s['mf_cover']}/{s['n']}  "
                  f"HF-only {s['hf_cover']}/{s['n']}")
            print(f"  NON-decisive: L9 x ratio {s['l9_ratio_rmse']:.1%}  plain L9 {s['l9_rmse']:.1%}  "
                  f"(model 1 trained on {s['n_l9_model1']} L9 points)")
            print(f"{'rpm':>6}{'theta':>6}{'obs':>11}{'MF':>11}{'MF 95% interval':>26}{'HF-only':>11}"
                  f"{'L8xratio':>11}{'L9xratio*':>11}{'L9*':>11}")
            for r in rows:
                print(f"{r['rpm']:6g}{r['theta']:6g}{r['obs']:11.4g}{r['mf_pred']:11.4g}"
                      f"{'[%.4g, %.4g]' % (r['mf_lo'], r['mf_hi']):>26}{r['hf_pred']:11.4g}"
                      f"{r['ratio_pred']:11.4g}{r['l9ratio_pred']:11.4g}{r['l9_pred']:11.4g}")
        print("(* = non-decisive reference columns)")
        verdict = h2_verdict(stats)
    else:
        verdict = dict(verdict=f"PENDING ({len(hf_use)}/{N_REQUIRED} HF points)", failed=[])
    print(f"\nVERDICT: {verdict['verdict']}", verdict["failed"] or "")
    print(f"cost note: {len(l9_pts)} L9 runs used (model 1 training set); "
          f"{len(hf_use)} L10 runs scored.")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(dict(verdict=verdict, n_hf=len(hf_use), hf_missing=hf_missing, l9_missing=l9_missing,
                   lf_missing=lf_missing, n_l9_used=len(l9_pts), n_model1_train=len(l9_pts),
                   cost_note=f"{len(l9_pts)} L9 runs used", results=res), open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
