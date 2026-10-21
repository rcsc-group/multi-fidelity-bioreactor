"""Hydrodynamic mfbml true test set (pre-registered, diary 2026-10-19): 6 NEW L10 runs at theta 5.5 (rpm 17.5, 20, 22.5, 27.5, 30, 32.5).

TRAIN = the 13 L10 runs (test_h1_mfbml.hf_points), L8 (40 runs, incl. theta 5.5) and L9 (16 + the new theta-5.5 warm starts).
TEST = the new L10 runs, no refit. Chain L8 KRR -> L9 GP -> L10 GP on ln(QoI), L9 variance propagated (h3_band_check recipe).
PASS if, for tau_95 and EDR, relative RMSE <= 10% and the 95% interval covers >= 5 of 6 (all available runs if fewer exist).
Variant C: the L9 runs at theta 5.5 are removed. Baselines: L9 x ratio, L8 x ratio (ratios from the 13 training runs).
Usage: uv run python scripts/test_h4_testset.py
"""
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from scripts.h3_band_check import DELTA, Offset, sd_log  # noqa: E402
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from scripts.test_h1_mfbml import BOUNDS, LF_THETAS, RPMS, hf_points, hf_run, lf_run, load_points, train_lf2d  # noqa: E402
from scripts.test_h2_mfbml import l9_points, l9_run  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

QS = ("tau_95_t", "tau_mean_t", "edr_mean_t")
TEST_RPM = (17.5, 20, 22.5, 27.5, 30, 32.5)
TH = 5.5


def main():
    warnings.filterwarnings("ignore")
    lf_pts, _ = load_points([(r, t) for r in RPMS for t in LF_THETAS], lf_run)
    l9_pts, _ = load_points(l9_points(), l9_run)
    hf_pts, _ = load_points(hf_points(), hf_run)
    test = {}
    for r in TEST_RPM:
        h = qois(f"x8_l10_rpm{r:g}_th{TH:g}") if (ROOT / "runs" / f"x8_l10_rpm{r:g}_th{TH:g}" / "results.json").exists() else None
        q9 = qois(f"x8_l9_rpm{r:g}_th{TH:g}") if (ROOT / "runs" / f"x8_l9_rpm{r:g}_th{TH:g}" / "results.json").exists() else None
        if h and q9:
            test[r] = (h, q9)
    print(f"test points available: {sorted(test)} of {len(TEST_RPM)}")
    out = {}
    hk = sorted(hf_pts, key=lambda k: (k[1], k[0]))
    X = np.array(hk)
    Xt = np.array([[r, TH] for r in sorted(test)])
    for q in QS:
        k8 = sorted(lf_pts, key=lambda k: (k[1], k[0]))
        lf = train_lf2d(np.array(k8), np.log([lf_pts[k][q] for k in k8]))
        l9 = {k: v[q] for k, v in l9_pts.items()}
        l9_new = {(r, TH): test[r][1][q] for r in test}
        yh = np.log([hf_pts[k][q] for k in hk])
        obs = np.array([test[r][0][q] for r in sorted(test)])
        row = {}
        for v in "AC":
            d9 = dict(l9) if v == "C" else {**l9, **l9_new}
            k9 = sorted(d9, key=lambda k: (k[1], k[0]))
            m1 = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(np.array(k9), np.log([d9[k] for k in k9]))
            ad = Offset(m1)
            m2 = FullBayesKRRLRGP(BOUNDS, lf_model=ad, basis="linear").train(X, yh)
            mean = np.ravel(m2.predict_mean(Xt))
            qq = m2.quantile(Xt, [0.025, 0.975])
            ad.off = DELTA
            rho = (np.ravel(m2.predict_mean(Xt)) - mean) / DELTA
            ad.off = 0.0
            hw = np.sqrt(((qq[:, 1] - qq[:, 0]) / 2) ** 2 + (1.96 * rho * sd_log(m1, Xt)) ** 2)
            pred = np.exp(mean)
            row[v] = dict(rmse=float(np.sqrt(np.mean((pred / obs - 1) ** 2))), cover=int(np.sum((mean - hw <= np.log(obs)) & (np.log(obs) <= mean + hw))),
                          half_width_pct=float(100 * hw.mean()), pred=pred.tolist())
        # ratio baselines from the 13 training runs
        y9 = np.array([l9[k] for k in hk]); y8 = np.array([lf_pts[k][q] for k in hk]); y10 = np.exp(yh)
        p9 = np.exp(np.ravel(m1.predict_mean(Xt))) if False else np.array([l9_new[(r, TH)] for r in sorted(test)])
        p8 = np.array([lf_pts[(float(r), TH)][q] if (float(r), TH) in lf_pts else lf_pts[(r, TH)][q] for r in sorted(test)])
        r9, r8 = np.exp(np.mean(np.log(y10 / y9))), np.exp(np.mean(np.log(y10 / y8)))
        row["l9_ratio_rmse"] = float(np.sqrt(np.mean((p9 * r9 / obs - 1) ** 2)))
        row["l8_ratio_rmse"] = float(np.sqrt(np.mean((p8 * r8 / obs - 1) ** 2)))
        out[q] = row
        print(f"{q}: A (with L9 at the point) RMSE {row['A']['rmse']:.1%} cover {row['A']['cover']}/{len(test)} half-width {row['A']['half_width_pct']:.0f}%"
              f" | C (without) RMSE {row['C']['rmse']:.1%} cover {row['C']['cover']}/{len(test)} | L9 x ratio {row['l9_ratio_rmse']:.1%}  L8 x ratio {row['l8_ratio_rmse']:.1%}")
    n = len(test)
    need = 5 if n >= 6 else max(1, math.ceil(n * 5 / 6))
    ok = all(out[q]["A"]["rmse"] <= 0.10 and out[q]["A"]["cover"] >= need for q in ("tau_95_t", "edr_mean_t"))
    out["verdict"] = ("PASS" if ok else "FAIL") + ("" if n >= 6 else f" (interim, {n} of 6 test points)")
    print("verdict:", out["verdict"])
    (ROOT / "experiments/multifidelity/test_h4_testset.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
