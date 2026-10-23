"""Check H3 interval width (diary 2026-10-19, pre-registered): variants A (published), B (L9 variance propagated), C (B and the L9 run at the held-out point removed).
Usage: uv run python scripts/h3_band_check.py
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from scripts.test_h1_mfbml import BOUNDS, LF_THETAS, RPMS, hf_points, hf_run, lf_run, load_points, train_lf2d  # noqa: E402
from scripts.test_h2_mfbml import MeanAdapter, l9_points, l9_run  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

QS = tuple(sys.argv[1].split(",")) if len(sys.argv) > 1 else ("tau_95_t", "tau_mean_t", "edr_mean_t")   # e.g. tau_amp_t,edr_max_t (Kim statistics)
DELTA = 0.02


class Offset(MeanAdapter):
    off = 0.0

    def predict(self, X):
        return super().predict(X) + self.off


def sd_log(m, X):
    mu, sc, nu, w = m.predict_dist(X)
    return np.sqrt(w @ (sc ** 2) + w @ (mu ** 2) - (w @ mu) ** 2)


def main():
    warnings.filterwarnings("ignore")
    lf_pts, _ = load_points([(r, t) for r in RPMS for t in LF_THETAS], lf_run)
    l9_pts, _ = load_points(l9_points(), l9_run)
    hf_pts, _ = load_points(hf_points(), hf_run)
    hk = sorted(hf_pts, key=lambda k: (k[1], k[0]))
    X = np.array(hk)
    out = {}
    for q in QS:
        k8 = sorted(lf_pts, key=lambda k: (k[1], k[0]))
        lf = train_lf2d(np.array(k8), np.log([lf_pts[k][q] for k in k8]))
        yh = np.log([hf_pts[k][q] for k in hk])
        res = {v: [] for v in "ABC"}
        rhos = []
        for i, k in enumerate(hk):
            tr = np.arange(len(hk)) != i
            xt = X[i:i + 1]
            for v in "ABC":
                k9 = [p for p in sorted(l9_pts, key=lambda p: (p[1], p[0])) if not (v == "C" and p == k)]
                m1 = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(np.array(k9), np.log([l9_pts[p][q] for p in k9]))
                ad = Offset(m1)
                m2 = FullBayesKRRLRGP(BOUNDS, lf_model=ad, basis="linear").train(X[tr], yh[tr])
                mean = float(m2.predict_mean(xt)[0])
                lo, hi = m2.quantile(xt, [0.025, 0.975])[0]
                if v != "A":
                    ad.off = DELTA
                    rho = (float(m2.predict_mean(xt)[0]) - mean) / DELTA
                    ad.off = 0.0
                    hw = np.sqrt(((hi - lo) / 2) ** 2 + (1.96 * rho * sd_log(m1, xt)[0]) ** 2)
                    lo, hi = mean - hw, mean + hw
                    rhos.append(rho)
                res[v].append(dict(obs=float(yh[i]), mean=mean, lo=float(lo), hi=float(hi)))
        row = {}
        for v in "ABC":
            r = res[v]
            err = np.array([np.exp(x["mean"]) / np.exp(x["obs"]) - 1 for x in r])
            row[v] = dict(rmse=float(np.sqrt(np.mean(err ** 2))), cover=int(sum(x["lo"] <= x["obs"] <= x["hi"] for x in r)),
                          half_width_pct=float(100 * np.mean([(x["hi"] - x["lo"]) / 2 for x in r])))
        row["rho_mean"] = float(np.mean(rhos))
        row["points_B"] = [dict(rpm=k[0], theta=k[1], **res["B"][i]) for i, k in enumerate(hk)]
        out[q] = row
        print(q, " | ".join(f"{v}: RMSE {row[v]['rmse']:.1%} cover {row[v]['cover']}/13 half-width {row[v]['half_width_pct']:.1f}%" for v in "ABC"), f"| rho {row['rho_mean']:.2f}")
    cover_c = [out[q]["C"]["cover"] for q in QS]
    ratio = [out[q]["C"]["rmse"] / out[q]["A"]["rmse"] for q in QS]
    verdict = "TOO NARROW" if min(cover_c) <= 10 or max(ratio) > 2 else "HONEST" if min(cover_c) >= 11 and max(ratio) <= 1.5 else "MIXED"
    out["verdict"] = verdict
    print("verdict:", verdict, "coverage C", cover_c, "RMSE ratio C/A", [round(r, 2) for r in ratio])
    (ROOT / "experiments/multifidelity" / ("h3_band_check.json" if len(sys.argv) < 2 else "h3_band_check_kim.json")).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
