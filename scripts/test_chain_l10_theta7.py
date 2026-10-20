"""Chain test L8 -> L9 -> L10 for kLa_1T_10 at theta 7 (pre-registered, diary 2026-10-20), leave-one-out over the L10 runs.

L10 runs: t10_l10_rpm{17.5, 20, 22.5, 25, 27.5, 30}_th7 + v_l10_rpm32.5_th7. L9: fig9_l9 runs (+ e9 22.4, 22.6, 29.9, 30.1; the e9 runs at 17.5/25/32.5 are
bit-near copies of fig9 and are excluded). Variant A: L9 model sees all L9 runs. Variant C: the L9 runs within 0.25 rpm of the held-out rpm are removed too.
Bands: L9 variance propagated (scripts/h3_band_check.py recipe). Baselines: L9 x constant ratio, L10-only GP, L8 x constant ratio.
USEFUL if variant A RMS ln error <= 0.20 and below the L9 x ratio baseline; reported: variant C, coverage.
Usage: uv run python scripts/test_chain_l10_theta7.py
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
from scripts.test_h1_mfbml import BOUNDS, train_lf2d  # noqa: E402
from scripts.test_s89_l8_to_l9 import collect  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

Q = "kLa_1T_10"
L10 = ["t10_l10_rpm17.5_th7", "t10_l10_rpm20_th7", "t10_l10_rpm22.5_th7", "t10_l10_rpm25_th7", "t10_l10_rpm27.5_th7", "t10_l10_rpm30_th7", "v_l10_rpm32.5_th7"]


def main():
    warnings.filterwarnings("ignore")
    l8 = collect(["t6_l8_*", "kmix_l8_*", "t5_l8_*", "e8_l8_*", "w_l8_*", "dose_l8_*"], Q)
    lf = train_lf2d(np.array([[p[0], p[1]] for p in l8]), np.log([p[2] for p in l8]))
    l9 = [p for p in collect(["fig9_l9_rpm*", "e9_l9_rpm*_th7"], Q) if p[1] == 7.0 and p[3] not in ("e9_l9_rpm17.5_th7", "e9_l9_rpm25_th7", "e9_l9_rpm32.5_th7")]
    X9 = np.array([[p[0], p[1]] for p in l9])
    y9 = np.log([p[2] for p in l9])
    l10 = []
    for r in L10:
        f = ROOT / "runs" / r / "results.json"
        if f.exists():
            v = json.loads(f.read_text()).get(Q)
            p = json.loads((ROOT / "runs" / r / "params.json").read_text())
            if v and v > 0:
                l10.append((p["omega_b"] * 60 / (2 * math.pi), 7.0, float(v), r))
    X10 = np.array([[p[0], p[1]] for p in l10])
    y10 = np.log([p[2] for p in l10])
    n = len(l10)
    f9 = np.ravel(lf.predict(X10))
    res = {v: [] for v in "AC"}
    base = []
    for i in range(n):
        tr = np.arange(n) != i
        xt = X10[i:i + 1]
        for v in "AC":
            keep = np.abs(X9[:, 0] - X10[i, 0]) > 0.25 if v == "C" else np.ones(len(X9), bool)
            m1 = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(X9[keep], y9[keep])
            ad = Offset(m1)
            m2 = FullBayesKRRLRGP(BOUNDS, lf_model=ad, basis="linear").train(X10[tr], y10[tr])
            mean = float(m2.predict_mean(xt)[0])
            lo, hi = m2.quantile(xt, [0.025, 0.975])[0]
            ad.off = DELTA
            rho = (float(m2.predict_mean(xt)[0]) - mean) / DELTA
            ad.off = 0.0
            hw = math.sqrt(((hi - lo) / 2) ** 2 + (1.96 * rho * sd_log(m1, xt)[0]) ** 2)
            res[v].append(dict(rpm=float(X10[i, 0]), obs=float(y10[i]), mean=mean, lo=mean - hw, hi=mean + hw))
        m1a = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(X9, y9)
        f9i = float(m1a.predict_mean(xt)[0])
        f9tr = np.ravel(m1a.predict_mean(X10[tr]))
        ratio9 = f9i + float(np.mean(y10[tr] - f9tr))
        ratio8 = f9[i] + float(np.mean(y10[tr] - f9[tr]))
        h = FullBayesKRRLRGP(BOUNDS, lf_model=None, basis="ordinary").train(X10[tr], y10[tr])
        base.append(dict(l9ratio=ratio9, l8ratio=ratio8, l10only=float(h.predict_mean(xt)[0])))
    rms = lambda pred: float(np.sqrt(np.mean([(p - o) ** 2 for p, o in zip(pred, y10)])))   # noqa: E731
    out = dict(n_L10=n, n_L9=len(l9), n_L8=len(l8))
    for v in "AC":
        out[v] = dict(rms=rms([r["mean"] for r in res[v]]), cover=int(sum(r["lo"] <= r["obs"] <= r["hi"] for r in res[v])),
                      half_width=float(np.mean([(r["hi"] - r["lo"]) / 2 for r in res[v]])), rows=res[v])
    out["baselines"] = dict(l9_ratio=rms([b["l9ratio"] for b in base]), l8_ratio=rms([b["l8ratio"] for b in base]), l10_only=rms([b["l10only"] for b in base]))
    verdict = "USEFUL" if out["A"]["rms"] <= 0.20 and out["A"]["rms"] < out["baselines"]["l9_ratio"] else "NOT USEFUL"
    out["verdict"] = verdict
    print(f"L10 runs {n}, L9 {len(l9)}, L8 {len(l8)}")
    print(f"A (L9 at the point available): RMS ln {out['A']['rms']:.3f}, coverage {out['A']['cover']}/{n}, half-width {100 * out['A']['half_width']:.0f}%")
    print(f"C (L9 nearby removed):         RMS ln {out['C']['rms']:.3f}, coverage {out['C']['cover']}/{n}, half-width {100 * out['C']['half_width']:.0f}%")
    print("baselines RMS ln:", {k: round(v, 3) for k, v in out["baselines"].items()}, "->", verdict)
    for r in res["A"]:
        print(f"   rpm {r['rpm']:5.1f} obs {math.exp(r['obs']):6.2f}  chain {math.exp(r['mean']):6.2f} [{math.exp(r['lo']):5.1f}, {math.exp(r['hi']):5.1f}]")
    (ROOT / "experiments/multifidelity/test_chain_l10_theta7.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
