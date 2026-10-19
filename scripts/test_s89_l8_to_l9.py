"""Test S89 (pre-registered, diary 2026-10-19): predict kLa_1T_10 and dtmix_0.50 at L9 (theta 7) from L8, leave-one-out over the L9 runs.

LF = L8 KRR on ln(QoI) over all L8 runs with the QoI; HF = L9 theta-7 runs (fig9_l9 old cold runs + e9_l9 same-protocol runs).
Compared: MF (FullBayesKRRLRGP), L9-only GP, L8 KRR x one constant ratio, plain L8 KRR. Metric: RMS of ln(pred/obs); 95% coverage.
USEFUL if MF RMS <= 0.25 and below the L9-only and ratio baselines; MARGINAL if below both but > 0.25; else NOT USEFUL.
Held-out unit = all L9 runs within 0.25 rpm of the scored run (first version, plain leave-one-out with a duplicate run, gave MF 0.143 / 0.119).
Usage: uv run python scripts/test_s89_l8_to_l9.py [--old-only]
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
from scripts.test_h1_mfbml import BOUNDS, train_lf2d  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

QS = ("kLa_1T_10", "dtmix_0.50")
GROUP = 0.25


def collect(patterns, q):
    pts = []
    for pat in patterns:
        for d in sorted((ROOT / "runs").glob(pat)):
            if "_ext" in d.name:   # continuation of the same run: a duplicate point (found 2026-10-19)
                continue
            f = d / "results.json"
            if not f.exists():
                continue
            v = json.loads(f.read_text()).get(q)
            if v is None or not (v > 0):
                continue
            p = json.loads((d / "params.json").read_text())
            pts.append((p["omega_b"] * 60 / (2 * math.pi), p["theta_max"][0], float(v), d.name))
    return pts


def main():
    warnings.filterwarnings("ignore")
    old_only = "--old-only" in sys.argv
    out = {}
    for q in QS:
        l8 = collect(["t6_l8_*", "kmix_l8_*", "t5_l8_*", "e8_l8_*", "w_l8_*", "dose_l8_*"], q)
        l8 = [p for p in l8 if abs(p[0] - round(p[0], 3)) < 1e-9 or True]
        X8 = np.array([[p[0], p[1]] for p in l8])
        lf = train_lf2d(X8, np.log([p[2] for p in l8]))
        l9 = collect(["fig9_l9_rpm*"] + ([] if old_only else ["e9_l9_rpm*_th7"]), q)
        l9 = [p for p in l9 if p[1] == 7.0]
        X9 = np.array([[p[0], p[1]] for p in l9])
        y9 = np.log([p[2] for p in l9])
        f8 = np.ravel(lf.predict(X9))            # L8 KRR (ln) at the L9 points
        rows = []
        for i in range(len(l9)):
            tr = np.abs(X9[:, 0] - X9[i, 0]) > GROUP   # leave the whole rpm group out: runs within 0.25 rpm are near-replicates (found 2026-10-19)
            xt = X9[i:i + 1]
            m = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(X9[tr], y9[tr])
            mean = float(m.predict_mean(xt)[0])
            lo, hi = m.quantile(xt, [0.025, 0.975])[0]
            h = FullBayesKRRLRGP(BOUNDS, lf_model=None, basis="ordinary").train(X9[tr], y9[tr])
            hm = float(h.predict_mean(xt)[0])
            ratio = f8[i] + float(np.mean(y9[tr] - f8[tr]))
            rows.append(dict(run=l9[i][3], rpm=l9[i][0], obs=float(y9[i]), mf=mean, lo=float(lo), hi=float(hi), hf=hm, ratio=ratio, l8=float(f8[i])))
        rms = lambda k: float(np.sqrt(np.mean([(r[k] - r["obs"]) ** 2 for r in rows])))   # noqa: E731
        cover = sum(r["lo"] <= r["obs"] <= r["hi"] for r in rows)
        mf, hf, ra, pl = rms("mf"), rms("hf"), rms("ratio"), rms("l8")
        verdict = "USEFUL" if mf <= 0.25 and mf < hf and mf < ra else "MARGINAL" if mf < hf and mf < ra else "NOT USEFUL"
        out[q] = dict(n_L8=len(l8), n_L9=len(l9), rms=dict(mf=mf, hf=hf, ratio=ra, plain_l8=pl), cover=cover, verdict=verdict, rows=rows)
        print(f"{q}: L8 runs {len(l8)}, L9 runs {len(l9)}; RMS ln error MF {mf:.3f}  L9-only {hf:.3f}  L8 x ratio {ra:.3f}  plain L8 {pl:.3f}; "
              f"coverage {cover}/{len(l9)} -> {verdict}")
        for r in sorted(rows, key=lambda r: r["rpm"]):
            print(f"   rpm {r['rpm']:6.2f} obs {math.exp(r['obs']):7.2f}  MF {math.exp(r['mf']):7.2f} [{math.exp(r['lo']):6.1f}, {math.exp(r['hi']):6.1f}]  L9-only {math.exp(r['hf']):7.2f}  L8 {math.exp(r['l8']):7.2f}")
    (ROOT / f"experiments/multifidelity/test_s89_l8_to_l9{'_oldonly' if old_only else ''}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
