"""What do the cheap levels tell us about kLa(C* = 0.10) at L10? Zero compute.
(1) Rank agreement: Spearman correlation of L8 (and L9 at theta 7) with L10 over the L10 points; does L8 find the best L10 condition?
(2) Learning curve, two levels: L8 KRR -> L10 GP (FullBayesKRRLRGP) vs L10-only GP vs L8 x constant ratio; n L10 training runs, 15 random draws, error = RMS of ln(pred/obs) on held-out L10 runs.
Usage: uv run python scripts/kla_cheap_levels.py
"""
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from scripts.test_h1_mfbml import BOUNDS, train_lf2d  # noqa: E402
from scripts.test_s89_l8_to_l9 import collect  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

Q = "kLa_1T_10"


def l10_points():
    pts = []
    for d in sorted((ROOT / "runs").glob("t10_l10_*")) + [ROOT / "runs" / "v_l10_rpm32.5_th7"]:
        f = d / "results.json"
        if not f.exists():
            continue
        v = json.loads(f.read_text()).get(Q)
        if v and v > 0:
            p = json.loads((d / "params.json").read_text())
            pts.append((round(p["omega_b"] * 60 / (2 * math.pi), 3), p["theta_max"][0], float(v)))
    return pts


def main():
    warnings.filterwarnings("ignore")
    l8 = collect(["t6_l8_*", "kmix_l8_*", "t5_l8_*", "e8_l8_*", "w_l8_*", "dose_l8_*"], Q)
    lf = train_lf2d(np.array([[p[0], p[1]] for p in l8]), np.log([p[2] for p in l8]))
    l9 = {(round(p[0], 2), p[1]): p[2] for p in collect(["fig9_l9_rpm*"], Q) if "_ext" not in p[3]}
    pts = l10_points()
    X = np.array([[p[0], p[1]] for p in pts])
    y = np.log([p[2] for p in pts])
    f8 = np.ravel(lf.predict(X))
    print(f"L10 kLa points: {len(pts)}")
    rho, pv = spearmanr(f8, y)
    print(f"(1) Spearman L8 vs L10 over {len(pts)} points: {rho:.2f} (p {pv:.3f}); best L10 point {pts[int(np.argmax(y))][:2]}, best L8-fit point {pts[int(np.argmax(f8))][:2]}")
    sel = [(math.exp(y[i]), l9[(round(X[i, 0], 2), 7.0)]) for i in range(len(pts)) if X[i, 1] == 7.0 and (round(X[i, 0], 2), 7.0) in l9]
    r9, _ = spearmanr([a for a, b in sel], [b for a, b in sel])
    print(f"    Spearman L9 vs L10 at theta 7 ({len(sel)} points): {r9:.2f}")
    rng = np.random.default_rng(1)
    out = {}
    for n in (3, 4, 5, 6, 8, 10, 12):
        if n >= len(pts):
            continue
        acc = {"chain_L8_L10": [], "L10_only": [], "L8_ratio": []}
        for _ in range(15):
            tr = np.zeros(len(pts), bool)
            tr[rng.choice(len(pts), n, replace=False)] = True
            te = ~tr
            m = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(X[tr], y[tr])
            h = FullBayesKRRLRGP(BOUNDS, lf_model=None, basis="ordinary").train(X[tr], y[tr])
            acc["chain_L8_L10"].append(np.sqrt(np.mean((np.ravel(m.predict_mean(X[te])) - y[te]) ** 2)))
            acc["L10_only"].append(np.sqrt(np.mean((np.ravel(h.predict_mean(X[te])) - y[te]) ** 2)))
            acc["L8_ratio"].append(np.sqrt(np.mean((f8[te] + np.mean(y[tr] - f8[tr]) - y[te]) ** 2)))
        out[n] = {k: float(np.mean(v)) for k, v in acc.items()}
        print(f"(2) n = {n:2d} L10 runs: " + "  ".join(f"{k} {v:.2f}" for k, v in out[n].items()))
    (ROOT / "experiments/multifidelity/kla_cheap_levels.json").write_text(json.dumps(dict(spearman_L8=float(rho), spearman_L9=float(r9), curve=out, n_points=len(pts)), indent=1))


if __name__ == "__main__":
    main()
