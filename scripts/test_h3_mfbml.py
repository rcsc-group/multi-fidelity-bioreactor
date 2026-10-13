"""H3 evaluation (pre-registered, diary.md): H2 with every GP stage on log(QoI).

Rule (quoted): "H3: H2 exactly (same models, data, LOO, decisive baselines, bar), except every GP stage
works on log(QoI): LF KRR on log L8, model 1 on log L9, model 2 and the HF-only baseline on log L10.
Point prediction = exp(posterior mean of log), interval = exp(2.5/97.5% log quantiles); RMSE in relative
terms as before; the L8 x ratio baseline is unchanged. PASS = the H2 bar on all 13 points AND, as an extra
condition because the transform was chosen after the 8-point H2 interim, rel RMSE <= 10% for tau_95 and
EDR on the 5 points no model has seen yet."
PENDING while fewer than 13 HF points exist (same as H1/H2). The L8 x ratio baseline and the
non-decisive L9 columns stay on the raw scale. Implementation: loo_h2(..., log=True) in test_h2_mfbml.

Usage: uv run python scripts/test_h3_mfbml.py
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

from scripts.test_h1_mfbml import (  # noqa: E402
    CRITERIA_QOIS, LF_THETAS, N_REQUIRED, QOIS, RMSE_MAX, RPMS, h1_verdict, hf_points, hf_run, lf_run,
    load_points, rel_rmse_ratio, train_lf2d,
)
from scripts.test_h2_mfbml import l9_points, l9_run, loo_h2  # noqa: E402

OUT = ROOT / "experiments" / "multifidelity" / "test_h3_mfbml.json"
# HF points no model had seen when the log transform was chosen (after the 8-point H2 interim)
UNSEEN = ((17.5, 4.0), (25.0, 4.0), (32.5, 4.0), (17.5, 9.0), (25.0, 9.0))
N_UNSEEN = 5


def h3_verdict(stats: dict, unseen: dict | None) -> dict:
    """h1_verdict(stats) (the H2 bar) AND unseen[q]['mf_rmse'] <= 10% with n == 5 for tau_95, EDR.
    PENDING whenever the H2 bar is PENDING (fewer than 13 HF points)."""
    v = h1_verdict(stats)
    if v["verdict"].startswith("PENDING"):
        return v
    failed = list(v["failed"])
    for q in CRITERIA_QOIS:
        u = unseen[q]
        if u["n"] != N_UNSEEN:
            failed.append((q, f"unseen n=={N_UNSEEN}"))
        if not u["mf_rmse"] <= RMSE_MAX:
            failed.append((q, "unseen rmse<=10%"))
    return dict(verdict="FAIL" if failed else "PASS", failed=failed)


def unseen_stats(rows) -> dict:
    """MF relative RMSE over the LOO rows whose (rpm, theta) is in UNSEEN."""
    sub = [r for r in rows if (float(r["rpm"]), float(r["theta"])) in UNSEEN]
    if not sub:
        return dict(n=0, mf_rmse=float("nan"))
    return dict(n=len(sub), mf_rmse=rel_rmse_ratio([r["mf_pred"] for r in sub], [r["obs"] for r in sub]))


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
    res, stats, unseen = {}, {}, {}
    if len(hf_use) >= 3 and len(lf_pts) >= 8 and len(l9_pts) >= 3:
        keys8 = sorted(lf_pts, key=lambda k: (k[1], k[0]))
        for q in QOIS:
            lf = train_lf2d(np.array(keys8), np.log([lf_pts[k][q] for k in keys8]))   # LF KRR on log L8
            s, rows = loo_h2(q, lf, l9_pts, lf_pts, hf_use, log=True)
            u = unseen_stats(rows)
            stats[q], unseen[q] = s, u
            res[q] = dict(stats=s, unseen=u, points=rows)
            print(f"\n{q}: n={s['n']}  rel RMSE  MF {s['mf_rmse']:.1%}  HF-only {s['hf_rmse']:.1%}  "
                  f"L8xratio {s['ratio_rmse']:.1%}  coverage MF {s['mf_cover']}/{s['n']}  "
                  f"HF-only {s['hf_cover']}/{s['n']}")
            print(f"  UNSEEN subset: MF rel RMSE {u['mf_rmse']:.1%} (n={u['n']})")
            print(f"  NON-decisive: L9 x ratio {s['l9_ratio_rmse']:.1%}  plain L9 {s['l9_rmse']:.1%}  "
                  f"(model 1 trained on {s['n_l9_model1']} L9 points)")
            print(f"{'rpm':>6}{'theta':>6}{'obs':>11}{'MF':>11}{'MF 95% interval':>26}{'HF-only':>11}"
                  f"{'L8xratio':>11}{'L9xratio*':>11}{'L9*':>11}")
            for r in rows:
                print(f"{r['rpm']:6g}{r['theta']:6g}{r['obs']:11.4g}{r['mf_pred']:11.4g}"
                      f"{'[%.4g, %.4g]' % (r['mf_lo'], r['mf_hi']):>26}{r['hf_pred']:11.4g}"
                      f"{r['ratio_pred']:11.4g}{r['l9ratio_pred']:11.4g}{r['l9_pred']:11.4g}")
        print("(* = non-decisive reference columns)")
        verdict = h3_verdict(stats, unseen)
    else:
        verdict = dict(verdict=f"PENDING ({len(hf_use)}/{N_REQUIRED} HF points)", failed=[])
    print(f"\nVERDICT: {verdict['verdict']}", verdict["failed"] or "")
    print(f"cost note: {len(l9_pts)} L9 runs used (model 1 training set); "
          f"{len(hf_use)} L10 runs scored.")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(dict(verdict=verdict, n_hf=len(hf_use), hf_missing=hf_missing, l9_missing=l9_missing,
                   lf_missing=lf_missing, n_l9_used=len(l9_pts), unseen=list(UNSEEN),
                   cost_note=f"{len(l9_pts)} L9 runs used", results=res), open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
