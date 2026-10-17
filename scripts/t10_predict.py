"""Test T10 predictions, written BEFORE any new L10 kLa data exist (diary 2026-10-17).

Ratio baseline: log kLa_1T_10 at L10 = log R + g8(rpm, theta) - g8(32.5, 7), with g8 a full quadratic in (rpm, theta) fitted by
least squares to all L8 runs that have kLa_1T_10, and log R = log of the L10 anchor at 32.5 rpm / theta 7 (mean of
l10c_rpm32.5_seg1 and v_l10_rpm32.5_th7, both release-after-spinup L10 runs).
Writes experiments/multifidelity/t10_predictions.json.
Usage: uv run python scripts/t10_predict.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
POINTS = [(r, 7.0) for r in (17.5, 25.0, 32.5)] + [(r, t) for t in (4.0, 9.0) for r in (17.5, 25.0, 32.5)]


def design(r, t):
    r, t = np.asarray(r, float), np.asarray(t, float)
    return np.column_stack([np.ones_like(r), r, t, r * r, r * t, t * t])


def main():
    rows = [x for x in csv.DictReader((ROOT / "experiments/multifidelity/kla_mix_dataset.csv").open()) if x["level"] == "8"]
    r = np.array([float(x["rpm"]) for x in rows])
    t = np.array([float(x["theta"]) for x in rows])
    y = np.array([float(x["kLa_1T_10"]) if x["kLa_1T_10"] not in ("", "nan") else np.nan for x in rows])
    ok = np.isfinite(y) & (y > 0)
    c, *_ = np.linalg.lstsq(design(r[ok], t[ok]), np.log(y[ok]), rcond=None)
    res = np.log(y[ok]) - design(r[ok], t[ok]) @ c
    s8 = float(res.std(ddof=6))
    anchors = [json.loads((ROOT / "runs" / n / "results.json").read_text())["kLa_1T_10"]
               for n in ("l10c_rpm32.5_seg1", "v_l10_rpm32.5_th7")]
    anchor = float(np.mean(anchors))
    g = lambda rr, tt: float((design([rr], [tt]) @ c)[0])
    pred = {f"{p[0]:g}/{p[1]:g}": float(np.exp(np.log(anchor) + g(*p) - g(32.5, 7.0))) for p in POINTS}
    out = dict(n_L8=int(ok.sum()), s8_residual_log=s8, anchor_L10=anchor, anchor_runs=anchors, coef=c.tolist(), prediction_L10=pred)
    (ROOT / "experiments/multifidelity/t10_predictions.json").write_text(json.dumps(out, indent=1))
    print(f"L8 runs {ok.sum()}, residual sd (log) {s8:.3f}, L10 anchor {anchor:.2f} (runs {anchors})")
    for k, v in pred.items():
        print(f"  rpm/theta {k}: predicted L10 kLa_1T_10 {v:.1f} 1/h")


if __name__ == "__main__":
    main()
