"""Test T10 report (pre-registered, diary 2026-10-17): L10 kLa_1T_10 map vs the frozen L8-based predictions (experiments/multifidelity/t10_predictions.json).
T10-A: RMS of ln(observed / predicted) over the new points: PASS <= 0.25, FAIL >= 0.40, else MARGINAL. Reported per point; also kLa_1T_25 / dtmix_0.50 where reached.
Usage: uv run python scripts/test_t10_report.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RUNS = {"17.5/7": "t10_l10_rpm17.5_th7", "25/7": "t10_l10_rpm25_th7", "17.5/4": "t10_l10_rpm17.5_th4", "25/4": "t10_l10_rpm25_th4",
        "32.5/4": "t10_l10_rpm32.5_th4", "17.5/9": "t10_l10_rpm17.5_th9", "25/9": "t10_l10_rpm25_th9", "32.5/9": "t10_l10_rpm32.5_th9"}


def val(run, q):
    f = ROOT / "runs" / run / "results.json"
    if not f.exists():
        return np.nan
    v = json.loads(f.read_text()).get(q, np.nan)
    return float(v) if v is not None and v > 0 else np.nan


def main():
    pred = json.loads((ROOT / "experiments/multifidelity/t10_predictions.json").read_text())["prediction_L10"]
    err = []
    print(f"{'rpm/theta':>9} {'predicted':>9} {'observed':>9} {'ln ratio':>9}   kLa_25  dtmix_0.50  dtmix_0.95")
    for k, run in RUNS.items():
        o = val(run, "kLa_1T_10")
        e = np.log(o / pred[k]) if np.isfinite(o) else np.nan
        err.append(e)
        print(f"{k:>9} {pred[k]:9.1f} {o:9.2f} {e:+9.3f}   {val(run, 'kLa_1T_25'):6.1f}  {val(run, 'dtmix_0.50'):9.1f}  {val(run, 'dtmix_0.95'):9.1f}")
    e = np.array(err)
    ok = np.isfinite(e)
    rms = float(np.sqrt(np.mean(e[ok] ** 2)))
    verdict = "PASS" if rms <= 0.25 else "FAIL" if rms >= 0.40 else "MARGINAL"
    print(f"T10-A: n {ok.sum()} of 8, RMS ln error {rms:.3f} -> {verdict} (mean ln error {e[ok].mean():+.3f})")
    (ROOT / "experiments/multifidelity/test_t10_report.json").write_text(json.dumps(dict(n=int(ok.sum()), rms=rms, verdict=verdict, ln_ratio=e.tolist()), indent=1))


if __name__ == "__main__":
    main()
