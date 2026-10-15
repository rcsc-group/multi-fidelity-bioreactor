"""Test E0 (pre-registered, diary 2026-10-15): is a single-run L8 QoI (theta 7) = smooth trend in rpm + independent
scatter? Quadratic trend on the 21 runs of the 0.625-rpm grid 20-32.5 (T5 + kmix), residual sd s, lag-1
autocorrelation rho1, +-0.1 rpm pair |differences| vs 2 s / sqrt(pi). All QoIs in log.
NOISE-LIKE if |rho1| <= 0.43 and median pair |diff| in [0.5, 2] x 2 s / sqrt(pi); STRUCTURED if rho1 > 0.43; else MIXED.
Usage: uv run python scripts/test_e0_scatter.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois  # noqa: E402

PRIMARY = ("kLa_1T_25", "dtmix_0.95")
REPORTED = ("kLa_1T_10", "dtmix_0.50")
CONTROL = ("tau_95_t", "edr_mean_t")
PAIRS = ((22.4, 22.6), (29.9, 30.1))
RHO_BAR = 1.96 / np.sqrt(21)
L10_CORE_H = 935.0


def value(run, q):
    if q.endswith("_t"):
        h = qois(run)
        return h[q] if h else np.nan
    f = ROOT / "runs" / run / "results.json"
    return float(json.loads(f.read_text()).get(q, np.nan)) if f.exists() else np.nan


def pair_run(r):
    return f"t5_l8_rpm{r:g}_th7"


def analyse(rpm, y, pair_diffs, deg):
    c = np.polyfit(rpm, y, deg)
    e = y - np.polyval(c, rpm)
    s = float(np.std(e, ddof=deg + 1))
    rho1 = float(np.corrcoef(e[:-1], e[1:])[0, 1])
    expect = 2 * s / np.sqrt(np.pi)
    med = float(np.median(pair_diffs))
    if rho1 > RHO_BAR:
        reading = "STRUCTURED"
    elif abs(rho1) <= RHO_BAR and 0.5 * expect <= med <= 2 * expect:
        reading = "NOISE-LIKE"
    else:
        reading = "MIXED"
    return dict(s=s, rho1=rho1, pair_med=med, pair_expect=expect, reading=reading,
                n_for_5pct=(s / 0.05) ** 2, residuals=e.tolist())


def main():
    ids = json.loads((ROOT / "experiments/multifidelity/test_t5_roughness.json").read_text())["run_ids"]
    grid = sorted(float(r) for r in ids if abs(float(r) / 0.625 - round(float(r) / 0.625)) < 1e-9
                  and 20 <= float(r) <= 32.5)
    assert len(grid) == 21, grid
    out = {}
    print(f"rho1 bar {RHO_BAR:.2f}; n grid runs {len(grid)}")
    for q in PRIMARY + REPORTED + CONTROL:
        y = np.log([value(ids[f"{r:g}"], q) for r in grid])
        pd = [abs(np.log(value(pair_run(b), q)) - np.log(value(pair_run(a), q))) for a, b in PAIRS]
        if not (np.all(np.isfinite(y)) and np.all(np.isfinite(pd))):
            print(f"{q}: missing values, skipped")
            continue
        res = {deg: analyse(np.array(grid), y, pd, deg) for deg in (2, 3)}
        out[q] = res
        a = res[2]
        tag = "PRIMARY" if q in PRIMARY else "reported" if q in REPORTED else "control"
        print(f"{q:12s} [{tag}] s {a['s']:.3f}  rho1 {a['rho1']:+.2f}  pair |d| med {a['pair_med']:.3f} "
              f"(iid expect {a['pair_expect']:.3f})  -> {a['reading']}   | cubic: rho1 {res[3]['rho1']:+.2f} "
              f"{res[3]['reading']}  | runs for +-5%: {a['n_for_5pct']:.0f} "
              f"(L10 ~{a['n_for_5pct'] * L10_CORE_H:,.0f} core-h per design point)")
    (ROOT / "experiments/multifidelity/test_e0_scatter.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
