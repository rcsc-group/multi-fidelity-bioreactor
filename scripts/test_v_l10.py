"""Test V-L10 (pre-registered, diary 2026-10-12): D_eff/D_phys of the tracer at L10 from the variance budget.

Rule (quoted): "L10 NUMERICAL-DOMINATED if median r >= 10 in the chi 0.1-0.4 window; PHYSICAL if <= 2. Reported,
not decisive: the dx^1.5 hypothesis predicts r(L10)/r(L9) ~ 0.35 at equal chi (32.5 rpm L9 bins: 173 / 189 / 236
-> L10 ~ 60-85)." Run: v_l10_rpm32.5_th7 (restart from settled L10, tracer released 1 cycle later).
L8 and L9 (scal_l{8,9}_rpm32.5_th7, from test_v_variance_budget.json) are binned in the same chi window, so all
three levels are compared at equal chi.
Usage: uv run python scripts/test_v_l10.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.test_v_variance_budget import series  # noqa: E402

RUN = "v_l10_rpm32.5_th7"
REF = {8: "scal_l8_rpm32.5_th7", 9: "scal_l9_rpm32.5_th7"}  # binned from test_v_variance_budget.json
EXTRA = {6: "v_l6_rpm32.5_th7", 7: "v_l7_rpm32.5_th7"}     # same protocol as RUN; reported, not decisive
WINDOW = (0.1, 0.4)
BINS = [(0.1, 0.2), (0.2, 0.3), (0.3, 0.4)]
OUT = ROOT / "experiments/multifidelity/test_v_l10.json"


def binned(chi, r):
    chi, r = np.asarray(chi, float), np.asarray(r, float)
    ok = np.isfinite(r)
    med = lambda a, b: float(np.median(r[ok & (chi >= a) & (chi < b)])) if (ok & (chi >= a) & (chi < b)).any() \
        else float("nan")  # noqa: E731
    return dict(window=med(*WINDOW), bins=[med(a, b) for a, b in BINS],
                n_window=int((ok & (chi >= WINDOW[0]) & (chi < WINDOW[1])).sum()))


def verdict(m):
    return "NUMERICAL-DOMINATED" if m >= 10 else "PHYSICAL" if m <= 2 else "MIXED"


def main():
    ref = json.loads((ROOT / "experiments/multifidelity/test_v_variance_budget.json").read_text())["runs"]
    s = series(RUN)
    res = {10: dict(run=RUN, chi_max=float(np.nanmax(s["chi"])), **binned(s["chi"], s["r"]))}
    for L, run in EXTRA.items():
        if (ROOT / "runs" / run / "results.json").exists():
            e = series(run)
            res[L] = dict(run=run, chi_max=float(np.nanmax(e["chi"])), **binned(e["chi"], e["r"]))
    for L, run in REF.items():
        res[L] = dict(run=run, **binned(ref[run]["chi"], ref[run]["r"]))
    print(f"{'level':>5} {'run':24s} {'median r, chi 0.1-0.4':>22s}  bins " + "  ".join(f"{a}-{b}" for a, b in BINS))
    for L in sorted(res):
        d = res[L]
        print(f"{L:5d} {d['run']:24s} {d['window']:16.1f} (n={d['n_window']:3d})  "
              + "  ".join(f"{x:7.1f}" for x in d["bins"]))
    v = verdict(res[10]["window"])
    ratio = res[10]["window"] / res[9]["window"]
    print(f"L10 max chi reached {res[10]['chi_max']:.3f}; r(L10)/r(L9) in window {ratio:.2f} (dx^1.5 predicts 0.35)")
    print("VERDICT (L10):", v)
    OUT.write_text(json.dumps(dict(verdict=v, ratio_l10_l9=ratio, window=WINDOW, bins=BINS,
                                   levels={str(k): d for k, d in res.items()}), indent=1))


if __name__ == "__main__":
    main()
