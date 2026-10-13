"""Zero-compute screens B (L9 middle level) and D (dx-continuous extrapolation to L10), pre-registered in diary
2026-10-13. Data: experiments/multifidelity/hydro_dataset_v2.json (theta 7, L6-L10, 9 rpm 17.5-37.5).
Usage: uv run python scripts/screen_bd.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
QS = ("tau_95_t", "edr_mean_t", "tau_mean_t")
RPMS = [17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]


def three_point(dx, q):
    """Exact q = q_inf + C dx^p through 3 points (dx ratio 2). Returns (q_inf, C, p) or None if not monotone."""
    d1, d2 = q[0] - q[1], q[1] - q[2]  # coarse-mid, mid-fine
    if d1 == 0 or d2 == 0 or np.sign(d1) != np.sign(d2):
        return None
    p = np.log(d1 / d2) / np.log(2.0)
    if p <= 0:
        return None
    C = d2 / (dx[1] ** p - dx[2] ** p)
    return q[2] - C * dx[2] ** p, C, p


def main():
    data = json.loads((ROOT / "experiments/multifidelity/hydro_dataset_v2.json").read_text())
    v = {(int(x["level"]), float(x["rpm"])): x for x in data}
    res = {"B": {}, "D": {}}
    print("Screen B: RMS residual of linear-in-rpm fit of L10/Lk")
    for q in QS:
        rr = {}
        for k in (8, 9):
            y = np.array([v[(10, r)][q] / v[(k, r)][q] for r in RPMS])
            rr[k] = float(np.std(y - np.polyval(np.polyfit(RPMS, y, 1), RPMS)))
            res["B"].setdefault(q, {})[f"ratio_L10_L{k}"] = y.tolist()
        res["B"][q].update(r8=rr[8], r9=rr[9])
        print(f"  {q}: r_8 {rr[8]:.4f}  r_9 {rr[9]:.4f}  r_9/r_8 {rr[9] / rr[8]:.2f}")
    b_ok = all(res["B"][q]["r9"] <= 0.5 * res["B"][q]["r8"] for q in ("tau_95_t", "edr_mean_t"))
    res["B"]["verdict"] = "WORTH BUILDING" if b_ok else "NOT SUPPORTED"
    print("  Screen B reading:", res["B"]["verdict"])

    print("Screen D: 3-point fit through L7, L8, L9 -> predict L10")
    dx = np.array([0.25 / 2 ** L for L in (7, 8, 9)])
    dx10 = 0.25 / 2 ** 10
    for q in QS:
        errs, pers, ps = [], [], []
        for r in RPMS:
            qq = np.array([v[(L, r)][q] for L in (7, 8, 9)])
            q10 = v[(10, r)][q]
            fit = three_point(dx, qq)
            pers.append(abs(qq[2] / q10 - 1))
            if fit is None:
                errs.append(np.nan)
                ps.append(np.nan)
                continue
            qi, C, p = fit
            errs.append(abs((qi + C * dx10 ** p) / q10 - 1))
            ps.append(p)
        errs, pers = np.array(errs), np.array(pers)
        n_ok = int(np.isfinite(errs).sum())
        med = float(np.nanmedian(np.where(np.isfinite(errs), errs, np.inf))) if n_ok else np.inf
        res["D"][q] = dict(err=errs.tolist(), persistence=pers.tolist(), p=ps, n_monotone=n_ok,
                           median_err_nonmonotone_as_inf=med, median_persistence=float(np.median(pers)))
        print(f"  {q}: monotone at {n_ok}/9 rpm; p = " + ", ".join(f"{x:.2f}" if np.isfinite(x) else "-" for x in ps))
        print(f"     |pred/L10-1|: " + ", ".join(f"{e:.3f}" if np.isfinite(e) else "n/a" for e in errs)
              + f"  median (n/a = inf) {med:.3f}; persistence |L9/L10-1| median {np.median(pers):.3f}")
    d_ok = all(res["D"][q]["median_err_nonmonotone_as_inf"] <= 0.05 and
               res["D"][q]["median_err_nonmonotone_as_inf"] < res["D"][q]["median_persistence"]
               for q in ("tau_95_t", "edr_mean_t"))
    res["D"]["verdict"] = "WORTH BUILDING" if d_ok else "NOT SUPPORTED"
    print("  Screen D reading:", res["D"]["verdict"])
    (ROOT / "experiments/multifidelity/screen_bd.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
