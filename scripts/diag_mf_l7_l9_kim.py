"""H2: does an L7 sweep + 3 L9 points predict L9 on Kim's six numbers?

Written before the L7 data existed (2026-10-01), so the criteria are fixed in
advance. Targets, Kim's: kLa_1T_{10,25,50} (whole-period kLa, h^-1) and
dtmix_{0.50,0.75,0.95} (s), on the rpm sweep at 7 deg.
  LF = kmix_l7_rpm*  (L7, ~64x cheaper per cycle than L9)
  HF = fig9_l9_rpm*  (L9)
Protocol (same as diag_mf_early_late): 3 HF training points, both endpoints
plus one interior rpm, every choice; the rest are held out. Model: the Fig 13 MF
(LOO-KRR on LF + heteroscedastic universal kriging, REML), HF sigma = 5% of y.
Baselines: GP on the 3 HF points only; LF times the mean HF/LF ratio.
PASS for a target: MF mean held-out relative RMSE below BOTH baselines.
H2 holds if it passes for most of the six targets.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.diag_mf_early_late import one_case       # noqa: E402

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
TARGETS = ["kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "dtmix_0.50", "dtmix_0.75", "dtmix_0.95"]


def val(run, key):
    f = ROOT / "runs" / run / "results.json"
    return json.load(open(f)).get(key, math.nan) if f.exists() else math.nan


def main() -> None:
    verdicts = {}
    for key in TARGETS:
        rows = [(r, val(f"kmix_l7_rpm{r:g}", key), val(f"fig9_l9_rpm{r:g}", key)) for r in RPMS]
        a = np.array([q for q in rows if all(np.isfinite(q))], float)
        if len(a) < 6:
            print(f"\n== {key}: only {len(a)} rpm with both levels")
            continue
        x, lf, hf = a[:, 0], a[:, 1], a[:, 2]
        print(f"\n== {key}   n={len(x)}  corr(log L7, log L9) = "
              f"{np.corrcoef(np.log(lf), np.log(hf))[0, 1]:+.2f}  L9/L7 ratio "
              f"{(hf / lf).min():.2f}-{(hf / lf).max():.2f}")
        res = {}
        for mid in x[1:-1]:
            train = np.isin(x, [x[0], mid, x[-1]])
            for k, v in one_case(x, lf, hf, train).items():
                res.setdefault(k, []).append(v[0])
        means = {k: float(np.mean(v)) for k, v in res.items()}
        for k, v in means.items():
            print(f"  {k:11s} mean held-out rel. RMSE {v:.0%}")
        verdicts[key] = means["MF"] < min(means["HF-only GP"], means["LF x ratio"])
        print(f"  -> {'PASS' if verdicts[key] else 'FAIL'}")
    print(f"\nH2: {sum(verdicts.values())}/{len(verdicts)} targets pass")


if __name__ == "__main__":
    main()
