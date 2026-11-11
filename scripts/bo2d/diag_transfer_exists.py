"""Does ANY transfer function T exist with L9 = T(L8)? (2026-11-11, user: rank correlation misses non-monotone T.)
Shape-free test: if L9 is a single-valued function of L8, sites with near-equal L8 have near-equal L9 (within noise).
Leave-one-out nearest-neighbour in L8: predict ln L9 at site i by the ln L9 of the site whose ln L8 is closest
(k = 1 and the mean of k = 2). Compare with the mean-only predictor and with the L9 noise floor (campaign.noise_floor).
If T exists and the L8 sampling is dense enough, NN error -> noise floor. Also prints the closest-L8 pairs.
Sets: pair12 (pr8_l8 vs bo1_l9, same site), theta = 7 line (kmix L8 vs fig9 L9, one family each). kLa_1T_10 and tau_95.
  uv run python scripts/bo2d/diag_transfer_exists.py
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts" / "campaign2d")]
import campaign as C  # noqa: E402
from scripts.bo2d.pair_l8 import SITES  # noqa: E402

KEY = {"kLa_1T_10": "kLa_1T_10", "tau_95": "tau_95_qss"}


def r(name):
    p = ROOT / "runs" / name / "results.json"
    return json.loads(p.read_text()) if p.exists() else {}


def sets(q):
    k = KEY[q]
    a = [(r(f"pr8_l8_rpm{x[0]:.2f}_th{x[1]:.2f}").get(k), r(f"bo1_l9_rpm{x[0]:.2f}_th{x[1]:.2f}").get(k), f"{x[0]:.1f}/{x[1]:.1f}") for x in SITES]
    b = []
    for rpm in ("15", "17.5", "20", "22.5", "25", "27.5", "30", "32.5", "35", "37.5"):
        b.append((r(f"kmix_l8_rpm{rpm}").get(k), r(f"fig9_l9_rpm{rpm}").get(k), f"{rpm}/7"))
    return {"pair12": a, "theta7": b}


for q in ("kLa_1T_10", "tau_95"):
    floor = C.noise_floor(q)
    for name, P in sets(q).items():
        P = [p for p in P if p[0] and p[1] and p[0] > 0 and p[1] > 0]
        if len(P) < 5:
            print(q, name, "n", len(P), "too few"); continue
        f = np.log([p[0] for p in P]); y = np.log([p[1] for p in P]); lab = [p[2] for p in P]
        e1, e2, e0 = [], [], []
        for i in range(len(y)):
            o = np.array([j for j in range(len(y)) if j != i])
            d = np.abs(f[o] - f[i]); idx = o[np.argsort(d)]
            e1.append(y[idx[0]] - y[i]); e2.append(y[idx[:2]].mean() - y[i]); e0.append(y[o].mean() - y[i])
        rms = lambda e: float(np.sqrt(np.mean(np.square(e))))  # noqa: E731
        print(f"{q:10s} {name:7s} n={len(y):2d} NN1 {rms(e1):.3f}  NN2 {rms(e2):.3f}  mean {rms(e0):.3f}  noise floor {floor:.3f}")
        o = np.argsort(f)
        print("   sorted by L8:  " + "  ".join(f"{lab[i]} L8 {np.exp(f[i]):.3g} L9 {np.exp(y[i]):.3g}" for i in o))
