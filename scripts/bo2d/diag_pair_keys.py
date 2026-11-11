"""Root cause hunt (2026-11-10 22:40): at the 12 pair sites (pr8_l8 vs bo1_l9, same x), which stored quantities
track across levels? For every numeric key in results.json common to both: Spearman rho(L8, L9), median L8/L9 ratio.
Gate ignored (all 12 used). Also kLa-type keys from tr_oxy.dat are already in results.json.
  uv run python scripts/bo2d/diag_pair_keys.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.bo2d.pair_l8 import SITES  # noqa: E402


def rid(l, x, tag):
    return f"{tag}_l{l}_rpm{x[0]:.2f}_th{x[1]:.2f}"


def num(d):
    return {k: float(v) for k, v in d.items() if isinstance(v, (int, float)) and not isinstance(v, bool)}


R8, R9 = [], []
for x in SITES:
    a = ROOT / "runs" / rid(8, x, "pr8") / "results.json"
    b = ROOT / "runs" / rid(9, x, "bo1") / "results.json"
    R8.append(num(json.loads(a.read_text()))); R9.append(num(json.loads(b.read_text())))
keys = sorted(set.intersection(*(set(r) for r in R8 + R9)))
rows = []
for k in keys:
    a = np.array([r[k] for r in R8]); b = np.array([r[k] for r in R9])
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 8 or np.std(b[ok]) == 0:
        continue
    rho = spearmanr(a[ok], b[ok]).correlation
    ratio = np.median(a[ok] / b[ok]) if np.all(b[ok] != 0) else np.nan
    rows.append((rho, k, ok.sum(), ratio, np.std(np.log(np.abs(b[ok]) + 1e-30))))
for rho, k, n, ratio, sd in sorted(rows, key=lambda r: -abs(r[0]) if np.isfinite(r[0]) else 0):
    print(f"{k:32s} n={n:2d} rho {rho:+.2f}  L8/L9 {ratio:8.3g}  sd ln L9 {sd:.2f}")
th = np.array([x[1] for x in SITES]); rpm = np.array([x[0] for x in SITES])
for k in ("kLa_1T_10",):
    if k in keys:
        a = np.array([r[k] for r in R8]); b = np.array([r[k] for r in R9])
        print("\nsite  rpm theta  L8  L9  L8/L9")
        for i in np.argsort(rpm):
            print(f"{rpm[i]:6.2f} {th[i]:4.2f} {a[i]:6.1f} {b[i]:6.1f} {a[i] / b[i]:5.2f}")
