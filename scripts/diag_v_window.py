"""Is D_eff/D_phys (test V) flat in chi below the pre-registered window? Median r in chi bins, per scal run, plus the
time after release at which each bin starts. Decides how long a recorded L10 run must be (diagnostic).
Usage: uv run python scripts/diag_v_window.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BINS = [(0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.75), (0.75, 0.95)]


def main():
    runs = json.loads((ROOT / "experiments/multifidelity/test_v_variance_budget.json").read_text())["runs"]
    print(f"{'run':24s} " + " ".join(f"{f'{a}-{b}':>14s}" for a, b in BINS))
    for name, d in runs.items():
        if not name.startswith("scal"):
            continue
        chi, r, t = (np.array(d[k], dtype=float) for k in ("chi", "r", "t_s"))
        cells = []
        for a, b in BINS:
            m = (chi >= a) & (chi < b) & np.isfinite(r)
            cells.append(f"{np.median(r[m]):6.0f} @{t[m][0]:5.1f}s" if m.any() else f"{'n/a':>14s}")
        print(f"{name:24s} " + " ".join(f"{c:>14s}" for c in cells))


if __name__ == "__main__":
    main()
