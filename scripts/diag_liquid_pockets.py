"""Does L9 form detached liquid pockets that L8 does not? (diagnostic for the test K L9 positive-control failure)

Hypothesis (diary 2026-10-11): at s = 1e5 the L9 chi(t) stalls at 0.85-0.90 for ~3 s, then resumes, unlike L7/L8.
A liquid pocket cut off from the main body (film, drop, trapped slug) that holds tracer would do that.
Measure on the s = 1 scal runs at the same rpm (32.5) and protocol, L8 vs L9 (liquid geometry does not depend on s):
per snapshot, liquid = (f > 0.5) & (cs > 0.5); connected components (8-connectivity); report the liquid area fraction
outside the largest component (detached fraction) and, with tracer, the share of the variance integral of c'^2 in it.
Usage: uv run python scripts/diag_liquid_pockets.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.chi_scale import locate, read_snapshot  # noqa: E402

RUNS = {8: "scal_l8_rpm32.5_th7", 9: "scal_l9_rpm32.5_th7"}
STRIDE = 5


def pocket_stats(s):
    w = s["cs"] * s["f"]
    liq = (s["f"] > 0.5) & (s["cs"] > 0.5)
    lab, n = ndimage.label(liq, structure=np.ones((3, 3)))
    if n == 0:
        return None
    sizes = ndimage.sum(liq, lab, index=np.arange(1, n + 1))
    main = 1 + int(np.argmax(sizes))
    det = liq & (lab != main)
    c = s["c"]
    cm = np.sum(w * c) / np.sum(w)
    v = w * (c - cm) ** 2
    return dict(n=n, det_area=float(det.sum() / liq.sum()), det_var=float(v[det].sum() / v.sum()),
                n_big=int(np.sum(sizes[np.arange(n) != main - 1] >= 16)))


def main():
    for L, r in RUNS.items():
        fd = locate(r)["fields_dir"]
        files = sorted(fd.glob("snap_*.bin"))[::STRIDE]
        rows = []
        for p in files:
            st = pocket_stats(read_snapshot(p, ["cs", "f", "c"]))
            if st:
                rows.append(st)
        da = np.array([x["det_area"] for x in rows])
        dv = np.array([x["det_var"] for x in rows])
        nb = np.array([x["n_big"] for x in rows])
        print(f"L{L} {r}: {len(rows)} snapshots (every {STRIDE}th of {len(files) * STRIDE})")
        print(f"   detached liquid area fraction: median {np.median(da):.4f}, 95th pct {np.percentile(da, 95):.4f}, "
              f"max {da.max():.4f}")
        print(f"   detached share of variance:    median {np.median(dv):.4f}, 95th pct {np.percentile(dv, 95):.4f}, "
              f"max {dv.max():.4f}")
        print(f"   pockets >= 16 cells: median {np.median(nb):.0f}, max {nb.max()}; "
              f"snapshots with any: {np.mean(nb > 0):.2f}")


if __name__ == "__main__":
    main()
