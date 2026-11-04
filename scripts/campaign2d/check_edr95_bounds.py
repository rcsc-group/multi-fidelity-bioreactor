"""Per-sample bound check for edr_95 (2026-11-04).

For every output sample i of a run directory:
    edr_mean(i) = shear_stress.dat column 9 (domain-mean EDR, nondim)
    edr_95(i)   = edr_pct.dat column 2      (95th percentile of the local EDR, nondim)
    edr_max(i)  = edr_pct.dat column 3      (spatial maximum of the local EDR, nondim)
Checks edr_mean <= edr_95 <= edr_max at every sample. Exit code 1 if any sample violates it.

Usage: uv run python scripts/campaign2d/check_edr95_bounds.py <run_dir>
"""
import sys
from pathlib import Path

import numpy as np


def load(path: Path) -> np.ndarray:
    rows = [l.split() for l in path.read_text().splitlines() if l.strip() and not l.strip().startswith("i")]
    return np.array([[float(x) for x in r] for r in rows])


def main() -> int:
    run = Path(sys.argv[1])
    ss = load(run / "shear_stress.dat")
    ed = load(run / "edr_pct.dat")
    n = min(len(ss), len(ed))
    if not np.array_equal(ss[:n, 0], ed[:n, 0]):
        print("FAIL: sample index mismatch between shear_stress.dat and edr_pct.dat")
        return 1
    mean, p95, mx = ss[:n, 9], ed[:n, 2], ed[:n, 3]
    tol = 1e-12
    v_lo = np.where(mean > p95 + tol)[0]
    v_hi = np.where(p95 > mx + tol)[0]
    print(f"samples: {n}")
    print(f"edr_mean > edr_95 at {len(v_lo)} samples" + (f", max excess {np.max(mean[v_lo] - p95[v_lo]):.3e}" if len(v_lo) else ""))
    print(f"edr_95 > edr_max at {len(v_hi)} samples")
    print(f"ratio edr_95/edr_mean over samples with mean > 0: min {np.min(p95[mean > 0] / mean[mean > 0]):.3f}, "
          f"max {np.max(p95[mean > 0] / mean[mean > 0]):.3f}")
    print(f"ratio edr_max/edr_95 over samples with edr_95 > 0: min {np.min(mx[p95 > 0] / p95[p95 > 0]):.3f}, "
          f"max {np.max(mx[p95 > 0] / p95[p95 > 0]):.3f}")
    return 1 if (len(v_lo) or len(v_hi)) else 0


if __name__ == "__main__":
    sys.exit(main())
