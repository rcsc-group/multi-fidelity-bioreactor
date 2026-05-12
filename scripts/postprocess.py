"""Extract kLa from BioReactor output files.

Input files (in run_dir):
  tr_oxy.dat:          i t oxy_liq_sum oxy_liq_sum2 ...  (12 cols)
  vol_frac_interf.dat: i t f_liq_sum f_liq_interf posY_max posY_min

C* (dimensionless dissolved O2 saturation) = oxy_liq_sum / f_liq_sum_mean
  where oxy_liq = f*oxy/(f*alpha+(1-f))  [liquid-phase oxygen, normalised: 1 at sat.]
  and   f_liq   = f                       [true liquid volume, not (1-cs)*f]

First-order kinetic model:  dC*/dt = kLa (1 - C*)
Exact solution:             C*(t) = 1 - exp(-kLa * t)
Linearised:                 ln(1 - C*(t)) = -kLa * t + const

kLa is estimated via a moving window of 5 consecutive time points.
kLa_10, kLa_25, kLa_50 are the kLa values at the time C* first reaches
10%, 25%, 50% saturation.  NaN is returned when data are insufficient.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


MIN_WINDOW = 5   # minimum number of time points for a valid kLa estimate


def _load_col(path: Path, col: int) -> np.ndarray:
    """Load a single column from a space-delimited dat file (skips header 'i ...')."""
    lines = [
        l for l in path.read_text().splitlines()
        if l.strip() and not l.strip().startswith("i")
    ]
    return np.array([float(l.split()[col]) for l in lines])


def _compute_c_star(run_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (t, C*) arrays from tr_oxy.dat and vol_frac_interf.dat."""
    tr_path  = run_dir / "tr_oxy.dat"
    vf_path  = run_dir / "vol_frac_interf.dat"

    t           = _load_col(tr_path, 1)     # simulation time
    oxy_liq_sum = _load_col(tr_path, 2)     # total dissolved O2 in liquid
    f_liq_sum   = _load_col(vf_path, 2)     # liquid volume (~ constant by VOF conservation)

    f_mean = f_liq_sum.mean()
    if f_mean <= 0:
        raise ValueError("f_liq_sum mean is zero — VOF field may be empty")

    c_star = oxy_liq_sum / f_mean           # dimensionless dissolved O2 (0 → 1)
    return t, c_star


def _kla_at_threshold(t: np.ndarray, c_star: np.ndarray, threshold: float) -> float:
    """Return kLa estimated by linear fit of ln(1-C*) vs t in a 5-point window
    centred on the first time C* crosses `threshold`.  Returns NaN if:
      - C* never reaches the threshold
      - fewer than MIN_WINDOW points are available around the crossing
    """
    if len(t) < MIN_WINDOW:
        return math.nan

    # find index where C* first crosses threshold
    idx = np.argmax(c_star >= threshold)
    if c_star[idx] < threshold:
        return math.nan      # threshold never reached

    # 5-point window centred on crossing (clipped to array bounds)
    half = MIN_WINDOW // 2
    lo = max(0, idx - half)
    hi = min(len(t), lo + MIN_WINDOW)
    lo = max(0, hi - MIN_WINDOW)

    t_win  = t[lo:hi]
    cs_win = c_star[lo:hi]

    # avoid log(0) or log of negative (C* may briefly exceed 1 due to VOF noise)
    cs_clip = np.clip(cs_win, 0.0, 1.0 - 1e-10)
    y = np.log(1.0 - cs_clip)

    # linear fit: y = -kLa * t + const
    slope, _ = np.polyfit(t_win, y, 1)
    return float(-slope)


def main(run_dir: str) -> dict:
    """Compute kLa at 10%, 25%, 50% saturation and write results.json.

    Parameters
    ----------
    run_dir : str
        Path to a completed BioReactor run directory containing
        tr_oxy.dat and vol_frac_interf.dat.

    Returns
    -------
    dict with keys kLa_10, kLa_25, kLa_50 (float, NaN if unavailable).
    """
    path = Path(run_dir)
    try:
        t, c_star = _compute_c_star(path)
    except (FileNotFoundError, ValueError):
        results = {"kLa_10": math.nan, "kLa_25": math.nan, "kLa_50": math.nan}
        (path / "results.json").write_text(json.dumps(results, indent=2))
        return results

    results = {
        "kLa_10": _kla_at_threshold(t, c_star, 0.10),
        "kLa_25": _kla_at_threshold(t, c_star, 0.25),
        "kLa_50": _kla_at_threshold(t, c_star, 0.50),
    }
    (path / "results.json").write_text(json.dumps(results, indent=2))
    return results


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: postprocess.py <run_dir>", file=sys.stderr)
        sys.exit(1)
    res = main(sys.argv[1])
    print(json.dumps(res, indent=2))
