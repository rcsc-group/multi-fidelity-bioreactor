"""Can an early-saturation run predict a late-saturation quantity? (pilot, no new runs)

Proposal: low fidelity = run until C* = 10%, high fidelity = run until C* = 95%,
and learn the 95% quantity from the cheap runs. No run reaches 95%, so the
same question is asked one step down, where both ends exist:
  oxygen : kLa at C* = 10%  ->  kLa at 25% / 50%     (kLa_1T_*, whole-period fit)
  mixing : dtmix at chi = 0.50  ->  dtmix at chi = 0.95
over the L9 rpm sweep (7 deg) and the L9 angle sweep (32.5 rpm).

What makes the LF informative (Yi et al.'s transfer): a strong, smooth
relation across conditions. Reported: Pearson r on logs, Spearman rank r, and
the spread of the ratio HF/LF. If the ratio were constant, LF alone would do.
If r is weak, the transfer has nothing to learn from.
Also: how far each run actually got (max C*), and t(C*) in rocking cycles.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _compute_c_star, _t_scales  # noqa: E402

SWEEPS = {"rpm, 7 deg": [f"fig9_l9_rpm{r:g}" for r in (15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5)],
          "angle, 32.5 rpm": [f"fig10_l9_th{a:g}" for a in (2, 3, 4, 5, 6, 7)]}
PAIRS = [("kLa_1T_10", "kLa_1T_25"), ("kLa_1T_10", "kLa_1T_50"),
         ("dtmix_0.50", "dtmix_0.95"), ("dtmix_0.75", "dtmix_0.95"),
         # does bulk MIXING predict the late transfer rate better than the
         # early rate does? (inverse mixing time vs kLa; only the sign of the
         # expected relation is assumed)
         ("inv_dtmix_0.95", "kLa_1T_50"), ("inv_dtmix_0.50", "kLa_1T_50"),
         ("inv_dtmix_0.95", "kLa_1T_25")]


def cycles_to(run, level):
    d = ROOT / "runs" / run
    p = json.load(open(d / "params.json"))
    _, T = _t_scales(p)
    t, c = _compute_c_star(d)
    i0 = int(np.argmax(c > 1e-6))
    k = int(np.argmax(c >= level))
    return (t[k] - t[i0]) / T if c[k] >= level else math.nan, float(c.max())


def main() -> None:
    for name, runs in SWEEPS.items():
        res = {r: json.load(open(ROOT / "runs" / r / "results.json")) for r in runs}
        for d in res.values():
            for k in ("dtmix_0.50", "dtmix_0.95"):
                v = d.get(k, math.nan)
                d["inv_" + k] = 1.0 / v if v and v == v else math.nan
        print(f"\n== {name}")
        print(f"  {'run':18s} {'maxC*':>6} {'cyc->10%':>9} {'cyc->25%':>9} {'cyc->50%':>9}")
        for r in runs:
            c10, cmax = cycles_to(r, 0.10)
            c25, _ = cycles_to(r, 0.25)
            c50, _ = cycles_to(r, 0.50)
            print(f"  {r:18s} {cmax:6.2f} {c10:9.1f} {c25:9.1f} {c50:9.1f}")
        for lf, hf in PAIRS:
            x = np.array([res[r].get(lf, np.nan) for r in runs], float)
            y = np.array([res[r].get(hf, np.nan) for r in runs], float)
            m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
            if m.sum() < 4:
                print(f"  {lf:>10} -> {hf:<10}: only {m.sum()} pairs")
                continue
            rp = pearsonr(np.log(x[m]), np.log(y[m]))[0]
            rs = spearmanr(x[m], y[m])[0]
            ratio = y[m] / x[m]
            print(f"  {lf:>10} -> {hf:<10}: n={m.sum():2d}  r_log={rp:+.2f}  rank={rs:+.2f}  "
                  f"HF/LF ratio {ratio.min():.2f}-{ratio.max():.2f} (cv {ratio.std() / ratio.mean():.0%})")


if __name__ == "__main__":
    main()
