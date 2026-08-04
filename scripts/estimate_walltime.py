"""Estimate SLURM --time for a BioReactor job from this project's own job
history (995 runs/ directories with usable logstats.dat, spanning
2026-01 to present, filtered to runs that reached >80% of their target
t_end -- see /oscar/scratch/eaguerov/tmp/walltime_formula/ for the full
data-mining scripts this table was built from, diary.md 2026-08-04).

Methodology
-----------
For each historical run, `core_sec_per_t = wall_clock_s * ntasks / t_reached`
-- a per-fidelity throughput rate (core-seconds needed per non-dimensional
simulated time unit). This is used instead of a single global regression
across (fidelity, t_end, ntasks) because those three are confounded in the
historical data: higher-fidelity runs historically used systematically
shorter t_end AND more cores, which biases a naive multi-variable fit's
individual coefficients (both came out implausibly close to zero -- see
diary.md 2026-08-04 for the fit that revealed this). Normalizing by ntasks
per fidelity level sidesteps that confounding at first order (it assumes
roughly linear core scaling, which is an approximation, not verified
independently -- MPI communication overhead at high core counts isn't
accounted for).

Per-fidelity rate table (core-seconds per non-dim time unit):
    fidelity   n     median        p90
    3          26    1.58          6.16
    4          20    14.14         22.01
    5          177   1.83          15.94
    6          11    10.92         88.84
    7          522   66.98         320.02
    8          136   1257.31       3329.28
    9          14    19868.68      27483.33
    10         22    174521.74     341465.70

fidelity=4/5 look noisy (small/mixed samples, possibly early exploratory
smoke tests not representative of production cost) -- treat estimates at
those levels with extra caution. fidelity=6/9/10 have modest sample sizes
(11-22); fidelity=7/8 are the best-supported (136-522 samples).

Cross-validated 2026-08-04 against fresh, clean, single-condition
measurements from this same session (hydrodynamics-only, n_mix_cycles=500,
no oxygen/tracer/video): L6 (ntasks=1) measured 17.4 core-sec/t vs this
table's f6 median 10.9 (1.6x); L8 (ntasks=8) measured 1828 vs f8's 1257
(1.5x); L10 (ntasks=64, post-cold-start-transient marginal rate) measured
236800 vs f10's 174522 (1.4x). All three sit modestly above the historical
median but well inside the p90 band -- the default `--stat p90` already
covers this. (Note: computing logstats.dat's #Cells sanity check requires
multiplying by ntasks first -- Basilisk's grid->n is a PER-RANK local leaf
count in MPI mode, only written by rank 0; dividing instead of multiplying
briefly looked like a fidelity-convention change across this project's
history before this was caught.)

Usage
-----
    python scripts/estimate_walltime.py --fidelity 10 --t-end 20.05 --ntasks 64
    python scripts/estimate_walltime.py --fidelity 8 --t-end 20.05 --ntasks 8 --stat median
"""
from __future__ import annotations

import argparse
import math

# {fidelity: (median_core_sec_per_t, p90_core_sec_per_t, n_samples)}
_RATE_TABLE = {
    3:  (1.578,     6.160,     26),
    4:  (14.143,    22.011,    20),
    5:  (1.831,     15.942,    177),
    6:  (10.924,    88.845,    11),
    7:  (66.976,    320.023,   522),
    8:  (1257.313,  3329.284,  136),
    9:  (19868.678, 27483.335, 14),
    10: (174521.744, 341465.697, 22),
}

# Geometric-mean per-level growth factor across the well-supported range
# (f5-f10), used only to extrapolate OUTSIDE the table (fidelity<3 or >10).
_FALLBACK_GROWTH = 8.0  # close to the theoretical NN^3 (2D, uniform-grid) factor


def estimate_wall_seconds(fidelity: int, t_end: float, ntasks: int, stat: str = "p90") -> float:
    """Return estimated wall-clock seconds for a run at (fidelity, t_end, ntasks).

    stat: "median" (typical case) or "p90" (conservative, recommended for
    actually requesting SLURM walltime -- protects against the natural
    condition-to-condition variance this table's spread already shows).
    """
    idx = 0 if stat == "median" else 1
    if fidelity in _RATE_TABLE:
        rate = _RATE_TABLE[fidelity][idx]
    else:
        # Extrapolate from the nearest known fidelity using the fallback growth factor.
        nearest = min(_RATE_TABLE, key=lambda f: abs(f - fidelity))
        rate = _RATE_TABLE[nearest][idx] * (_FALLBACK_GROWTH ** (fidelity - nearest))
    return rate * t_end / ntasks


def format_hms(seconds: float) -> str:
    seconds = max(60, int(math.ceil(seconds)))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fidelity", type=int, required=True)
    ap.add_argument("--t-end", type=float, required=True, help="non-dimensional target simulation duration")
    ap.add_argument("--ntasks", type=int, required=True)
    ap.add_argument("--stat", choices=["median", "p90"], default="p90")
    ap.add_argument("--margin", type=float, default=1.0,
                     help="extra multiplicative safety factor on top of the chosen stat (default 1.0x)")
    args = ap.parse_args()

    wall_s = estimate_wall_seconds(args.fidelity, args.t_end, args.ntasks, args.stat) * args.margin
    print(f"Estimated wall-clock ({args.stat}, x{args.margin} margin): {wall_s:.0f}s")
    print(f"Recommended --time={format_hms(wall_s)}")

    if args.fidelity not in _RATE_TABLE:
        print(f"NOTE: fidelity={args.fidelity} extrapolated outside the historical table "
              f"(fitted range: {min(_RATE_TABLE)}-{max(_RATE_TABLE)}) -- treat with extra caution.")
    elif _RATE_TABLE[args.fidelity][2] < 15:
        print(f"NOTE: fidelity={args.fidelity} has only {_RATE_TABLE[args.fidelity][2]} historical "
              "samples -- wider real uncertainty than the p90 alone suggests.")


if __name__ == "__main__":
    main()
