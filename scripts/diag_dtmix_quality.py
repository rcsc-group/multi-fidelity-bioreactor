"""Quality of the mixing-time data (L6-L9 rpm sweeps, 7 deg): are the non-smooth points bad data or physics?

Question (user, 2026-10-08): the held-out L9 dtmix points of the H1 test look non-smooth in rpm; are they bad
data points rather than poor predictive power of the MF model?
Checks, per threshold chi in {0.50, 0.75, 0.95}:
  1. Roughness in rpm: r_i = log y(rpm_i) - mean(log y(rpm_{i-1}), log y(rpm_{i+1})) (interior points).
     Compared with the run-to-run noise in log units (replicate sd, releases 80/82/85, L6-L8 and L9 at 25 rpm).
     A point whose |r| is many noise sd is not explained by run noise.
  2. Does the same rpm deviate at every level and in Kim et al.? A dip that persists across levels and in Kim is
     a property of the flow; one present at one level only is suspect.
Writes experiments/multifidelity/dtmix_quality.json and prints tables.
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
CHIS = ("0.50", "0.75", "0.95")
RUN = {6: "kmix_l6_rpm{:g}", 7: "kmix_l7_rpm{:g}", 8: "kmix_l8_rpm{:g}", 9: "fig9_l9_rpm{:g}"}
KIM = pd.read_csv(ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_frequency.csv").set_index("RPM")


def val(run, key):
    f = ROOT / "runs" / run / "results.json"
    return json.load(open(f)).get(key, math.nan) if f.exists() else math.nan


def replicate_sd(level, rpm, key):
    """sd of log y over the base run and its replicates (rep_l{L}_rpm{r}_rel{80,82,85})."""
    vals = [val(RUN[level].format(rpm), key)]
    vals += [val(f"rep_l{level}_rpm{rpm:g}_rel{rel}", key) for rel in (80, 82, 85)]
    v = np.log([x for x in vals if x == x and x > 0])
    return (float(np.std(v, ddof=1)), len(v)) if len(v) >= 2 else (math.nan, len(v))


def rough(logy):
    r = np.full(len(logy), np.nan)
    r[1:-1] = logy[1:-1] - 0.5 * (logy[:-2] + logy[2:])
    return r


def main():
    out = {}
    for chi in CHIS:
        key = f"dtmix_{chi}"
        print(f"\n=== dtmix chi {chi} (s); roughness r = log y - mean of neighbours")
        tab = {lev: np.array([val(RUN[lev].format(r), key) for r in RPMS]) for lev in RUN}
        tab["Kim"] = np.array([KIM.loc[r, f"dtmix_strict_{float(chi):g}"] if r in KIM.index else math.nan
                               for r in RPMS])
        R = {k: rough(np.log(v)) for k, v in tab.items()}
        noise = {}
        for lev in (6, 7, 8, 9):
            for rpm in (17.5, 25, 32.5):
                sd, n = replicate_sd(lev, rpm, key)
                if n >= 2:
                    noise[f"L{lev}@{rpm:g}"] = (sd, n)
        hdr = "rpm   " + "".join(f"{('L' + str(k)) if k != 'Kim' else 'Kim':>16s}" for k in tab)
        print(hdr)
        for i, rpm in enumerate(RPMS):
            print(f"{rpm:5g} " + "".join(f"{tab[k][i]:8.1f} ({R[k][i]:+.2f})" for k in tab))
        print("replicate sd of log y [n]: " + ", ".join(f"{k} {s:.3f} [{n}]" for k, (s, n) in noise.items()))
        # robust roughness scale per series: median |r| over interior points
        print("median |r| per series: " + ", ".join(f"{k}: {np.nanmedian(np.abs(R[k])):.2f}" for k in tab))
        out[chi] = dict(rpm=RPMS, values={str(k): v.tolist() for k, v in tab.items()},
                        roughness={str(k): v.tolist() for k, v in R.items()},
                        replicate_sd={k: s for k, (s, n) in noise.items()})
    json.dump(out, open(ROOT / "experiments/multifidelity/dtmix_quality.json", "w"), indent=1)


if __name__ == "__main__":
    main()
