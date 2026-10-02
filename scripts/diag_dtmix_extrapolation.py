"""Per-rpm Eca & Hoekstra extrapolation of dtmix (chi 0.50/0.75/0.95) from L6-L9.

Target is the PHYSICAL mixing time (h -> 0), not L9. Levels: kmix_l6 (cold),
kmix_l7 (warm, verified within 2% of the cold L7 at 32.5 rpm), kmix_l8 (cold),
fig9_l9 (cold). Prints phi0 +- U(phi0-scale) beside Kim, plus the "cheap" variant
L6-L8 only (3 levels: no E&H, so a p=2 Richardson triplet) to see whether the
expensive L9 is needed at all.
"""
import json, math, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.eca_hoekstra import numerical_uncertainty  # noqa: E402

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
RUNS = {9: "fig9_l9_rpm{:g}", 8: "kmix_l8_rpm{:g}", 7: "kmix_l7_rpm{:g}", 6: "kmix_l6_rpm{:g}"}
KIM = pd.read_csv(ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_frequency.csv").set_index("RPM")


def val(run, key):
    f = ROOT / "runs" / run / "results.json"
    return json.load(open(f)).get(key, math.nan) if f.exists() else math.nan


def table():
    rows = []
    for thr in (0.50, 0.75, 0.95):
        for r in RPMS:
            q = np.array([val(RUNS[N].format(r), f"dtmix_{thr:.2f}") for N in (9, 8, 7, 6)])
            kim = float(KIM.loc[r, f"dtmix_strict_{thr:g}"])
            row = dict(thr=thr, rpm=r, L9=q[0], L8=q[1], L7=q[2], L6=q[3], kim=kim,
                       phi0=math.nan, U=math.nan, p=None, est="")
            if np.all(np.isfinite(q)):
                h = np.array([2.0, 4.0, 8.0, 16.0])
                res = numerical_uncertainty(h, q, index=0)
                # U is the 95% bound for L9 around phi_exact; for phi0 itself use
                # the same U (its distance to L9 is already inside it)
                row.update(phi0=res["phi0"], U=res["U"], p=res["p"], est=res["estimator"])
            rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    t = table()
    pd.set_option("display.width", 200)
    t["phi0/Kim"] = t.phi0 / t.kim
    t["L9/Kim"] = t.L9 / t.kim
    t["U/phi0"] = t.U / t.phi0
    print(t.round(3).to_string(index=False))
