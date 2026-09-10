"""Median relative error of tau_mean_max (Ours <tau>) at L9 vs Kim's
published <tau> curve, across the 9-point rpm sweep.

Usage:
    uv run python scripts/compute_l9_median_relerr.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent.parent
KIM_CSV = HERE / "experiments/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv"
RUNS_DIR = HERE / "runs"

kim = pd.read_csv(KIM_CSV, skiprows=[1])
kim["RPM"] = pd.to_numeric(kim["RPM"])
kim = kim.sort_values("RPM").set_index("RPM")

L9_RUN_IDS = {
    17.5: "l9_sweep_rpm17.5", 20.0: "20a14369", 22.5: "255e0b87",
    25.0: "bce29aa5", 27.5: "9ec56180", 30.0: "60d09a80",
    32.5: "a34fc4d4", 35.0: "a281a16f", 37.5: "0f0ad3ea",
}

rows = []
for rpm, rid in sorted(L9_RUN_IDS.items()):
    ours = json.loads((RUNS_DIR / rid / "results.json").read_text())["tau_mean_max"]
    ref = kim.loc[rpm, "tau_liq_mean"]
    relerr = abs(ours - ref) / ref
    rows.append((rpm, ours, ref, relerr))
    print(f"  {rpm:5.1f} rpm: ours={ours:.5f} Pa  kim={ref:.5f} Pa  relerr={100*relerr:.1f}%")

relerrs = np.array([r[3] for r in rows])
print(f"\nmedian relative error (tau_mean_max, L9 vs Kim <tau>): {100*np.median(relerrs):.1f}%")
print(f"mean relative error: {100*np.mean(relerrs):.1f}%")
