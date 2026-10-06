"""Why is the h -> 0 epistemic band wide? Compare it with the data's own level-to-level differences.

For each qoi and rpm (L8-L10 fit, experiments/multifidelity/gcbml_hydro_v3_twy2_L8_10.json):
  d98  = |f9 - f8| / f10, d109 = |f10 - f9| / f10   (relative level differences)
  rich = d109 / (2^p - 1)                          Richardson estimate of the remaining error at L10, p = fitted median
  half = (epi_q975 - epi_q025) / 2 / m              half-width of the 95% epistemic band, relative
Run: uv run --project /oscar/data/dharri15/eaguerov/Github/gcbml python scripts/hydro_band_diag.py
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
rows = json.load(open(ROOT / "experiments/multifidelity/hydro_dataset_v2.json"))
fits = json.load(open(ROOT / "experiments/multifidelity/gcbml_hydro_v3_twy2_L8_10.json"))
for fit in fits:
    key = fit["qoi"]
    p = fit["p0_q"][1]
    rpm = np.array(fit["rpm"])
    half = (np.array(fit["epi_q975"]) - np.array(fit["epi_q025"])) / 2 / np.array(fit["m"])
    print(f"\n{key}: median p = {p:.2f}")
    print(" rpm   f8       f9       f10     |f9-f8|/f10 |f10-f9|/f10  Richardson  band half-width")
    for r in sorted({r["rpm"] for r in rows if r["level"] == 10}):
        f = {lev: next((x[key] for x in rows if x["rpm"] == r and x["level"] == lev), np.nan) for lev in (8, 9, 10)}
        d98, d109 = abs(f[9] - f[8]) / f[10], abs(f[10] - f[9]) / f[10]
        h = half[np.argmin(abs(rpm - r))]
        print(f"{r:5.1f} {f[8]:.5f} {f[9]:.5f} {f[10]:.5f}   {d98:6.1%}      {d109:6.1%}     {d109 / (2**p - 1):6.1%}      {h:6.1%}")
