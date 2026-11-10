"""Held-out RMSE includes the L9 run-to-run noise, which no model can remove. Reducible error = sqrt(RMSE^2 - floor^2)
(floor = C.noise_floor, the measured repeat-run sd in ln). Median over reps, per qoi x n9, for the uniform candidates."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "campaign2d"))
import campaign as C  # noqa: E402
d = pd.read_csv("scripts/bo2d/state/mf_ablation.csv")
fl = {q: C.noise_floor(q) for q in d.qoi.unique()}
print("noise floors (ln):", {k: round(v, 3) for k, v in fl.items()})
d["excess"] = np.sqrt(np.clip(d.rmse ** 2 - d.qoi.map(fl) ** 2, 0, None))
t = d[d.variant.isin(["B1", "U1", "B2-lin", "B3-lin", "U2-lin", "B2-quad"])].groupby(["qoi", "variant", "n9"]).excess.median().unstack("n9").round(3)
print(t.to_string())
r = t.div(t.xs("B1", level="variant"), level="qoi").round(2)
print("\nreducible error relative to L9-only (B1):\n" + r.to_string())
