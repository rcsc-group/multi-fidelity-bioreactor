"""How far over the cap are the 'infeasible' recommendations of diag_mf_recommend? (2026-11-10)
Hypothesis: the chain recommends ON the constraint boundary (its tau band is narrow), so tiny errors flip feasibility;
the excess is then below the tau_95 noise floor. Prints the truth excess ln(tau9/cap) of every infeasible row, per arm,
and the noise floor. Usage: uv run python scripts/bo2d/diag_infeas_margin.py"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

os.environ["BO_OBJ"] = "mix"
sys.path.insert(0, str(Path(__file__).parent))
import synth_bench as S  # noqa: E402

S.setup("real")
d = pd.read_csv(Path(__file__).parent / "state" / "diag_mfrec.csv").dropna(subset=["rpm"])
d["excess"] = S.F[9][1](d[["rpm", "theta"]].to_numpy(float)) - np.log(0.014)
print(f"tau_95 noise floor (ln): {S.C.noise_floor('tau_95'):.3f}")
bad = d[~d.feasible.astype(bool)]
print(bad.groupby(["obj", "arm"]).excess.describe().round(3).to_string())
for tol in (0.02, 0.05, 0.10):
    ok = d.excess <= tol
    print(f"feasible within ln-tolerance {tol}: " + ", ".join(f"{k}: {v:.2f}" for k, v in ok.groupby([d.obj, d.arm]).mean().items()))
