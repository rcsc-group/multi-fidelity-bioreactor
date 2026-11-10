"""Table of diag_mf_recommend results (scripts/bo2d/state/diag_mfrec.csv). Infeasible recommendation = failure."""
import numpy as np
import pandas as pd
d = pd.read_csv("scripts/bo2d/state/diag_mfrec.csv")
d["r"] = np.where(d.feasible, d.regret.clip(lower=0), np.nan)
g = d.groupby(["obj", "arm", "n9"]).agg(core=("core_h", "first"), med=("r", "median"), feas=("feasible", "mean"),
                                         hit10=("r", lambda s: (s <= 0.10).mean()), hit20=("r", lambda s: (s <= 0.20).mean()), n=("r", "size"))
print(g.round(2).to_string())
