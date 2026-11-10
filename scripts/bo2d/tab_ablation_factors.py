"""Factor effects of the model ablation (mf_ablation.csv), paired by (qoi, n9, rep) and by the other two factors.
For each factor (fully Bayesian vs ML; 3 vs 2 levels; quadratic vs linear): median RMSE ratio (B/A) and fraction of
paired cases where B is better; per qoi and n9. Also one uniform method score: mean over qoi and n9 of RMSE / L9-only (B1)."""
import numpy as np
import pandas as pd
d = pd.read_csv("scripts/bo2d/state/mf_ablation.csv")
w = d.pivot_table(index=["qoi", "n9", "rep"], columns="variant", values="rmse")
c = d.pivot_table(index=["qoi", "n9", "rep"], columns="variant", values="cover90")
pairs = {"fully Bayesian / ML": [(f"U{l}-{b}", f"B{l}-{b}") for l in (2, 3) for b in ("lin", "quad")] + [("U1", "B1")],
         "3 levels / 2 levels": [(f"{h}2-{b}", f"{h}3-{b}") for h in "UB" for b in ("lin", "quad")],
         "quadratic / linear": [(f"{h}{l}-lin", f"{h}{l}-quad") for h in "UB" for l in (2, 3)]}
for name, pp in pairs.items():
    r = pd.concat([np.log(w[b] / w[a]) for a, b in pp], axis=0)
    g = r.groupby(level=["qoi", "n9"])
    out = pd.DataFrame({"median ratio": np.exp(g.median()), "B better": g.apply(lambda s: (s < 0).mean())}).round(2)
    print(f"\n== {name} (ratio < 1 = second option better)\n" + out.unstack("n9").to_string())
cv = c.groupby(level="qoi").median().T.round(2)
print("\n== median 90% coverage per variant (target 0.90)\n" + cv.to_string())
s = (w.div(w["B1"], axis=0)).groupby(level=["qoi", "n9"]).median()
print("\n== uniform-method score: median RMSE / L9-only(B1), then mean over qoi x n9 (lower = better)\n" +
      s.mean().sort_values().round(3).to_string())
print("\n== same, per qoi (mean over n9)\n" + s.groupby(level="qoi").mean().T.round(2).to_string())
