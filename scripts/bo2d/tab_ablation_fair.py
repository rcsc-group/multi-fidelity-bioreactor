"""Summary of the fair ablation (mf_ablation_fair.py outputs logs/mf_abl_fair_*.csv).
Median held-out off-line L9 RMSE (ln) per qoi x n9 x on7 x l8line x variant, and the ratio to B1 (L9 only).
  uv run python scripts/bo2d/tab_ablation_fair.py"""
import glob
import pandas as pd

rows = []
for f in glob.glob("logs/mf_abl_fair_*.csv"):
    rows += [l.strip().split(",") for l in open(f) if l[:3] in ("kLa", "dtm", "tau")]
d = pd.DataFrame(rows, columns="qoi,n9,on7,l8line,rep,variant,rmse,cover90,seconds,ntest".split(","))
d[["rmse", "cover90"]] = d[["rmse", "cover90"]].astype(float)
d.to_csv("scripts/bo2d/state/mf_ablation_fair.csv", index=False)
g = d.groupby(["qoi", "n9", "on7", "l8line", "variant"]).rmse.median().unstack("variant")
print("median RMSE (ln):"); print(g.round(3).to_string())
print("\nratio to L9 only (B1):"); print(g.div(g["B1"], axis=0).round(2).to_string())
print("\nreps per cell:", d.groupby(["qoi", "n9", "on7", "l8line"]).rep.nunique().to_dict())
