"""Dump every L9 row we have (archive + BO arms) for kLa_1T_10, tau_95, dtmix_0.50, with the feasibility flag.
Usage: uv run python scripts/bo2d/dump_l9.py"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import loop as L  # noqa: E402

C = L.C
cap = json.loads(L.CFG.read_text())["cap"] if L.CFG.exists() else None
print("cap (config):", cap)
out = {}
for q in ("kLa_1T_10", "tau_95", "dtmix_0.50"):
    try:
        db = C.load_db(q)
        db = db[db.level == 9][["rpm", "theta", "run", "y"]]
    except Exception as e:  # noqa: BLE001
        print("load_db", q, "failed:", e); db = pd.DataFrame(columns=["rpm", "theta", "run", "y"])
    bo = L.rows(q, ["bo1_", "bo4_"], {9})[["rpm", "theta", "run", "y"]]
    out[q] = pd.concat([db, bo]).drop_duplicates("run").set_index("run")
d = out["kLa_1T_10"].rename(columns={"y": "kLa"})
d["tau95"] = out["tau_95"].y
d["dtmix50"] = out["dtmix_0.50"].y
d = d.sort_values(["theta", "rpm"])
if cap is not None:
    d["feas"] = d.tau95 <= cap
pd.set_option("display.width", 200)
print(d.to_string())
print("n =", len(d), " kLa range", d.kLa.min(), d.kLa.max())
for th in sorted(d.theta.round(1).unique()):
    s = d[d.theta.round(1) == th]
    if len(s) >= 3:
        print(f"theta {th}: Spearman(rpm, kLa) = {s.rpm.corr(s.kLa, method='spearman'):+.2f}  n={len(s)}")
print("Spearman(theta, kLa) all:", round(d.theta.corr(d.kLa, method="spearman"), 2),
      " (rpm, kLa):", round(d.rpm.corr(d.kLa, method="spearman"), 2))
