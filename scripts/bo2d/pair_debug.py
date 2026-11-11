"""Debug pair_score: which bo1 L9 rows load per QoI at the 12 pair sites, and which pr8 L8 rows."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import loop as LP  # noqa: E402
from scripts.bo2d.pair_l8 import SITES  # noqa: E402

for q in ("kLa_1T_10",):
    d9 = LP.rows(q, ["bo1_"], {9}); d8 = LP.rows(q, ["pr8_"], {8})
    print(q, "bo1 L9 rows", len(d9), "pr8 L8 rows", len(d8))
    print(d9[["rpm", "theta", "run", "y"]].to_string())
    print(d8[["rpm", "theta", "run", "y"]].to_string())
for s in SITES:
    d = LP.C.ROOT / "runs" / f"bo1_l9_rpm{s[0]:.2f}_th{s[1]:.2f}"
    print(s, d.name, d.exists(), (d / "results.json").exists() if d.exists() else "")
