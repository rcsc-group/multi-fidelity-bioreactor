"""Truth snapshot v2 (2026-11-11 08:40, user go for a zero-compute rerun): v1 snapshot minus the 8 c2d L9 runs
(released after 12 cycles on the L9 grid; the wave needs ~27; measured kLa bias +2..+21 % vs proper 30-cycle runs at
4 near sites, diag_c2d_vs_proper.py), plus every proper L9 run finished since the freeze (bo1_/bo4_, same gate as v1
via loop.rows). L7 / L8 rows unchanged. Writes state/truth_snapshot_v2_<qoi>.csv (read with env BO_TRUTH=_v2).
  uv run python scripts/bo2d/freeze_truth_v2.py
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import loop as LP  # noqa: E402

S = Path(__file__).parent / "state"
for q in ("kLa_1T_10", "dtmix_0.25", "tau_95"):
    d = pd.read_csv(S / f"truth_snapshot_{q}.csv")
    drop = (d.level == 9) & d.run.str.startswith("c2d_")
    new = LP.rows(q, ["bo1_", "bo4_"], {9})[["level", "rpm", "theta", "run", "y"]]
    new = new[~new.run.isin(d.run)]
    out = pd.concat([d[~drop], new], ignore_index=True)
    out.to_csv(S / f"truth_snapshot_v2_{q}.csv", index=False)
    print(f"{q}: L9 {int((d.level == 9).sum())} -> {int((out.level == 9).sum())} (dropped c2d {int(drop.sum())}, added {len(new)})")
