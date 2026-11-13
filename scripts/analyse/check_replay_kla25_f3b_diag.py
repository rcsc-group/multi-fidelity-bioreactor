"""Task 3 diagnosis (2026-11-04), read-only: why does the 22-cycle K=1 replay (f3b_loop1) give +51 % on kLa_1T_25
while the 150-cycle chain replay rc8_rpm32.5_rep gives +1 %? Compares C*(cycle) and the t_cross of C*=0.25 for
f3b_coupled, f3b_loop1 and the database kmix_l8_rpm32.5 / rc8_rpm32.5_rep. Prints the replay params that differ.
Usage: uv run python scripts/analyse/check_replay_kla25_f3b_diag.py
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402

runs = {"f3b_coupled": rg.run_dir("f3b_coupled"), "f3b_loop1": rg.run_dir("f3b_loop1"),
        "db_kmix32.5": ROOT / "runs" / "kmix_l8_rpm32.5", "rc8_32.5_rep": rg.RUNS / "rc8_rpm32.5_rep"}
probe = [5, 10, 15, 17, 20, 22]
for name, d in runs.items():
    params = json.loads((d / "params.json").read_text())
    t, c = pp._compute_c_star(d)
    _, T = pp._t_scales(params)
    cyc = t / T
    k = int(np.argmax(c >= 0.25))
    row = [f"{p}:{np.interp(p, cyc, c):.4f}" for p in probe]
    print(f"{name:14s} t_cross25 {cyc[k]:6.2f} cyc  C* at cycle " + " ".join(row) + f"  first t {cyc[0]:.2f}")
    keys = ("_binary", "n_rec_periods", "t_checkpoint", "replay_loop", "restart_continue", "n_mix_cycles", "t_end")
    print("   ", {k2: params.get(k2) for k2 in keys})
