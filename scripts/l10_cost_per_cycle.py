"""Measured L10 cost per rocking cycle from logstats.dat (wall clock vs simulation time), 32 ranks. Usage: uv run python scripts/l10_cost_per_cycle.py"""
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path[:0] = ["."]
from scripts.postprocess import _t_scales  # noqa: E402

RANKS = 32
for run in sorted(p.name for p in Path("runs").glob("t10_l10_rpm*_th*")) + ["v_l10_rpm32.5_th7"]:
    if run.endswith("_ext"):
        continue
    d = Path("runs") / run
    try:
        p = json.load(open(d / "params.json")); _, T = _t_scales(p)
        rows = [(float(m.group(1)), float(m.group(2))) for m in (re.search(r"t: (\S+) .*Wall clock time \(s\): (\S+)", l) for l in open(d / "logstats.dat")) if m]
        t, w = np.array(rows).T
        c = t / T
        sel = c > c.max() - 5   # last 5 cycles
        cyc_cost = (w[sel][-1] - w[sel][0]) / (c[sel][-1] - c[sel][0]) * RANKS / 3600
        print(f"{run:28s} cycles done {c.max():6.1f}  core-h per cycle (last 5 cycles) {cyc_cost:6.1f}")
    except Exception as e:  # noqa: BLE001
        print(run, "skip", repr(e)[:60])
