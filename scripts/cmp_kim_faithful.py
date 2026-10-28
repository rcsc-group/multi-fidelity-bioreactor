"""Compare the Kim-style replay (kim_l8_rpm<r>_th<theta>) with the faithful replay (lf_l8_rpm<r>_th<theta>), same record protocol (scripts/replay_kim.py).
Prints the QoIs and, from tr_oxy.dat, the tracer total c2_liq_sum and the oxygen total oxy_liq_sum at chosen replay cycles after the first row.
Usage: uv run python scripts/cmp_kim_faithful.py [rpm theta]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402

rpm, theta = (float(sys.argv[1]), float(sys.argv[2])) if len(sys.argv) == 3 else (32.5, 7.0)
QOIS = ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "dtmix_0.50", "dtmix_0.75")
runs = {"faithful": f"lf_l8_rpm{rpm:g}_th{theta:g}", "kim": f"kim_l8_rpm{rpm:g}_th{theta:g}"}
data = {}
for k, rid in runs.items():
    d = rg.RUNS / rid
    res = json.loads((d / "results.json").read_text())
    _, T = pp._t_scales(json.loads((d / "params.json").read_text()))
    a = np.loadtxt(d / "tr_oxy.dat", skiprows=1)
    data[k] = (res, T, a)
print("QoI            faithful        kim   kim/faithful")
for q in QOIS:
    f, m = data["faithful"][0].get(q), data["kim"][0].get(q)
    print(f"{q:12s} {f:10.4f} {m:10.4f} {m / f:10.3f}" if f and m else f"{q:12s} {f} {m}")
cols = {"oxygen": 2, "tracer": 8}   # oxy_liq_sum, c2_liq_sum
print("\nratio of totals to the first row (t = release):  cycle  " + "  ".join(f"{k}:{n}" for k in runs for n in cols))
cyc = (0.5, 1, 2, 5, 10, 15, 20)
for c in cyc:
    row = []
    for k in runs:
        _, T, a = data[k]
        t0 = a[0, 1]
        i = int(np.argmin(np.abs(a[:, 1] - (t0 + c * T))))
        row += [a[i, cols[n]] / a[0, cols[n]] for n in cols]
    print(f"                                            {c:5g}  " + "  ".join(f"{v:9.4f}" for v in row))
for k in runs:
    _, T, a = data[k]
    print(f"{k}: rows {len(a)}, replay span {(a[-1, 1] - a[0, 1]) / T:.2f} cycles")
