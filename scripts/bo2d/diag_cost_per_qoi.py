"""Cycles after release each QoI needs, per BO L9 run (user 2026-11-09: is kLa or t_mix cheaper?).
Oxygen: cycles until C* = 0.10 / 0.25 (+1 period for the 1T fit). Tracer: cycles until chi = 0.25 / 0.50.
Usage: uv run python scripts/bo2d/diag_cost_per_qoi.py [glob]"""
import json
import sys
from pathlib import Path

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.postprocess import _compute_c_star, _t_scales  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
pat = sys.argv[1] if len(sys.argv) > 1 else "bo1_l9_*"
print(f"{'run':26s} {'C*.10':>6s} {'C*.25':>6s} {'chi.25':>6s} {'chi.50':>6s}   (cycles after release; +1 for 1T fit on C*)")
for p in sorted((ROOT / "runs").glob(pat)):
    rel = p / "mix_release.dat"
    if not rel.exists():
        continue
    pp = json.loads((p / "params.json").read_text()); T = _t_scales(pp)[1]
    r = np.atleast_2d(np.loadtxt(rel, skiprows=1))[0]
    t, c = _compute_c_star(p)
    a = np.loadtxt(p / "tr_oxy.dat", skiprows=1); m = a[:, 1] > r[1]
    chi = 1 - (a[m, 14] / a[m, 12] - (a[m, 13] / a[m, 12]) ** 2) / (r[4] / r[2] - (r[3] / r[2]) ** 2)
    def first(tt, v, th):
        i = np.argmax(v >= th); return f"{(tt[i] - r[1]) / T:6.1f}" if v[i] >= th else "   nan"
    print(f"{p.name:26s} {first(t, c, .10)} {first(t, c, .25)} {first(a[m, 1], chi, .25)} {first(a[m, 1], chi, .50)}")
