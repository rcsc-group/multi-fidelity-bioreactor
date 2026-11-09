"""chi(t) after release in every BO L9 run: cycles simulated after release, max chi reached, and time to
chi = 0.10 / 0.25 / 0.50 in cycles. Question (user 2026-11-09): why is dtmix NaN in most BO L9 runs?
Usage: uv run python scripts/bo2d/diag_chi_bo.py [glob]"""
import json
import sys
from pathlib import Path

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.postprocess import _t_scales  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
pat = sys.argv[1] if len(sys.argv) > 1 else "bo1_l9_*"
print(f"{'run':28s} {'T':>6s} {'cyc_post':>8s} {'chi_max':>7s} {'c@.10':>6s} {'c@.25':>6s} {'c@.50':>6s}  sig2max")
for p in sorted((ROOT / "runs").glob(pat)):
    tr, rel = p / "tr_oxy.dat", p / "mix_release.dat"
    if not (tr.exists() and rel.exists()):
        continue
    pp = json.loads((p / "params.json").read_text())
    T = _t_scales(pp)[1]
    a = np.loadtxt(tr, skiprows=1)
    r = np.atleast_2d(np.loadtxt(rel, skiprows=1))[0]
    s2m = r[4] / r[2] - (r[3] / r[2]) ** 2
    m = a[:, 1] > r[1]
    t, b0, b1, b2 = a[m, 1], a[m, 12], a[m, 13], a[m, 14]
    chi = 1 - (b2 / b0 - (b1 / b0) ** 2) / s2m
    out = []
    for th in (0.10, 0.25, 0.50):
        i = np.argmax(chi >= th)
        out.append(f"{(t[i] - r[1]) / T:6.1f}" if chi[i] >= th else "   nan")
    print(f"{p.name:28s} {T:6.3f} {(t[-1] - r[1]) / T:8.1f} {chi.max():7.3f} {' '.join(out)}  {s2m:.4f}")
