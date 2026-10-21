"""Why do kLa_50 (5-sample fit) and kLa_1T_50 differ between a continuation and the reference?
Recomputes both with the CURRENT postprocess code on both runs, and reports C* at the 0.5 crossing.
Usage: uv run python scripts/diag_continue_kla50.py REF B
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402

runs = {n: ROOT / "runs" / n for n in sys.argv[1:3]}
cs = {}
for n, d in runs.items():
    t, c = pp._compute_c_star(d)
    cs[n] = (t, c)
    p = pp._load_params(d)
    _, T = pp._t_scales(p)
    i = int(np.argmax(c >= 0.5))
    print(f"{n}: crossing 0.5 at t={t[i]:.4f} (cycle {t[i] / T:.2f}); kLa_50 (5pt) "
          f"{pp._kla_5pt_at_threshold(t, c, 0.5):.6f} [nondim]; 1T "
          f"{pp._kla_period_at_threshold(t, c, 0.5, T):.6f}; C* length {len(t)}")
(ta, ca), (tb, cb) = cs.values()
n = min(len(ta), len(tb))
print("same time grid:", np.allclose(ta[:n], tb[:n]), "; max |dC*| over all rows:", np.max(abs(ca[:n] - cb[:n])))
i = int(np.argmax(ca >= 0.5))
print("C* ref / B around the ref crossing, and the 5 samples used by the 5pt fit:")
for k in range(i - 2, i + 3):
    print(f"  t={ta[k]:.3f} ref {ca[k]:.6f} B {cb[k]:.6f}")
j = int(np.argmax(cb >= 0.5))
print(f"B crossing index {j} vs ref {i} (B t={tb[j]:.3f})")
