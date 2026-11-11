"""How biased are the c2d L9 runs (release after 12 cycles on the L9 grid) vs proper L9 runs (30 cycles, bo1/bo4)?
(2026-11-11 08:30). For each c2d_l9 site: nearest proper L9 run (scaled distance), and kLa_1T_10, tau_95_qss of both.
  uv run python scripts/bo2d/diag_c2d_vs_proper.py
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2] / "runs"


def load(g):
    out = []
    for d in sorted(ROOT.glob(g)):
        if "failed" in d.name or not (d / "results.json").exists():
            continue
        p = json.loads((d / "params.json").read_text()); r = json.loads((d / "results.json").read_text())
        out.append((d.name, p["omega_b"] * 60 / (2 * np.pi), p["theta_max"][0], p.get("n_mix_cycles"), r.get("kLa_1T_10"), r.get("tau_95_qss")))
    return out


c2d = load("c2d_l9_*"); ok = load("bo1_l9_*") + load("bo4_l9_*")
P = np.array([[o[1] / 22.5, o[2] / 5.0] for o in ok])
for c in c2d:
    dist = np.linalg.norm(P - [c[1] / 22.5, c[2] / 5.0], axis=1); i = int(np.argmin(dist)); o = ok[i]
    print(f"{c[0]:26s} rel {c[3]}  kLa {c[4]:6.2f} tau {1000 * c[5]:5.2f} mPa | nearest {o[0]:24s} rel {o[3]} d {dist[i]:.3f} "
          f"kLa {o[4]:6.2f} tau {1000 * o[5]:5.2f} mPa | kLa ratio {c[4] / o[4]:.2f}")
