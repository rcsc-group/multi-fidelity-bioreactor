"""Step 2 (lead, 2026-11-05), L8 only: tau_95_t of every fidelity-8 run with rpm in [19, 28], listed per theta with the
run id, so that a step or a peak near 23.5 rpm can be read against the neighbours at the same theta. Read-only.
Compute job only."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois  # noqa: E402

rows = []
for d in sorted((ROOT / "runs").glob("*")):
    pf = d / "params.json"
    if not (pf.exists() and (d / "shear_stress.dat").exists()):
        continue
    p = json.loads(pf.read_text())
    if int(p.get("fidelity", -1)) != 8:
        continue
    rpm = round(float(p["omega_b"]) * 60.0 / (2 * np.pi), 3)
    if not (19.0 <= rpm <= 28.0):
        continue
    try:
        t = qois(d.name)["tau_95_t"]
    except Exception:
        continue
    rows.append((round(float(p["theta_max"][0]), 3), rpm, t, d.name))
for th in sorted({r[0] for r in rows}):
    sub = sorted([r for r in rows if r[0] == th], key=lambda r: r[1])
    print(f"theta {th}: n={len(sub)}")
    for r in sub:
        print(f"   rpm {r[1]:7.3f}  tau_95_t {r[2]:.5f}  {r[3]}")
