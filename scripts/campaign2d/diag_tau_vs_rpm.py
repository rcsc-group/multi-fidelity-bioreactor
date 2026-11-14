"""Step 2 of the sloshing test (lead, 2026-11-05): tau_95 against rpm for every run folder, by theta group.
Read-only. Uses the hydro window of scripts/hydro_dataset_v2.qois (the database definition), so all runs are on one
definition. Prints, per theta group and per 1-rpm bin, the median tau_95 (ln sd and count), and the step/peak near the
sloshing candidates 23.5 rpm (n=3 tongue) and 35.2 rpm (primary mode, campaign.py-independent, see diary).
Compute job only.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois  # noqa: E402

GROUPS = [("theta 2-3", 2.0, 3.0), ("theta 4", 3.9, 4.1), ("theta 5.5", 5.4, 5.6), ("theta 7", 6.9, 7.1), ("theta 9", 8.9, 9.1)]


def main():
    rows = []
    for d in sorted((ROOT / "runs").glob("*")):
        pf = d / "params.json"
        if not (pf.exists() and (d / "shear_stress.dat").exists()):
            continue
        try:
            p = json.loads(pf.read_text())
            q = qois(d.name)
        except Exception:
            continue
        rpm = round(float(p["omega_b"]) * 60.0 / (2 * np.pi), 3)
        rows.append((float(p["theta_max"][0]), rpm, int(p.get("fidelity", -1)), q["tau_95_t"], d.name))
    print(f"runs with a tau_95 window: {len(rows)}")
    for name, lo, hi in GROUPS:
        sub = [r for r in rows if lo <= r[0] <= hi and r[3] > 0]
        print(f"\n== {name}: n={len(sub)}")
        if not sub:
            continue
        rp = np.array([r[1] for r in sub])
        tv = np.array([r[3] for r in sub])
        for b0 in np.arange(15.0, 38.0, 1.0):
            m = (rp >= b0) & (rp < b0 + 1.0)
            if m.sum() == 0:
                continue
            print(f"   rpm {b0:5.1f}-{b0 + 1:5.1f}: n={m.sum():3d}  median tau_95_t = {np.median(tv[m]):.5f}  "
                  f"ln sd = {np.std(np.log(tv[m]), ddof=1) if m.sum() > 1 else float('nan'):.3f}")


if __name__ == "__main__":
    main()
