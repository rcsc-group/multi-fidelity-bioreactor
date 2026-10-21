"""kLa at four C* thresholds (0.10, 0.25, 0.50, 0.75) for every L8 and L9 run with a full C*(t) record, with the SAME whole-period fit as
postprocess.py (_kla_period_at_threshold). Writes experiments/multifidelity/kla_thresholds.csv (level, run, rpm, theta, kLa_10, kLa_25, kLa_50, kLa_75 in 1/h).
Usage: uv run python scripts/kla_thresholds_table.py
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402

THR = (0.10, 0.25, 0.50, 0.75)
SETS = {8: ["t6_l8_*", "kmix_l8_*", "t5_l8_*", "e8_l8_*", "w_l8_*"], 9: ["fig9_l9_rpm*", "e9_l9_rpm22.4_th7", "e9_l9_rpm22.6_th7", "e9_l9_rpm29.9_th7", "e9_l9_rpm30.1_th7"]}


def main():
    rows = []
    for level, pats in SETS.items():
        for pat in pats:
            for d in sorted((ROOT / "runs").glob(pat)):
                if "_ext" in d.name or not (d / "tr_oxy.dat").exists():
                    continue
                try:
                    p = json.loads((d / "params.json").read_text())
                    t, c = pp._compute_c_star(d)
                    T_bio, T_nd = pp._t_scales(p)
                except Exception:
                    continue
                k = [pp._kla_period_at_threshold(t, c, x, T_nd) * 3600.0 / T_bio for x in THR]
                rows.append((level, d.name, p["omega_b"] * 60 / (2 * math.pi), p["theta_max"][0], *k))
    with (ROOT / "experiments/multifidelity/kla_thresholds.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["level", "run", "rpm", "theta", "kLa_10", "kLa_25", "kLa_50", "kLa_75"])
        w.writerows(rows)
    fin = [r for r in rows if all(math.isfinite(x) for x in r[4:])]
    print(f"{len(rows)} runs, {len(fin)} with all four thresholds; L8 {sum(1 for r in fin if r[0] == 8)}, L9 {sum(1 for r in fin if r[0] == 9)}")


if __name__ == "__main__":
    main()
