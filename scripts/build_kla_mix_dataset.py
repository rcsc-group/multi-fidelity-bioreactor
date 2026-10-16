"""Collect the kLa and mixing-time values of every multifidelity run into one CSV.

Sets: L8 design (kmix_l8 / t6_l8), L8 T5 grid + pairs (theta 7), L9 (fig9_l9 / x8_l9), L10 (l10_fig13a_xlevel / x8_l10).
Values are read from runs/<run>/results.json (no recomputation).
Usage: uv run python scripts/build_kla_mix_dataset.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.test_h1_mfbml import LF_THETAS, RPMS, hf_points, hf_run, lf_run  # noqa: E402
from scripts.test_h2_mfbml import l9_points, l9_run  # noqa: E402

COLS = ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "kLa_10", "kLa_25", "kLa_50", "dtmix_0.50", "dtmix_0.75", "dtmix_0.95")


def row(level, set_name, rpm, theta, run):
    f = ROOT / "runs" / run / "results.json"
    d = json.loads(f.read_text()) if f.exists() else {}
    return dict(level=level, set=set_name, rpm=rpm, theta=theta, run=run, **{c: d.get(c, "") for c in COLS})


def main():
    rows = []
    for r in RPMS:
        for t in LF_THETAS:
            rows.append(row(8, "design", r, float(t), lf_run(r, t)))
    ids = json.loads((ROOT / "experiments/multifidelity/test_t5_roughness.json").read_text())["run_ids"]
    for r, run in sorted(ids.items(), key=lambda kv: float(kv[0])):
        rows.append(row(8, "t5_grid", float(r), 7.0, run))
    for r, t in l9_points():
        rows.append(row(9, "design", r, t, l9_run(r, t)))
    for r, t in hf_points():
        rows.append(row(10, "design", float(r), t, hf_run(r, t)))
    out = ROOT / "experiments/multifidelity/kla_mix_dataset.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["level", "set", "rpm", "theta", "run", *COLS])
        w.writeheader()
        w.writerows(rows)
    print("wrote", out, len(rows), "rows")
    for lev in (10, 9):
        print(f"\nL{lev}")
        print(f"{'rpm':>5} {'th':>4} {'kLa10':>7} {'kLa25':>7} {'kLa50':>7} {'dt50':>7} {'dt75':>7} {'dt95':>7}")
        for x in rows:
            if x["level"] == lev:
                v = [x[c] for c in ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "dtmix_0.50", "dtmix_0.75", "dtmix_0.95")]
                print(f"{x['rpm']:5g} {x['theta']:4g} " + " ".join(f"{float(a):7.1f}" if a != "" else f"{'--':>7}" for a in v))
    miss = [x["run"] for x in rows if x["kLa_1T_25"] == ""]
    print("\nruns without kLa_1T_25:", miss)


if __name__ == "__main__":
    main()
