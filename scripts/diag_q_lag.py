"""Descriptive companion of test Q (NOT pre-registered; added 2026-10-11 after seeing the Q results, disclosed).

Test Q's floor rule raised the PERIODIC bar above the 10 % APERIODIC line at 15-20 rpm, so the pre-registered rule
lost its power there. This check does not need the floor: the phase-matched relative L2 difference e(lag) is
computed at lag = 1, 2 and 3 periods with the same function (scripts/test_q_periodicity.periodicity_error).
The interpolation error is the same at every lag (same fractional frame offset multiplied by the lag, see the
printed offsets), while a drift or chaotic departure grows with the lag. A limit cycle gives e(3T)/e(T) ~ 1.

Usage: uv run python scripts/diag_q_lag.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from test_q_periodicity import load_run, periodicity_error  # noqa: E402



def main():
    runs = sorted((p.name for p in (ROOT / "runs").glob("flow_l8_rpm*_th7")),
                  key=lambda r: float(r.split("rpm")[1].split("_")[0]))
    out = []
    print(f"{'run':22s} {'fpp':>6s} " + " ".join(f"{'e(' + str(k) + 'T)':>9s} {'off':>5s}" for k in (1, 2, 3))
          + f" {'e3/e1':>6s}")
    for r in runs:
        d = load_run(r)
        if d is None:
            print(f"{r:22s} missing")
            continue
        row = {"run": r, "rpm": d["rpm"], "fpp": d["fpp"]}
        cells = []
        for k in (1, 2, 3):
            _, el, _ = periodicity_error(d["t"], d["U"], d["mask"], k * d["T"])
            off = (k * d["fpp"]) % 1.0
            off = min(off, 1 - off)  # fractional frame offset of t + kT, folded to [0, 0.5]
            row[f"e{k}"] = float(np.median(el)) if len(el) else None
            row[f"off{k}"] = off
            cells.append(f"{row[f'e{k}']:9.3e} {off:5.2f}" if row[f"e{k}"] is not None else f"{'n/a':>9s} {off:5.2f}")
        ratio = row["e3"] / row["e1"] if row.get("e3") and row.get("e1") else None
        row["e3_over_e1"] = ratio
        print(f"{r:22s} {d['fpp']:6.1f} " + " ".join(cells) + (f" {ratio:6.2f}" if ratio else "    n/a"))
        out.append(row)
    (ROOT / "experiments/multifidelity/diag_q_lag.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
