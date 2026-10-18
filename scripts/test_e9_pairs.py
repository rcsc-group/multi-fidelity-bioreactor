"""Test E9 (pre-registered, diary 2026-10-16): does the +-0.1 rpm jump of dtmix_0.95 seen at L8 persist at L9?

Pairs 22.4/22.6 and 29.9/30.1 rpm, theta 7, level 9 (e9_l9, T5 binary and protocol) vs level 8 (t5_l8).
PERSISTS if |d log dtmix_0.95| >= 0.10 in both pairs; SMOOTHER if <= 0.05 in both; else MIXED. Same rule reported for kLa_1T_10, dtmix_0.50.
Usage: uv run python scripts/test_e9_pairs.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAIRS = (("22.4", "22.6"), ("29.9", "30.1"))
QS = ("dtmix_0.95", "kLa_1T_10", "dtmix_0.50", "kLa_1T_25")


def val(run, q):
    f = ROOT / "runs" / run / "results.json"
    v = json.loads(f.read_text()).get(q, math.nan) if f.exists() else math.nan
    return float(v) if v is not None else math.nan


def rule(d):
    return "PERSISTS" if all(x >= 0.10 for x in d) else "SMOOTHER" if all(x <= 0.05 for x in d) else "MIXED"


def main():
    out = {}
    for q in QS:
        row = {}
        for lev, pre in ((8, "t5_l8_rpm{}_th7"), (9, "e9_l9_rpm{}_th7")):
            vals = [(val(pre.format(a), q), val(pre.format(b), q)) for a, b in PAIRS]
            d = [abs(math.log(x) - math.log(y)) for x, y in vals]
            row[lev] = dict(values=vals, dlog=d)
        out[q] = row
        print(f"{q:11s} L8 d {[round(x, 3) for x in row[8]['dlog']]}  L9 d {[round(x, 3) for x in row[9]['dlog']]}"
              f"  values L9 {[tuple(round(v, 1) for v in p) for p in row[9]['values']]}  L9 rule: "
              f"{rule([x for x in row[9]['dlog'] if math.isfinite(x)]) if all(math.isfinite(x) for x in row[9]['dlog']) else 'incomplete'}")
    (ROOT / "experiments/multifidelity/test_e9_pairs.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
