"""L9 protocol check (pre-registered, diary 2026-10-19): same-protocol L9 runs e9_l9_rpm{17.5, 25, 32.5}_th7 vs the older cold runs fig9_l9_rpm*.
PROTOCOL-NEUTRAL if kLa_1T_10 and dtmix_0.50 agree within 15% at >= 2 of 3 rpm.
Usage: uv run python scripts/test_l9_protocol.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def val(run, q):
    f = ROOT / "runs" / run / "results.json"
    v = json.loads(f.read_text()).get(q) if f.exists() else None
    return float(v) if v else np.nan


def main():
    ok = {q: 0 for q in ("kLa_1T_10", "dtmix_0.50")}
    for rpm in (17.5, 25, 32.5):
        line = f"{rpm:5g} rpm: "
        for q in ok:
            new, old = val(f"e9_l9_rpm{rpm:g}_th7", q), val(f"fig9_l9_rpm{rpm:g}", q)
            d = new / old - 1
            ok[q] += abs(d) <= 0.15
            line += f"{q} new {new:7.2f} old {old:7.2f} ({100 * d:+5.1f}%)   "
        for q in ("kLa_1T_25", "dtmix_0.95"):
            new, old = val(f"e9_l9_rpm{rpm:g}_th7", q), val(f"fig9_l9_rpm{rpm:g}", q)
            line += f"{q} {new:7.1f} vs {old:7.1f}   "
        print(line)
    verdict = "PROTOCOL-NEUTRAL" if all(v >= 2 for v in ok.values()) else "NOT NEUTRAL"
    print("agree within 15% at (of 3 rpm):", ok, "->", verdict)


if __name__ == "__main__":
    main()
