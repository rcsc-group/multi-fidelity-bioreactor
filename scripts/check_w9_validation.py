"""Pre-registered warm-L9 validation (diary 2026-10-21, W9): w9_l9 (L8 -> L9 warm start, release after 12 cycles) vs the cold L9 fig9 runs at theta 7.
PASS: kLa_1T_10 and dtmix_0.50 within 15 % at >= 2 of 3 rpm. Reported: dtmix_0.95, kLa_1T_25. Usage: uv run python scripts/check_w9_validation.py"""
import json
from pathlib import Path

R = Path("runs")
QS = ("kLa_1T_10", "dtmix_0.50", "kLa_1T_25", "dtmix_0.95")
ok = {q: 0 for q in QS[:2]}
for rpm in ("17.5", "25", "32.5"):
    try:
        w = json.load(open(R / f"w9_l9_rpm{rpm}_th7/results.json")); c = json.load(open(R / f"fig9_l9_rpm{rpm}/results.json"))
    except Exception as e:  # noqa: BLE001
        print(rpm, "missing", e); continue
    line = []
    for q in QS:
        a, b = w.get(q), c.get(q)
        if a and b and a == a and b == b:
            d = a / b - 1; line.append(f"{q} warm {a:.4g} cold {b:.4g} ({100*d:+.1f} %)")
            if q in ok and abs(d) <= 0.15: ok[q] += 1
        else:
            line.append(f"{q} n/a")
    print(rpm, " | ".join(line))
print("within 15 % count:", ok, "->", "PASS" if all(v >= 2 for v in ok.values()) else "FAIL")
