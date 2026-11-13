"""Summarise the control block of experiments/multifidelity/check_replay_kla25_L{LEVEL}.json (read-only).
Usage: uv run python scripts/analyse/check_replay_kla25_control.py 8
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
d = json.loads((ROOT / "experiments" / "multifidelity" / f"check_replay_kla25_L{sys.argv[1]}.json").read_text())
Q = ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "dtmix_0.50")
for block in ("dns_recomputed_vs_stored_pct", "replay_recomputed_vs_stored_pct"):
    worst = {q: max(abs(c[block][q]) for c in d["control"] if c[block][q] == c[block][q]) for q in Q}
    print(block, {q: round(v, 4) for q, v in worst.items()})
print("stored (Oct-1 DB, Oct-31 replay) replay/DNS-1 %:")
for c in d["control"]:
    print(f"{c['rpm']:5g}", "  ".join(f"{q} {c['stored_dev_pct'][q]:+6.1f}" for q in Q))
