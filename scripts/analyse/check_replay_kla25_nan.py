"""List NaN QoIs and window spans per rpm for an experiments/multifidelity/check_replay_kla25_L{LEVEL}.json (read-only).
Usage: uv run python scripts/analyse/check_replay_kla25_nan.py 9
"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
d = json.loads((ROOT / "experiments" / "multifidelity" / f"check_replay_kla25_L{sys.argv[1]}.json").read_text())
for r in d["rows"]:
    bad = [f"{s}:{q}" for s in ("dns", "replay") for q in ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "dtmix_0.50")
           if r[s][q] is None or (isinstance(r[s][q], float) and math.isnan(r[s][q]))]
    print(f"{r['rpm']:5g}", "NaN:", bad or "none", "| replay span", round(r["replay"]["span_cycles"], 1),
          "dns span", round(r["dns"]["span_cycles"], 1), "| replay c*_last", round(r["replay"]["c_star_last"], 3))
