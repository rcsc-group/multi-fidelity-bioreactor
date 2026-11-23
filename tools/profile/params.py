"""params.json for a profiling run: the bundle base case (32.5 rpm, theta 7) at a level, for N periods.

    params.py LEVEL PERIODS OUT_JSON
Diary 2026-11-23 pre-registration (profiling campaign).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.bundle import cases as bc  # noqa: E402


def profile_params(level: int, periods: float) -> dict:
    return bc.base(f"prof_l{level}", level=level) | {"t_end": periods * bc.period(bc.RPM_A)}


if __name__ == "__main__":
    lev, per, out = int(sys.argv[1]), float(sys.argv[2]), Path(sys.argv[3])
    out.write_text(json.dumps(profile_params(lev, per), indent=1))
