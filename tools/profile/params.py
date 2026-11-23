"""params.json for a profiling run: the bundle base case (32.5 rpm, theta 7) at a level, for N periods.

    params.py LEVEL PERIODS OUT_JSON [N_MIX_CYCLES]
Diary 2026-11-23 pre-registration (profiling campaign). N_MIX_CYCLES (default 1) sets the tracer/oxygen release.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.bundle import cases as bc  # noqa: E402


def profile_params(level: int, periods: float, n_mix: int = 1) -> dict:
    return bc.base(f"prof_l{level}", level=level) | {"t_end": periods * bc.period(bc.RPM_A), "n_mix_cycles": n_mix}


if __name__ == "__main__":
    lev, per, out = int(sys.argv[1]), float(sys.argv[2]), Path(sys.argv[3])
    n_mix = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    out.write_text(json.dumps(profile_params(lev, per, n_mix), indent=1))
