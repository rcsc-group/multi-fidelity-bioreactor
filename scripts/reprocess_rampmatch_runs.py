"""Recompute results.json for every run that used the ramp-matched binary,
after fixing postprocess.py's QSS-window start (diary.md 2026-09-03): it
was hardcoded to 3 rocking cycles (the fork's own default ramp), but the
ramp-matched binary's actual ramp lasts 8.8-18.7 cycles depending on RPM
(upstream's 30-physical-second linear ramp), so 19-53% of what was being
called "QSS" was still mid-ramp. Overwrites results.json in place for
each run listed.

Usage:
    uv run python scripts/reprocess_rampmatch_runs.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.postprocess import main as postprocess_main

RUNS_DIR = Path(__file__).parent.parent / "runs"

RUN_IDS = (
    [f"fig13a_rampmatch_rpm{rpm:g}" for rpm in [17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5]]
    + [f"fig13a_l6_rpm{rpm:g}" for rpm in [17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5]]
    + ["restart_bias_fresh_l6", "459e1b76", "d74d35b7", "a3070b72", "7beefe12"]
)

for run_id in RUN_IDS:
    run_dir = RUNS_DIR / run_id
    if not run_dir.exists():
        print(f"{run_id}: MISSING, skipped")
        continue
    result = postprocess_main(str(run_dir))
    print(f"{run_id}: tau_100_max={result['tau_100_max']:.6f}  tau_mean_max={result['tau_mean_max']:.6f}")
