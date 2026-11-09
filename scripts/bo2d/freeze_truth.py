"""Write the frozen truth snapshots (state/truth_snapshot_<qoi>.csv) once, so every benchmark job uses one truth.
Usage: uv run python scripts/bo2d/freeze_truth.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import synth_bench as S  # noqa: E402
for q in ("kLa_1T_10", "dtmix_0.25", "tau_95"):
    S._real_truth(q, 9)
