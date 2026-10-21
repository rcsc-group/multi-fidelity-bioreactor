"""Control for the looped-replay validation: coupled restart-and-release run from the same checkpoint (scripts/replay_lf.py).
Usage: uv run python scripts/submit_replay_lf_control.py 32.5 7"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import replay_lf as lf  # noqa: E402
from scripts import postprocess as pp  # noqa: E402

rpm, th = float(sys.argv[1]), float(sys.argv[2])
src = lf.SRC_RUNS / lf.source_id(rpm, th)
sp = json.loads((src / "params.json").read_text())
p = lf.control_params(sp, f"lfctl_l8_rpm{rpm:g}_th{th:g}", lf.checkpoint_time(src / "checkpoint.dump"))
print(lf._submit(p, src / "checkpoint.dump", "04:00:00"), p["_binary"], p["t_end"])
