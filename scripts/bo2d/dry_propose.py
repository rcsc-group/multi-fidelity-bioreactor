"""Dry run of _propose_obj for one arm and objective with an empty in-flight set (nothing saved, nothing submitted).
Usage: uv run python scripts/bo2d/dry_propose.py <arm> <kla|mix>"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import loop as L  # noqa: E402

arm, obj = sys.argv[1], sys.argv[2]
cfg = json.loads(L.CFG.read_text())
tasks = []
L._propose_obj(arm, obj, 2, np.zeros((0, 2)), int(cfg["q"]), np.log(cfg["cap"]), tasks)
print("rate units: 1/s" if obj == "mix" else "", [t["x"] for t in tasks])
