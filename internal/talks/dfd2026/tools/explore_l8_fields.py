"""Which snapshot frames of rec_l8 hold tracer (plane 7 = c2) and oxygen (plane 9); times in cycles."""
import glob
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO)]
from scripts.postprocess import _t_scales  # noqa: E402

RUN = REPO / "runs/rec_l8_rpm32.5_th7"
p = json.load(open(RUN / "params.json")); _, T = _t_scales(p)
files = sorted(glob.glob(str(RUN / "fields/snap_*.bin")))
for i in (0, 1, 100, 300, 500, 700, 900, 1100, len(files) - 1):
    f = files[i]
    t = float(np.fromfile(f, dtype=np.float64, count=1, offset=8)[0])
    mm = np.memmap(f, dtype=np.float32, mode="r", offset=16, shape=(10, 256, 256))
    print(i, round(t / T, 3), "c2 max", float(mm[7].max()), "c2 min", float(mm[7].min()), "oxy max", float(mm[9].max()))
