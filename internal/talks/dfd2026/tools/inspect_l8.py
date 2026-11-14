"""Inspect the dense L8 snapshot sets: frames, cycle span, oxygen and tracer ranges at a few frames. Usage: uv run python inspect_l8.py"""
import glob
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO)]
from scripts.postprocess import _t_scales  # noqa: E402

for run in ("rec_l8_rpm32.5_th7", "scal_l8_rpm32.5_th7", "ks_l8_rpm30.1_th7"):
    d = REPO / "runs" / run
    p = json.load(open(d / "params.json")); _, T = _t_scales(p)
    fs = sorted(glob.glob(str(d / "fields/snap_*.bin")))
    t = np.array([float(np.fromfile(f, dtype=np.float64, count=1, offset=8)[0]) for f in fs])
    print(run, "frames", len(fs), "cycles %.2f-%.2f" % (t.min() / T, t.max() / T), "dt/T %.4f" % np.median(np.diff(t) / T),
          "snap", p.get("snap_start_cycle"), p.get("snap_end_cycle"), "n_mix", p.get("n_mix_cycles"), "t_ck/T %.1f" % (p.get("t_checkpoint", 0) / T))
    for i in (0, len(fs) // 3, 2 * len(fs) // 3, len(fs) - 1):
        mm = np.memmap(fs[i], dtype=np.float32, mode="r", offset=16, shape=(10, 256, 256))
        liq = (mm[3] > 0.5) & (mm[4] > 0.5)
        print("   frame %4d cyc %.2f  oxy mean %.3f max %.3f | c mean %.3f | c1 %.3f c3 %.3f" % (i, t[i] / T, mm[9][liq].mean(), mm[9][liq].max(), mm[5][liq].mean(), mm[6][liq].mean(), mm[8][liq].mean()))
