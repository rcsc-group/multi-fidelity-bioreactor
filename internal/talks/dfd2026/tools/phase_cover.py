"""Phase coverage of the L10 snapshot set: phase = (t / T) mod 1 of every frame; gaps between sorted phases. Usage: uv run python phase_cover.py"""
import glob
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO)]
from scripts.postprocess import _t_scales  # noqa: E402

for run in ("v_l10_rpm32.5_th7", "scal_l9_rpm32.5_th7"):
    d = REPO / "runs" / run
    p = json.load(open(d / "params.json")); _, T = _t_scales(p)
    fs = sorted(glob.glob(str(d / "fields/snap_*.bin")))
    t = np.array([float(np.fromfile(f, dtype=np.float64, count=1, offset=8)[0]) for f in fs])
    c = t / T; sel = c >= c.min() + 3
    ph = np.sort(np.mod(c[sel], 1.0)); gaps = np.diff(np.r_[ph, ph[0] + 1])
    print(run, "frames", sel.sum(), "cycles %.1f-%.1f" % (c[sel].min(), c[sel].max()), "dt/T %.4f" % np.median(np.diff(c)),
          "max phase gap %.4f (%.1f frames of 96)" % (gaps.max(), gaps.max() * 96), "median gap %.4f" % np.median(gaps))
