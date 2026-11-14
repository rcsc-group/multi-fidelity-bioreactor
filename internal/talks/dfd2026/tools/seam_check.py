"""Loop seam check: mean |difference| between consecutive frames, including the wrap frame 95 -> 0. Usage: uv run python seam_check.py"""
import numpy as np
from PIL import Image
fr = [np.asarray(Image.open(f"frames/f_{k:04d}.png").convert("RGB"), dtype=float) for k in range(120)]
d = [np.abs(fr[(k + 1) % 120] - fr[k]).mean() for k in range(120)]
print("median consecutive diff %.3f, max (not wrap) %.3f, wrap 95->0 %.3f" % (np.median(d[:-1]), max(d[:-1]), d[-1]))
