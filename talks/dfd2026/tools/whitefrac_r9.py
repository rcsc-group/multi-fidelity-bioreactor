"""Fraction of near-white pixels in a slide render (every channel >= 240), whole slide and the lower half.
Usage: uv run python tools/whitefrac_r9.py <png>
"""
import sys
import numpy as np
from PIL import Image

path = sys.argv[1]
im = np.asarray(Image.open(path).convert("RGB")).astype(int)
h, w, _ = im.shape
white = np.all(im >= 240, axis=2)
print(f"{path}: {w}x{h}px, near-white fraction (whole slide) = {100 * white.mean():.1f} %")
print(f"lower half (rows {h // 2}..{h}): near-white fraction = {100 * white[h // 2:].mean():.1f} %")
bands = np.array_split(np.arange(h), 6)
print("by row band (top to bottom, % near-white):", [round(100 * white[b].mean(), 1) for b in bands])
