"""Crop Kim 2025 Fig. (p.7, img i-007-006) to the field of panel (a) only: no axes, no colour bar, no labels.
Usage: uv run python tools/crop_s3_kim.py SRC.png -> ../public/s3_kim.png
Fractions found by eye from the full image; re-check the output visually.
"""
import sys
from pathlib import Path

from PIL import Image

src = Image.open(sys.argv[1]).convert("RGB")
w, h = src.size
box = (int(0.135 * w), int(0.045 * h), int(0.915 * w), int(0.415 * h))
out = Path(__file__).resolve().parent.parent / "public" / "s3_kim.png"
src.crop(box).save(out)
print("wrote", out, src.crop(box).size, "from", (w, h))
