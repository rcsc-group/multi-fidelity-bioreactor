"""Where is the ~32% RMS mismatch actually concentrated? Direct spatial
difference of the late-window-averaged tau field, warm (L6->L7) minus
cold, rather than eyeballing two separately-colored panels at a scale
tuned for the raw magnitude.

Usage:
    uv run python scripts/analyze_cross_level_diff_field.py
"""
import struct
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

RUNS = Path(__file__).parent.parent / "runs"


def load_frame(path):
    with open(path, "rb") as fh:
        n, = struct.unpack("i", fh.read(4))
        t, = struct.unpack("d", fh.read(8))
        fh.read(16)  # Th, xh_nd
        f = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
        tau = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
    return n, t, f, tau


def late_avg(run, n_last=20):
    files = sorted((RUNS / run / "frames_tau").glob("*.bin"))[-n_last:]
    n, _, f0, _ = load_frame(files[0])
    tau_sum = np.zeros((n, n))
    f_sum = np.zeros((n, n))
    for p in files:
        _, _, f, tau = load_frame(p)
        tau_sum += tau
        f_sum += f
    return f_sum / len(files), tau_sum / len(files)


f_warm, tau_warm = late_avg("crosslevel_L6toL7_pilot_vid")
f_cold, tau_cold = late_avg("l7_coldstart_vid")

mask = (f_warm > 0.5) & (f_cold > 0.5)
diff = tau_warm - tau_cold
rms_warm = np.sqrt(np.mean(tau_warm[mask] ** 2))
rms_diff = np.sqrt(np.mean(diff[mask] ** 2))
print(f"liquid-region RMS(diff) / RMS(warm) = {100*rms_diff/rms_warm:.1f}%")

# Where is |diff| concentrated -- interface-adjacent vs bulk?
row_has_liquid = mask.any(axis=1)
rows = np.where(row_has_liquid)[0]
top, bot = rows.min(), rows.max()
band = max(2, (bot - top) // 8)  # top/bottom ~1/8 of the liquid depth
interface_mask = mask.copy()
interface_mask[top + band:bot - band, :] = False
bulk_mask = mask & ~interface_mask

print(f"RMS(diff) near interface (top/bottom {band} rows): "
     f"{np.sqrt(np.mean(diff[interface_mask]**2)):.5f}")
print(f"RMS(diff) in bulk:                                  "
     f"{np.sqrt(np.mean(diff[bulk_mask]**2)):.5f}")

# [FIX] First version resized the FULL (mostly-air, mostly-white) frame to
# a thin 120px strip, and scaled color by the 99th percentile of |diff| --
# a few interface-adjacent outliers stretched the scale so far that the
# actual liquid band (already a small fraction of the frame) compressed to
# a near-invisible sliver. Crop to the liquid band FIRST (like the other
# render script does), THEN scale color from what's actually being shown.
rows_with_liquid = np.where(mask.any(axis=1))[0]
r0, r1 = max(0, rows_with_liquid.min() - 2), min(mask.shape[0], rows_with_liquid.max() + 3)
diff_c, mask_c = diff[r0:r1], mask[r0:r1]

vmax = np.percentile(np.abs(diff_c[mask_c]), 95)  # tighter than 99th -- show the bulk
img = np.clip(diff_c / vmax, -1, 1)
rgb = np.zeros((*img.shape, 3), dtype=np.uint8)
rgb[..., 0] = np.where(img > 0, 255, ((1 + img) * 255)).astype(np.uint8)
rgb[..., 2] = np.where(img < 0, 255, ((1 - img) * 255)).astype(np.uint8)
rgb[..., 1] = (255 * (1 - np.abs(img))).astype(np.uint8)
rgb[~mask_c] = 255
h, w = rgb.shape[:2]
out = Image.fromarray(np.flipud(rgb)).resize((1000, max(60, h * 6)), Image.NEAREST)
out.save("/tmp/xlevel_diff_field.png")
print(f"cropped to liquid band: rows {r0}-{r1} of {mask.shape[0]}, vmax(95th pctile)={vmax:.5f}")
print("wrote /tmp/xlevel_diff_field.png")
