"""Stacked warm(top)/cold(bottom) video of the actual bag/liquid motion in
the LAB frame (rocking + translation), for visually catching any
unphysical motion in one relative to the other. Reuses render_videos.py's
lab-frame rotation logic (_render_lab/_make_mask), adapted to read
frames_tau's binary format (which carries Th/xh_nd in its header even
when the plain frames/ VOF-only files are empty, e.g. when t_mix >= t_end
gates movies_output() to write nothing).

VOF field only (f), not the tau colormap -- this is about the physical
motion, not the stress field.

Usage:
    uv run python scripts/render_bag_motion_compare.py <warm_dir> <cold_dir> <out.mp4>
"""
import argparse
import json
import math
import struct
import subprocess
import tempfile
from pathlib import Path

import imageio_ffmpeg
import numpy as np
from PIL import Image


def load_frame(path):
    with open(path, "rb") as fh:
        n, = struct.unpack("i", fh.read(4))
        t, = struct.unpack("d", fh.read(8))
        Th, = struct.unpack("d", fh.read(8))
        xh, = struct.unpack("d", fh.read(8))
        f = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
    return n, t, Th, xh, f


def make_mask(n, Ly, n_exp):
    coords = (np.arange(n) + 0.5) / n - 0.5
    X, Y = np.meshgrid(coords, coords)
    if n_exp >= 8.0:
        return np.abs(2 * Y / Ly) <= 1.0
    return (np.abs(2 * X) ** n_exp + np.abs(2 * Y / Ly) ** n_exp) <= 1.0


def render_lab(f, mask, Ly, Th, xh_nd, Th_max, out_w):
    n = f.shape[0]
    liq = np.flipud(f)
    msk = np.flipud(mask)
    # blue = liquid, white = air, both inside the bag; grey = outside
    rgb = np.zeros((n, n, 3), dtype=np.uint8)
    rgb[..., 0] = ((1 - liq) * 255).astype(np.uint8)
    rgb[..., 1] = ((1 - liq) * 255).astype(np.uint8)
    rgb[..., 2] = 255
    rgb[~msk] = 255

    eroded = np.zeros_like(msk)
    eroded[1:-1, 1:-1] = (msk[1:-1, 1:-1] & msk[:-2, 1:-1] & msk[2:, 1:-1]
                          & msk[1:-1, :-2] & msk[1:-1, 2:])
    rgb[msk & ~eroded] = 0

    pad = max(int(n * 0.6), 1)
    ps = n + 2 * pad
    canvas = np.full((ps, ps, 3), 255, dtype=np.uint8)
    canvas[pad:pad + n, pad:pad + n] = rgb
    img = Image.fromarray(canvas)
    img = img.rotate(math.degrees(Th), resample=Image.BICUBIC, expand=False,
                     fillcolor=(255, 255, 255))
    dx = int(xh_nd * n)
    if dx != 0:
        arr = np.array(img)
        shifted = np.full_like(arr, 255)
        if dx > 0:
            shifted[:, dx:] = arr[:, :ps - dx]
        else:
            shifted[:, :ps + dx] = arr[:, -dx:]
        img = Image.fromarray(shifted)

    y_env = abs(math.sin(Th_max)) * 0.5 + abs(math.cos(Th_max)) * Ly / 2
    half_h_lab = y_env * 1.25
    cy = ps // 2
    ch = max(2, (int(half_h_lab * 2 * n) // 2) * 2)
    r_top, r_bot = max(0, cy - ch // 2), min(ps, cy + ch // 2)
    img = img.crop((0, r_top, ps, r_bot))
    w, h = img.size
    new_h = max(2, (round(out_w * h / w) // 2) * 2)
    return img.resize((out_w, new_h), Image.LANCZOS)


def series(run_dir, out_w):
    params = json.loads((run_dir / "params.json").read_text())
    Ly = params["geometry"]["b"] / params["geometry"]["a"]
    n_exp = params["geometry"]["n"]
    Th_max = math.radians(params["theta_max"][0])
    files = sorted((run_dir / "frames_tau").glob("*.bin"))
    mask_cache = {}
    out = []
    for p in files:
        n, t, Th, xh, f = load_frame(p)
        if n not in mask_cache:
            mask_cache[n] = make_mask(n, Ly, n_exp)
        out.append(render_lab(f, mask_cache[n], Ly, Th, xh, Th_max, out_w))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("warm_dir")
    ap.add_argument("cold_dir")
    ap.add_argument("out_mp4")
    ap.add_argument("--width", type=int, default=700)
    ap.add_argument("--fps", type=float, default=10.0)
    args = ap.parse_args()

    warm = series(Path(args.warm_dir), args.width)
    cold = series(Path(args.cold_dir), args.width)
    n = min(len(warm), len(cold))
    print(f"warm: {len(warm)} frames, cold: {len(cold)} frames -> using {n}")

    gap = 8
    frames = []
    for i in range(n):
        w_h = warm[i].size[1]
        canvas = Image.new("RGB", (args.width, w_h + cold[i].size[1] + gap), (0, 0, 0))
        canvas.paste(warm[i], (0, 0))
        canvas.paste(cold[i], (0, w_h + gap))
        frames.append(canvas)

    with tempfile.TemporaryDirectory() as tmp:
        for i, img in enumerate(frames):
            img.save(f"{tmp}/f{i:06d}.png")
        ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
        subprocess.run([ffmpeg, "-y", "-r", str(args.fps), "-i", f"{tmp}/f%06d.png",
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", str(args.out_mp4)],
                       check=True, capture_output=True)
    print(f"wrote {args.out_mp4}")


if __name__ == "__main__":
    main()
