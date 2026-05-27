"""Render VOF videos from binary frame dumps produced by BioReactor-video.

Binary frame format (per file frames/frame_XXXXXX.bin):
  int32    n        grid size (n×n uniform grid)
  float64  t_nd     non-dim simulation time
  float64  Th       current tilt angle (radians)
  float64  xh_nd    horizontal displacement (non-dim, lab frame)
  float32  [n*n]    VOF field f, row-major, row 0 = bottom (y=Y0)

Usage:
  python scripts/render_videos.py <run_dir>

Outputs:
  <run_dir>/volume_fraction.mp4      body-frame video
  <run_dir>/volume_fraction_lab.mp4  lab-frame video
Cleans up <run_dir>/frames/ when done.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


# ── geometry ─────────────────────────────────────────────────────────────────

def compute_T_bio(params: dict) -> float:
    a  = params["geometry"]["a"]
    b  = params["geometry"]["b"]
    th = math.radians(params["theta_max"][0])
    T_per = 2 * math.pi / params["omega_b"]
    V     = a / 4 * (b + 0.5 * a * math.tan(th))
    U     = V / (b * 0.5) / T_per
    return a / U


# ── I/O ──────────────────────────────────────────────────────────────────────

def load_frame(path: Path) -> tuple[int, float, float, float, np.ndarray]:
    """Returns (n, t_nd, Th, xh_nd, field[n,n]) — row 0 = bottom of domain."""
    with open(path, "rb") as fh:
        (n,)   = struct.unpack("i", fh.read(4))
        (t,)   = struct.unpack("d", fh.read(8))
        (Th,)  = struct.unpack("d", fh.read(8))
        (xh,)  = struct.unpack("d", fh.read(8))
        data   = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
    return n, t, Th, xh, data


# ── rendering ────────────────────────────────────────────────────────────────

def _make_mask(n: int, Ly: float, n_exp: float) -> np.ndarray:
    """Boolean (n,n) mask — True inside superellipse. Row 0 = bottom (y=-0.5)."""
    coords = (np.arange(n) + 0.5) / n - 0.5   # [-0.5, 0.5)
    X, Y   = np.meshgrid(coords, coords)        # row = y index
    return (np.abs(2 * X) ** n_exp + np.abs(2 * Y / Ly) ** n_exp) <= 1.0


def _bwr(f: np.ndarray) -> np.ndarray:
    """Map f∈[0,1] → uint8 RGB using blue(0)→white(0.5)→red(1) colormap."""
    f = np.clip(f, 0.0, 1.0)
    r = np.where(f < 0.5, (f * 2 * 255).astype(np.uint8), np.uint8(255))
    g = np.where(f < 0.5, (f * 2 * 255).astype(np.uint8),
                 ((1.0 - f) * 2 * 255).astype(np.uint8))
    b = np.where(f < 0.5, np.uint8(255), ((1.0 - f) * 2 * 255).astype(np.uint8))
    return np.stack([r, g, b], axis=-1)


def _draw_label(img: Image.Image, text: str) -> Image.Image:
    """Stamp a physical-time label in the top-left corner."""
    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.load_default(size=28)
    except TypeError:
        font = ImageFont.load_default()
    bbox = draw.textbbox((8, 8), text, font=font)
    draw.rectangle([bbox[0] - 3, bbox[1] - 3, bbox[2] + 3, bbox[3] + 3],
                   fill=(255, 255, 255))
    draw.text((8, 8), text, fill=(0, 0, 0), font=font)
    return img


def _render_body(data: np.ndarray, mask: np.ndarray,
                 Ly: float, out_w: int, t_nd: float, T_bio: float) -> Image.Image:
    """Body-frame image cropped to bag + 20 % margin."""
    n   = data.shape[0]
    rgb = _bwr(np.flipud(data))          # row 0 = top (y=+0.5)
    msk = np.flipud(mask)
    rgb[~msk] = 255                      # exterior → white

    # Draw bag boundary: set pixels at mask edge to black
    eroded = np.zeros_like(msk)
    eroded[1:-1, 1:-1] = (msk[1:-1, 1:-1] &
                           msk[:-2, 1:-1] & msk[2:, 1:-1] &
                           msk[1:-1, :-2] & msk[1:-1, 2:])
    boundary = msk & ~eroded
    rgb[boundary] = 0

    # Vertical crop: show Ly*1.2 centred on y=0
    half_h = Ly / 2 * 1.2
    row_top = max(0, int((0.5 - half_h) * n))
    row_bot = min(n, int((0.5 + half_h) * n))
    rgb = rgb[row_top:row_bot, :, :]

    img = Image.fromarray(rgb)
    h, w = rgb.shape[:2]
    new_h = max(2, (round(out_w * h / w) // 2) * 2)
    return img.resize((out_w, new_h), Image.LANCZOS)


def _render_lab(data: np.ndarray, mask: np.ndarray,
                Ly: float, Th: float, xh_nd: float,
                Th_max: float, out_w: int) -> Image.Image:
    """Lab-frame image: body frame rotated by Th and translated by xh_nd."""
    n   = data.shape[0]
    rgb = _bwr(np.flipud(data))
    msk = np.flipud(mask)
    rgb[~msk] = 255

    # Boundary
    eroded = np.zeros_like(msk)
    eroded[1:-1, 1:-1] = (msk[1:-1, 1:-1] &
                           msk[:-2, 1:-1] & msk[2:, 1:-1] &
                           msk[1:-1, :-2] & msk[1:-1, 2:])
    rgb[msk & ~eroded] = 0

    # Pad so rotation doesn't clip corners
    pad = max(int(n * 0.6), 1)
    ps  = n + 2 * pad
    canvas = np.full((ps, ps, 3), 255, dtype=np.uint8)
    canvas[pad:pad + n, pad:pad + n] = rgb

    img = Image.fromarray(canvas)
    # Rotate: positive Th → bag tilts counterclockwise in lab frame
    img = img.rotate(math.degrees(Th), resample=Image.BICUBIC,
                     expand=False, fillcolor=(255, 255, 255))

    # Horizontal translation: xh_nd × n pixels
    dx = int(xh_nd * n)
    if dx != 0:
        arr     = np.array(img)
        shifted = np.full_like(arr, 255)
        if dx > 0:
            shifted[:, dx:] = arr[:, :ps - dx]
        else:
            shifted[:, :ps + dx] = arr[:, -dx:]
        img = Image.fromarray(shifted)

    # Crop to rotated-bag envelope + 20 % margin
    y_env = abs(math.sin(Th_max)) * 0.5 + abs(math.cos(Th_max)) * Ly / 2
    half_h_lab = y_env * 1.2
    cy    = ps // 2
    ch    = max(2, (int(half_h_lab * 2 * n) // 2) * 2)
    r_top = max(0, cy - ch // 2)
    r_bot = min(ps, cy + ch // 2)
    img   = img.crop((0, r_top, ps, r_bot))

    w, h  = img.size
    new_h = max(2, (round(out_w * h / w) // 2) * 2)
    return img.resize((out_w, new_h), Image.LANCZOS)


# ── video assembly ────────────────────────────────────────────────────────────

def _to_mp4(frames: list[Image.Image], out: Path, fps: float = 25.0) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        for i, img in enumerate(frames):
            img.save(f"{tmp}/f{i:06d}.png")
        first = frames[0]
        print(f"  frame size (WxH): {first.size[0]}x{first.size[1]}, "
              f"total frames: {len(frames)}")
        subprocess.run(
            ["ffmpeg", "-y", "-r", str(fps),
             "-i", f"{tmp}/f%06d.png",
             "-c:v", "libx264", "-pix_fmt", "yuv420p", str(out)],
            check=True,
        )


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    args = ap.parse_args()
    run_dir = Path(args.run_dir)

    params  = json.loads((run_dir / "params.json").read_text())
    L_bio   = params["geometry"]["a"]
    Ly      = params["geometry"]["b"] / L_bio
    n_exp   = params["geometry"]["n"]
    Th_max  = math.radians(params["theta_max"][0])
    T_bio   = compute_T_bio(params)

    frame_files = sorted((run_dir / "frames").glob("frame_*.bin"))
    if not frame_files:
        print(f"No frame files in {run_dir}/frames/ — nothing to render.")
        return

    # Realtime fps: 1 physical second of video = 1 physical second of simulation
    if len(frame_files) >= 2:
        _, t0_nd, _, _, _ = load_frame(frame_files[0])
        _, t1_nd, _, _, _ = load_frame(frame_files[1])
        dt_phys = (t1_nd - t0_nd) * T_bio
        fps = 1.0 / dt_phys if dt_phys > 0 else 25.0
    else:
        fps = 25.0
    print(f"  realtime fps: {fps:.2f}")

    print(f"Rendering {len(frame_files)} frames …")
    body_frames: list[Image.Image] = []
    lab_frames:  list[Image.Image] = []

    for path in frame_files:
        n, t_nd, Th, xh_nd, data = load_frame(path)
        mask  = _make_mask(n, Ly, n_exp)
        label = f"t = {t_nd * T_bio:.2f} s"
        body_frames.append(_draw_label(
            _render_body(data, mask, Ly, 1200, t_nd, T_bio), label))
        lab_frames.append(_draw_label(
            _render_lab(data, mask, Ly, Th, xh_nd, Th_max, 1200), label))

    print("Writing volume_fraction.mp4 …")
    _to_mp4(body_frames, run_dir / "volume_fraction.mp4", fps=fps)

    print("Writing volume_fraction_lab.mp4 …")
    _to_mp4(lab_frames, run_dir / "volume_fraction_lab.mp4", fps=fps)

    shutil.rmtree(run_dir / "frames")
    print("Done.")


if __name__ == "__main__":
    main()
