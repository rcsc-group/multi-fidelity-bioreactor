"""Render a shear-stress-field video (body/bag frame) from frames_tau/*.bin,
produced by BioReactor's movies_output_tau event (VIDEOS=1 build).

Binary frame format (per file frames_tau/frame_XXXXXX.bin):
  int32    n        grid size (n×n uniform grid)
  float64  t_nd     non-dim simulation time
  float64  Th       current tilt angle (radians)
  float64  xh_nd    horizontal displacement (non-dim, lab frame)
  float32  [n*n]    VOF field f, row-major, row 0 = bottom (y=Y0)
  float32  [n*n]    tau field (Pa-equivalent after *conv), same layout

Body frame only (non-inertial/bag frame) per user request: no rotation/
translation applied, matches _render_body() in scripts/render_videos.py.
Fixed color-axis limits across the whole video (global min/max over all
frames' liquid cells) -- no legend/colorbar, just the numeric range stamped
in the corner alongside the time label. A star marks the per-frame argmax
(tau_max) location.
"""
import argparse
import json
import math
import struct
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

try:
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    _HAVE_MPL = True
except ImportError:
    _HAVE_MPL = False


def compute_T_bio(params: dict) -> float:
    a = params["geometry"]["a"]
    b = params["geometry"]["b"]
    th = math.radians(params["theta_max"][0])
    T_per = 2 * math.pi / params["omega_b"]
    V = a / 4 * (b + 0.5 * a * math.tan(th))
    U = V / (b * 0.5) / T_per
    return a / U


def load_frame(path: Path):
    with open(path, "rb") as fh:
        (n,) = struct.unpack("i", fh.read(4))
        (t,) = struct.unpack("d", fh.read(8))
        (Th,) = struct.unpack("d", fh.read(8))
        (xh,) = struct.unpack("d", fh.read(8))
        f_data = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
        tau_data = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
    return n, t, Th, xh, f_data, tau_data


def _make_mask(n: int, Ly: float, n_exp: float) -> np.ndarray:
    coords = (np.arange(n) + 0.5) / n - 0.5
    X, Y = np.meshgrid(coords, coords)
    return (np.abs(2 * X) ** n_exp + np.abs(2 * Y / Ly) ** n_exp) <= 1.0


def _colormap(tau: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    norm = np.clip((tau - vmin) / max(vmax - vmin, 1e-30), 0.0, 1.0)
    if _HAVE_MPL:
        rgb = (cm.get_cmap("inferno")(norm)[..., :3] * 255).astype(np.uint8)
    else:
        # hand-rolled black -> red -> yellow -> white fallback
        r = np.clip(norm * 3.0, 0, 1)
        g = np.clip(norm * 3.0 - 1.0, 0, 1)
        b = np.clip(norm * 3.0 - 2.0, 0, 1)
        rgb = (np.stack([r, g, b], axis=-1) * 255).astype(np.uint8)
    return rgb


def _star_polygon(cx: float, cy: float, r_out: float, r_in: float, n_points: int = 5):
    pts = []
    for k in range(2 * n_points):
        ang = math.pi / 2 + k * math.pi / n_points
        r = r_out if k % 2 == 0 else r_in
        pts.append((cx + r * math.cos(ang), cy - r * math.sin(ang)))
    return pts


def _draw_label(img: Image.Image, lines: list[str]) -> Image.Image:
    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.load_default(size=26)
    except TypeError:
        font = ImageFont.load_default()
    text = "\n".join(lines)
    bbox = draw.multiline_textbbox((8, 8), text, font=font)
    draw.rectangle([bbox[0] - 4, bbox[1] - 4, bbox[2] + 4, bbox[3] + 4], fill=(255, 255, 255))
    draw.multiline_text((8, 8), text, fill=(0, 0, 0), font=font)
    return img


def _render_body_tau(f_data, tau_data, mask, Ly, out_w, vmin, vmax, t_nd, T_bio, conv):
    n = f_data.shape[0]
    liquid = f_data > 0.5
    valid = liquid & mask

    rgb = _colormap(np.flipud(tau_data) * conv, vmin, vmax)
    msk_disp = np.flipud(mask)
    valid_disp = np.flipud(valid)
    rgb[~msk_disp] = 255
    rgb[msk_disp & ~valid_disp] = 255  # gas cells inside bag -> white

    eroded = np.zeros_like(msk_disp)
    eroded[1:-1, 1:-1] = (msk_disp[1:-1, 1:-1] &
                           msk_disp[:-2, 1:-1] & msk_disp[2:, 1:-1] &
                           msk_disp[1:-1, :-2] & msk_disp[1:-1, 2:])
    boundary = msk_disp & ~eroded
    rgb[boundary] = 0

    half_h = Ly / 2 * 1.2
    row_top = max(0, int((0.5 - half_h) * n))
    row_bot = min(n, int((0.5 + half_h) * n))
    rgb = rgb[row_top:row_bot, :, :]

    img = Image.fromarray(rgb)
    h, w = rgb.shape[:2]
    new_h = max(2, (round(out_w * h / w) // 2) * 2)
    img = img.resize((out_w, new_h), Image.LANCZOS)

    # star at argmax(tau) within valid cells (original, non-flipped indices)
    star_xy = None
    if valid.any():
        tau_masked = np.where(valid, tau_data, -np.inf)
        row0, col0 = np.unravel_index(np.argmax(tau_masked), tau_masked.shape)
        row_disp = (n - 1 - row0) - row_top   # flipud + crop
        col_disp = col0
        if 0 <= row_disp < (row_bot - row_top):
            px = col_disp / n * out_w
            py = row_disp / (row_bot - row_top) * new_h
            star_xy = (px, py)

    if star_xy is not None:
        draw = ImageDraw.Draw(img)
        r_out = out_w * 0.02
        pts = _star_polygon(star_xy[0], star_xy[1], r_out, r_out * 0.42)
        draw.polygon(pts, fill=(0, 255, 0), outline=(0, 0, 0))

    return img


def _to_mp4(frames, out: Path, fps: float = 25.0) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        for i, img in enumerate(frames):
            img.save(f"{tmp}/f{i:06d}.png")
        print(f"  frame size: {frames[0].size}, total frames: {len(frames)}")
        subprocess.run(
            ["ffmpeg", "-y", "-r", str(fps), "-i", f"{tmp}/f%06d.png",
             "-c:v", "libx264", "-pix_fmt", "yuv420p", str(out)],
            check=True,
        )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    run_dir = Path(args.run_dir)

    params = json.loads((run_dir / "params.json").read_text())
    L_bio = params["geometry"]["a"]
    Ly = params["geometry"]["b"] / L_bio
    n_exp = params["geometry"]["n"]
    T_bio = compute_T_bio(params)
    rho_w = 1000.0
    U_bio = L_bio / T_bio
    conv = rho_w * U_bio ** 2   # nondim tau -> Pa

    frame_files = sorted((run_dir / "frames_tau").glob("frame_*.bin"))
    if not frame_files:
        print(f"No frame files in {run_dir}/frames_tau/ — nothing to render.")
        return

    loaded = [load_frame(p) for p in frame_files]
    n0, _, _, _, _, _ = loaded[0]
    mask = _make_mask(n0, Ly, n_exp)

    # fixed color-axis limits: global min/max over all frames' liquid cells (Pa)
    all_liquid_tau = []
    for n, t_nd, Th, xh, f_data, tau_data in loaded:
        valid = (f_data > 0.5) & mask
        if valid.any():
            all_liquid_tau.append(tau_data[valid] * conv)
    all_liquid_tau = np.concatenate(all_liquid_tau)
    vmin, vmax = float(all_liquid_tau.min()), float(all_liquid_tau.max())
    print(f"color-axis limits (Pa): [{vmin:.5f}, {vmax:.5f}]")

    if len(loaded) >= 2:
        dt_phys = (loaded[1][1] - loaded[0][1]) * T_bio
        fps = 1.0 / dt_phys if dt_phys > 0 else 25.0
    else:
        fps = 25.0
    print(f"realtime fps: {fps:.2f}")

    frames = []
    for n, t_nd, Th, xh, f_data, tau_data in loaded:
        img = _render_body_tau(f_data, tau_data, mask, Ly, 1200, vmin, vmax, t_nd, T_bio, conv)
        label = [f"t = {t_nd * T_bio:.2f} s", f"tau range: [{vmin:.4f}, {vmax:.4f}] Pa"]
        frames.append(_draw_label(img, label))

    out = Path(args.out) if args.out else run_dir / "shear_stress_body.mp4"
    print(f"Writing {out} ...")
    _to_mp4(frames, out, fps=fps)
    print("Done.")


if __name__ == "__main__":
    main()
