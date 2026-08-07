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
    """Must mirror BioReactor.c's solid() branch: at n_exp>=8 the C code
    builds a SHARP rectangle (a_nd=1.0 > domain half-width 0.5, so the
    x-constraint never binds -- only the y-extent/bag height is embedded).
    See scripts/render_videos.py's _make_mask for the full explanation;
    this used to apply the superellipse formula unconditionally, rendering
    rounded corners that don't exist in the real simulated geometry."""
    coords = (np.arange(n) + 0.5) / n - 0.5
    X, Y = np.meshgrid(coords, coords)
    if n_exp >= 8.0:
        return np.abs(2 * Y / Ly) <= 1.0
    return (np.abs(2 * X) ** n_exp + np.abs(2 * Y / Ly) ** n_exp) <= 1.0


def _colormap(tau: np.ndarray, vmin: float, vmax: float, gamma: float = 0.4) -> np.ndarray:
    """Gamma-compressed normalization: shear-stress fields are heavily
    right-skewed (99% of liquid cells sit under the 3rd percentile of the
    range set by a few near-wall peak cells) -- a plain linear norm against
    the true max renders as almost solid black. gamma<1 boosts contrast in
    the low/mid range at the cost of compressing the extreme peak, which is
    fine here since the peak is separately marked with the star."""
    norm = np.clip((tau - vmin) / max(vmax - vmin, 1e-30), 0.0, 1.0) ** gamma
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


def _load_font(size: int):
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        return ImageFont.load_default()


def _make_header(width: int, height: int, text: str) -> Image.Image:
    """Plain white strip with the time label -- rendered in its own space,
    above the bag image, so it never overlaps the field (the previous
    version stamped text directly onto the bag render and it collided with
    the bag's own top edge)."""
    img = Image.new("RGB", (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    font = _load_font(max(14, height // 2))
    bbox = draw.textbbox((0, 0), text, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    draw.text((10, (height - th) / 2 - bbox[1]), text, fill=(0, 0, 0), font=font)
    return img


def _make_colorbar(width: int, height: int, vmin: float, vmax_disp: float,
                    vmax_true: float, gamma: float) -> Image.Image:
    """Real gradient colorbar (not just text): a horizontal strip using the
    exact same colormap/gamma/vmax the field itself is rendered with, plus
    tick marks and value labels at 0 / mid / vmax_disp, and a note that the
    true max (near the star) saturates past the shown scale."""
    bar_h = max(18, height // 3)
    grad = np.linspace(vmin, vmax_disp, width).reshape(1, width)
    grad_rgb = _colormap(grad, vmin, vmax_disp, gamma=gamma)
    grad_rgb = np.repeat(grad_rgb, bar_h, axis=0)

    img = Image.new("RGB", (width, height), (255, 255, 255))
    img.paste(Image.fromarray(grad_rgb), (0, 4))
    draw = ImageDraw.Draw(img)
    draw.rectangle([0, 4, width - 1, 4 + bar_h - 1], outline=(0, 0, 0))

    font = _load_font(max(12, height // 5))
    tick_row_y = 4 + bar_h + 2
    ticks = [vmin, (vmin + vmax_disp) / 2, vmax_disp]
    tick_h = 0
    for tv in ticks:
        x = (tv - vmin) / max(vmax_disp - vmin, 1e-30) * (width - 1)
        draw.line([(x, 4), (x, 4 + bar_h - 1)], fill=(255, 255, 255) if tv not in (vmin, vmax_disp) else (0, 0, 0))
        label = f"{tv:.4f} Pa"
        bbox = draw.textbbox((0, 0), label, font=font)
        tw, tick_h = bbox[2] - bbox[0], bbox[3] - bbox[1]
        tx = min(max(0, x - tw / 2), width - tw)
        draw.text((tx, tick_row_y), label, fill=(0, 0, 0), font=font)

    note_row_y = tick_row_y + tick_h + 4
    note = f"clipped at 99th pct -- true max (star) = {vmax_true:.4f} Pa"
    draw.text((4, note_row_y), note, fill=(80, 80, 80), font=font)
    return img


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

    rgb = _colormap(np.flipud(tau_data) * conv, vmin, vmax, gamma=0.5)
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
    vmin, vmax_true = float(all_liquid_tau.min()), float(all_liquid_tau.max())
    # tau is heavily right-skewed (near-zero through the laminar core, only
    # spiking near the wall boundary layer) -- coloring against the true max
    # renders as almost solid black since ~99% of cells sit far below it.
    # Cap the DISPLAYED color scale at the 99th percentile instead (the true
    # peak, near the star, saturates to the colormap's brightest end, which
    # is itself a reasonable way to flag "this is the hottest region").
    vmax_disp = float(np.percentile(all_liquid_tau, 99))
    print(f"true tau range (Pa): [{vmin:.5f}, {vmax_true:.5f}]")
    print(f"displayed color scale (99th pct cap, Pa): [{vmin:.5f}, {vmax_disp:.5f}]")

    MIN_DURATION_S = 10.0
    realtime_fps = 25.0
    if len(loaded) >= 2:
        dt_phys = (loaded[1][1] - loaded[0][1]) * T_bio
        if dt_phys > 0:
            realtime_fps = 1.0 / dt_phys
    fps = min(realtime_fps, len(loaded) / MIN_DURATION_S)
    fps = max(fps, 1.0)
    print(f"realtime fps would be {realtime_fps:.2f} ({len(loaded)/realtime_fps:.1f}s) -- "
          f"using {fps:.2f} fps instead for a {len(loaded)/fps:.1f}s video")

    out_w = 1200
    header_h = 50
    colorbar_h = 96
    colorbar_img = _make_colorbar(out_w, colorbar_h, vmin, vmax_disp, vmax_true, gamma=0.5)

    frames = []
    for n, t_nd, Th, xh, f_data, tau_data in loaded:
        field_img = _render_body_tau(f_data, tau_data, mask, Ly, out_w, vmin, vmax_disp, t_nd, T_bio, conv)
        header_img = _make_header(out_w, header_h, f"t = {t_nd * T_bio:.2f} s")

        total_h = header_h + field_img.height + colorbar_h
        total_h += total_h % 2  # h264/yuv420p requires even height
        canvas = Image.new("RGB", (out_w, total_h), (255, 255, 255))
        canvas.paste(header_img, (0, 0))
        canvas.paste(field_img, (0, header_h))
        canvas.paste(colorbar_img, (0, header_h + field_img.height))
        frames.append(canvas)

    out = Path(args.out) if args.out else run_dir / "shear_stress_body.mp4"
    print(f"Writing {out} ...")
    _to_mp4(frames, out, fps=fps)
    print("Done.")


if __name__ == "__main__":
    main()
