"""Stacked warm(top)/cold(bottom) comparison video: shear-stress field,
with live tau_mean and u_rms legends, aligned by CYCLE (not raw t -- the
warm run's t starts at t_checkpoint, the cold run's at 0).

Reads frames_tau/*.bin (format: int32 n, float64 t, Th, xh_nd, then three
n*n float32 buffers f, tau_field [signed], ediss_field). Per frame:
  tau_mean = mean(|tau_field|) over liquid cells (f>0.5) -- reproduces the
             production tau_mean KPI (shear_stress.dat) to within its own
             timestep noise, verified directly (scripts/check_frame_format.py).
  u_rms    = not available per-cell in these frames (only f/tau/ediss are
             dumped); read from the run's own normf.dat
             (sqrt(ux_liq_rms^2+uy_liq_rms^2)) at the nearest output row to
             each frame's t. Labelled u_rms, not u_mean: normf.dat carries
             RMS columns, not a directional mean, and mislabelling it would
             repeat exactly the imprecision flagged earlier this session.

Usage:
    uv run python scripts/render_warmcold_compare.py <warm_run_dir> <cold_run_dir> <out.mp4>
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
from PIL import Image, ImageDraw, ImageFont


def t_per_nd(rpm, theta, L_bio=0.25, H_bio_full=2 * 0.03575):
    w = rpm * 2 * math.pi / 60.0
    T = 2 * math.pi / w
    V = L_bio / 4 * (H_bio_full + 0.5 * L_bio * math.tan(math.radians(theta)))
    U = V / (H_bio_full * 0.5) / T
    return T / (L_bio / U)


def load_frame(path: Path):
    with open(path, "rb") as fh:
        n, = struct.unpack("i", fh.read(4))
        t, = struct.unpack("d", fh.read(8))
        Th, = struct.unpack("d", fh.read(8))
        xh, = struct.unpack("d", fh.read(8))
        f = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
        tau = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
    return n, t, Th, xh, f, tau


def load_u_rms_series(run_dir: Path):
    d = np.loadtxt(run_dir / "normf.dat", skiprows=1)
    t = d[:, 1]
    u_rms = np.sqrt(d[:, 7] ** 2 + d[:, 11] ** 2)  # ux_liq_rms, uy_liq_rms
    return t, u_rms


def _make_mask(n, Ly, n_exp):
    coords = (np.arange(n) + 0.5) / n - 0.5
    X, Y = np.meshgrid(coords, coords)
    if n_exp >= 8.0:
        return np.abs(2 * Y / Ly) <= 1.0
    return (np.abs(2 * X) ** n_exp + np.abs(2 * Y / Ly) ** n_exp) <= 1.0


def _tau_colormap(tau, vmax):
    """Diverging blue-white-red colormap for SIGNED tau, clipped to +-vmax."""
    x = np.clip(tau / vmax, -1.0, 1.0)
    r = np.where(x > 0, 255, ((1 + x) * 255).astype(np.uint8)).astype(np.uint8)
    b = np.where(x < 0, 255, ((1 - x) * 255).astype(np.uint8)).astype(np.uint8)
    g = (255 * (1 - np.abs(x))).astype(np.uint8)
    return np.stack([r, g, b], axis=-1)


def _panel(n, f, tau, mask, Ly, out_w, vmax):
    rgb = _tau_colormap(np.flipud(tau), vmax)
    msk = np.flipud(mask)
    rgb[~msk] = 255
    fmask = np.flipud(f) <= 0.5
    rgb[msk & fmask] = 235  # air, inside bag: light grey (not part of colormap)

    eroded = np.zeros_like(msk)
    eroded[1:-1, 1:-1] = (msk[1:-1, 1:-1] & msk[:-2, 1:-1] & msk[2:, 1:-1]
                          & msk[1:-1, :-2] & msk[1:-1, 2:])
    rgb[msk & ~eroded] = 0

    half_h = Ly / 2 * 1.2
    r0 = max(0, int((0.5 - half_h) * n))
    r1 = min(n, int((0.5 + half_h) * n))
    rgb = rgb[r0:r1]
    img = Image.fromarray(rgb)
    h, w = rgb.shape[:2]
    new_h = max(2, (round(out_w * h / w) // 2) * 2)
    return img.resize((out_w, new_h), Image.LANCZOS)


def _font(size):
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        return ImageFont.load_default()


def _stamp(img, lines, origin):
    """Draw a white-boxed multiline label anchored at (x, y) = origin."""
    draw = ImageDraw.Draw(img)
    font = _font(26)
    text = "\n".join(lines)
    ox, oy = origin
    bbox = draw.multiline_textbbox((ox, oy), text, font=font)
    pad = 4
    x0, y0 = bbox[0] - pad, bbox[1] - pad
    x1, y1 = bbox[2] + pad, bbox[3] + pad
    draw.rectangle([x0, y0, x1, y1], fill=(255, 255, 255))
    draw.multiline_text((ox, oy), text, fill=(0, 0, 0), font=font)
    return img


def build_series(run_dir: Path, rpm: float, theta: float, out_w: int, vmax: float,
                 t0: float = None):
    params = json.loads((run_dir / "params.json").read_text())
    Ly = params["geometry"]["b"] / params["geometry"]["a"]
    n_exp = params["geometry"]["n"]
    T_ND = t_per_nd(rpm, theta)

    frame_files = sorted((run_dir / "frames_tau").glob("frame_*.bin"))
    t_u, u_rms = load_u_rms_series(run_dir)

    rows = []
    for path in frame_files:
        n, t, Th, xh, f, tau = load_frame(path)
        rows.append((t, n, f, tau))
    if t0 is None:
        t0 = rows[0][0]

    mask_cache = {}
    out = []
    for t, n, f, tau in rows:
        cyc = (t - t0) / T_ND
        tau_mean = float(np.abs(tau[f > 0.5]).mean()) if (f > 0.5).any() else 0.0
        u = float(np.interp(t, t_u, u_rms))
        if n not in mask_cache:
            mask_cache[n] = _make_mask(n, Ly, n_exp)
        img = _panel(n, f, tau, mask_cache[n], Ly, out_w, vmax)
        out.append((cyc, tau_mean, u, img))
    return out, t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("warm_dir")
    ap.add_argument("cold_dir")
    ap.add_argument("out_mp4")
    ap.add_argument("--rpm", type=float, required=True)
    ap.add_argument("--theta", type=float, default=7.0)
    ap.add_argument("--width", type=int, default=900)
    args = ap.parse_args()

    warm_dir, cold_dir = Path(args.warm_dir), Path(args.cold_dir)

    # Common color scale from the COLD run's converged (late-window) tau
    # range, so both panels share one colorbar and the visual difference
    # reflects the flow, not a rescaled colormap.
    _, _, _, _, _, tau0 = load_frame(sorted((cold_dir / "frames_tau").glob("*.bin"))[-1])
    vmax = float(np.percentile(np.abs(tau0), 99.5))

    warm, tw0 = build_series(warm_dir, args.rpm, args.theta, args.width, vmax)
    cold, tc0 = build_series(cold_dir, args.rpm, args.theta, args.width, vmax)

    n_frames = min(len(warm), len(cold))
    print(f"warm: {len(warm)} frames, cold: {len(cold)} frames -> using {n_frames}")

    frames = []
    for i in range(n_frames):
        cw, tauw, uw, imgw = warm[i]
        cc, tauc, uc, imgc = cold[i]
        w = imgw.size[0]
        canvas = Image.new("RGB", (w, imgw.size[1] + imgc.size[1] + 6), (0, 0, 0))
        canvas.paste(imgw, (0, 0))
        canvas.paste(imgc, (0, imgw.size[1] + 6))
        canvas = _stamp(canvas, [f"WARM START  cycle {cw:5.1f}",
                                 f"tau_mean = {tauw:.4e}", f"u_rms    = {uw:.4f}"],
                        (8, 8))
        canvas = _stamp(canvas, [f"COLD START  cycle {cc:5.1f}",
                                 f"tau_mean = {tauc:.4e}", f"u_rms    = {uc:.4f}"],
                        (8, imgw.size[1] + 6 + 8))
        frames.append(canvas)

    with tempfile.TemporaryDirectory() as tmp:
        for i, img in enumerate(frames):
            img.save(f"{tmp}/f{i:06d}.png")
        ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
        subprocess.run([ffmpeg, "-y", "-r", "10", "-i", f"{tmp}/f%06d.png",
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", str(args.out_mp4)],
                       check=True)
    print(f"wrote {args.out_mp4}")


if __name__ == "__main__":
    main()
