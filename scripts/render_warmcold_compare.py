"""Stacked warm(top)/cold(bottom) comparison video: signed shear-stress
field in Pa (with a labeled colorbar), live tau_mean/u_rms legend, and a
persistent caption stating the warm-start's actual initial condition --
aligned by CYCLE (not raw t, since the warm run's t starts at
t_checkpoint, the cold run's at 0).

Reads frames_tau/*.bin (int32 n, float64 t/Th/xh_nd, then three n*n
float32 buffers f, tau_field [signed, nondimensional], ediss_field).
tau_field and tau_mean are converted to Pa via tau_Pa = tau_nd * rho_w *
U_bio^2 (rho_w=1000 kg/m^3), matching scripts/postprocess.py's own
convention -- the project's y-axis in every other tau figure (e.g.
plot_fig13a_current.py) is already in Pa, and the raw nondimensional
values mean nothing to a reader without that context.

u_rms is read from the run's own normf.dat (sqrt(ux_liq_rms^2 +
uy_liq_rms^2), nearest output row to each frame's t) and left
NONdimensional -- normf.dat carries RMS columns, not a directional mean,
so it is labelled u_rms rather than u_mean.

Usage:
    uv run python scripts/render_warmcold_compare.py <warm_dir> <cold_dir> <out.mp4> \\
        --rpm 22.5 --source-rpm 17.5
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

RHO_W = 1.0e3  # kg/m^3


def t_per_nd(rpm, theta, L_bio=0.25, H_bio_full=2 * 0.03575):
    w = rpm * 2 * math.pi / 60.0
    T = 2 * math.pi / w
    V = L_bio / 4 * (H_bio_full + 0.5 * L_bio * math.tan(math.radians(theta)))
    U = V / (H_bio_full * 0.5) / T
    return T / (L_bio / U), L_bio / (V / (H_bio_full * 0.5) / T)  # (T_per_nd, T_bio)


def U_bio_of(rpm, theta, L_bio=0.25, H_bio_full=2 * 0.03575):
    w = rpm * 2 * math.pi / 60.0
    T = 2 * math.pi / w
    V = L_bio / 4 * (H_bio_full + 0.5 * L_bio * math.tan(math.radians(theta)))
    return V / (H_bio_full * 0.5) / T


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


def _diverging_rgb(x):
    """x in [-1, 1] -> blue(-1) - white(0) - red(+1), vectorized uint8 RGB."""
    x = np.clip(x, -1.0, 1.0)
    r = np.where(x > 0, 255, ((1 + x) * 255)).astype(np.uint8)
    b = np.where(x < 0, 255, ((1 - x) * 255)).astype(np.uint8)
    g = (255 * (1 - np.abs(x))).astype(np.uint8)
    return np.stack([r, g, b], axis=-1)


def _panel(f, tau_pa, mask, Ly, out_w, vmax_pa):
    n = f.shape[0]
    rgb = _diverging_rgb(np.flipud(tau_pa) / vmax_pa)
    msk = np.flipud(mask)
    rgb[~msk] = 255
    fmask = np.flipud(f) <= 0.5
    rgb[msk & fmask] = 232  # air inside the bag: light grey, outside the colormap

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


def _colorbar(total_height, vmax_pa, width=46):
    """Vertical colorbar, red at top (+vmax), blue at bottom (-vmax).

    The gradient itself is inset by `pad` on top and bottom so the tick
    labels at its two ends have room to draw without being clipped by the
    image bounds (the naive version drew text centred on y=0 and y=height-1,
    which is exactly where PIL silently clips it)."""
    pad = 12
    bar_h = total_height - 2 * pad
    grad = np.linspace(1, -1, bar_h).reshape(-1, 1)
    bar = np.repeat(_diverging_rgb(grad), width, axis=1)
    img = Image.new("RGB", (width + 78, total_height), (255, 255, 255))
    img.paste(Image.fromarray(bar), (0, pad))
    draw = ImageDraw.Draw(img)
    font = _font(20)
    for frac, val in [(0.0, vmax_pa), (0.5, 0.0), (1.0, -vmax_pa)]:
        y = pad + int(frac * (bar_h - 1))
        y = min(max(y, pad + 10), total_height - 11)  # keep text on-canvas
        draw.line([(width, pad + int(frac * (bar_h - 1))), (width + 6, pad + int(frac * (bar_h - 1)))],
                  fill=(0, 0, 0), width=2)
        draw.text((width + 10, y - 10), f"{val:+.2g}", fill=(0, 0, 0), font=font)
    draw.text((width + 10, total_height // 2 - 32), "Pa", fill=(0, 0, 0), font=font)
    return img


def _label_block(lines, font_size=20):
    """Small, semi-opaque label as its own image (composited, not stamped
    directly on the flow field, so it never obscures more than its own box)."""
    font = _font(font_size)
    tmp = Image.new("RGB", (10, 10))
    d = ImageDraw.Draw(tmp)
    text = "\n".join(lines)
    bbox = d.multiline_textbbox((0, 0), text, font=font)
    w, h = bbox[2] - bbox[0] + 12, bbox[3] - bbox[1] + 10
    img = Image.new("RGBA", (w, h), (255, 255, 255, 200))
    ImageDraw.Draw(img).multiline_text((6, 5), text, fill=(0, 0, 0), font=font)
    return img


def build_series(run_dir: Path, rpm: float, theta: float, out_w: int, vmax_pa: float,
                 tau_scale: float):
    params = json.loads((run_dir / "params.json").read_text())
    Ly = params["geometry"]["b"] / params["geometry"]["a"]
    n_exp = params["geometry"]["n"]
    T_ND, _ = t_per_nd(rpm, theta)

    frame_files = sorted((run_dir / "frames_tau").glob("frame_*.bin"))
    t_u, u_rms = load_u_rms_series(run_dir)
    t0 = None
    mask_cache = {}
    out = []
    for path in frame_files:
        n, t, Th, xh, f, tau = load_frame(path)
        if t0 is None:
            t0 = t
        cyc = (t - t0) / T_ND
        tau_pa = tau * tau_scale
        tau_mean_pa = float(np.abs(tau_pa[f > 0.5]).mean()) if (f > 0.5).any() else 0.0
        u = float(np.interp(t, t_u, u_rms))
        if n not in mask_cache:
            mask_cache[n] = _make_mask(n, Ly, n_exp)
        img = _panel(f, tau_pa, mask_cache[n], Ly, out_w, vmax_pa)
        out.append((cyc, tau_mean_pa, u, img))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("warm_dir")
    ap.add_argument("cold_dir")
    ap.add_argument("out_mp4")
    ap.add_argument("--rpm", type=float, required=True,
                    help="target rpm both runs are driven at")
    ap.add_argument("--source-rpm", type=float, required=True,
                    help="rpm the WARM run's initial condition was restarted from")
    ap.add_argument("--theta", type=float, default=7.0)
    ap.add_argument("--width", type=int, default=1000)
    ap.add_argument("--fps", type=float, default=8.0)
    args = ap.parse_args()

    warm_dir, cold_dir = Path(args.warm_dir), Path(args.cold_dir)
    tau_scale = RHO_W * U_bio_of(args.rpm, args.theta) ** 2

    _, _, _, _, _, tau0 = load_frame(sorted((cold_dir / "frames_tau").glob("*.bin"))[-1])
    vmax_pa = float(np.percentile(np.abs(tau0) * tau_scale, 99.5))

    warm = build_series(warm_dir, args.rpm, args.theta, args.width, vmax_pa, tau_scale)
    cold = build_series(cold_dir, args.rpm, args.theta, args.width, vmax_pa, tau_scale)
    n_frames = min(len(warm), len(cold))
    print(f"warm: {len(warm)} frames, cold: {len(cold)} frames -> using {n_frames}")
    print(f"colorbar range: +-{vmax_pa:.3g} Pa")

    panel_h = warm[0][3].size[1]
    total_h = panel_h * 2 + 8
    cbar = _colorbar(total_h, vmax_pa)
    gap = 10

    su = args.source_rpm / args.rpm
    caption_lines = [
        f"target: {args.rpm:g} rpm, theta_max = {args.theta:g} deg   |   "
        f"WARM initial condition: checkpoint restart from a converged "
        f"{args.source_rpm:g} rpm flow (su = {su:.3f})   |   "
        f"COLD initial condition: fluid at rest",
    ]
    cap_font = _font(20)
    tmp = Image.new("RGB", (10, 10))
    cap_h = ImageDraw.Draw(tmp).multiline_textbbox((0, 0), "\n".join(caption_lines),
                                                    font=cap_font)[3] + 16

    canvas_w = warm[0][3].size[0] + gap + cbar.size[0]
    frames = []
    for i in range(n_frames):
        cw, tauw, uw, imgw = warm[i]
        cc, tauc, uc, imgc = cold[i]

        canvas = Image.new("RGB", (canvas_w, cap_h + total_h), (255, 255, 255))
        d = ImageDraw.Draw(canvas)
        d.multiline_text((8, 8), "\n".join(caption_lines), fill=(0, 0, 0), font=cap_font)

        canvas.paste(imgw, (0, cap_h))
        canvas.paste(imgc, (0, cap_h + panel_h + 8))
        canvas.paste(cbar, (imgw.size[0] + gap, cap_h))

        lw = _label_block([f"WARM   cycle {cw:5.1f}", f"tau_mean = {tauw:6.3f} Pa",
                           f"u_rms = {uw:.3f} (nd)"])
        lc = _label_block([f"COLD   cycle {cc:5.1f}", f"tau_mean = {tauc:6.3f} Pa",
                           f"u_rms = {uc:.3f} (nd)"])
        canvas.paste(lw, (8, cap_h + 6), lw)
        canvas.paste(lc, (8, cap_h + panel_h + 14), lc)
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
