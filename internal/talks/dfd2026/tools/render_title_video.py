"""Title-slide background: ONE rocking period of the L10 flow (32.5 rpm, theta 7), lab frame, |u| in the water, seamless loop.
Source: runs/v_l10_rpm32.5_th7/fields/snap_*.bin (float32 planes ux, uy, omega, f, cs, ...; header int n, int nplanes, double t).
Frames inside one period are resampled to M evenly spaced phases (linear in time between neighbouring snapshots); the flow of the settled
cycles repeats to ~1 %, so frame M equals frame 0 and the loop shows no seam.
Usage: uv run python render_title_video.py still|video  [--cycle 250] [--m 96] [--wfrac 0.5]
"""
import argparse
import glob
import math
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO)]
from scripts.postprocess import _t_scales  # noqa: E402
import json  # noqa: E402

RUN = REPO / "runs/v_l10_rpm32.5_th7"
OUT = Path(__file__).resolve().parent.parent / "public"


def header(f):
    h = np.fromfile(f, dtype=np.int32, count=2)
    t = float(np.fromfile(f, dtype=np.float64, count=1, offset=8)[0])
    return int(h[0]), int(h[1]), t


def planes(f, n, npl, which=(0, 1, 3, 4)):
    mm = np.memmap(f, dtype=np.float32, mode="r", offset=16, shape=(npl, n, n))
    return [np.array(mm[k]) for k in which]


GROUND = np.array([0.043, 0.082, 0.140])   # deep navy
AIRC = GROUND   # no bag outline and no air tint: only the water is drawn


def field_rgb(kind, ux, uy, om, f, cs, n, state):
    """RGB (n, n, 3) of the chosen field in the water on a dark ground; state carries the fixed colour limits from the first frame."""
    inside = cs > 0.5; liquid = inside & (f > 0.5)
    core = ndimage.binary_erosion(liquid, iterations=3)
    if kind == "speed":
        q = np.hypot(ux, uy); lo, hi = 0.0, state.setdefault("hi", np.percentile(q[core], 99)); cm = plt.get_cmap("viridis"); v = (q - lo) / (hi - lo)
    elif kind == "edr":
        dx = 1.0 / n
        uxx = np.gradient(ux, dx, axis=1); uxy = np.gradient(ux, dx, axis=0); uyx = np.gradient(uy, dx, axis=1); uyy = np.gradient(uy, dx, axis=0)
        q = np.log10(2 * (uxx ** 2 + uyy ** 2) + (uxy + uyx) ** 2 + 1e-12)
        lo, hi = state.setdefault("lim", (np.percentile(q[core], 8), np.percentile(q[core], 99.7)))
        cm = plt.get_cmap("magma"); v = 0.24 + 0.76 * np.clip((q - lo) / (hi - lo), 0, 1)
    else:
        q = om; a = state.setdefault("a", np.percentile(np.abs(q[core]), 97)); cm = plt.get_cmap("RdBu_r"); v = 0.5 + 0.5 * np.clip(q / a, -1, 1)
    q_s = ndimage.gaussian_filter(np.where(core, v, np.nan_to_num(np.nanmean(v[core]))), 0.9)
    col = cm(np.clip(q_s, 0, 1))[..., :3]
    fsm = ndimage.gaussian_filter(liquid.astype(float), 1.3)
    ins = ndimage.gaussian_filter(inside.astype(float), 0.8)
    rgb = np.zeros((n, n, 3))
    for c3 in range(3):
        rgb[..., c3] = fsm * col[..., c3] + (1 - fsm) * AIRC[c3]
    return rgb, ins


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["still", "video"])
    ap.add_argument("--cycle", type=float, default=250.0)
    ap.add_argument("--m", type=int, default=96)
    ap.add_argument("--wfrac", type=float, default=0.5)
    ap.add_argument("--field", choices=["speed", "edr", "vort"], default="edr")
    ap.add_argument("--xc", type=float, default=-0.18, help="window centre, domain units")
    a = ap.parse_args()
    p = json.load(open(RUN / "params.json")); _, T = _t_scales(p)
    files = sorted(glob.glob(str(RUN / "fields/snap_*.bin")))
    meta = [header(f) for f in files]
    ts = np.array([m[2] for m in meta]); n, npl = meta[0][0], meta[0][1]
    th_max = math.radians(p["theta_max"][0] if isinstance(p["theta_max"], list) else p["theta_max"])
    t0 = a.cycle * T
    M = a.m
    cyc = ts / T; ok = cyc >= cyc.min() + 3          # skip the first 3 cycles of the window (settling)
    ph = np.mod(cyc, 1.0)
    # phase pooling: frame spacing dt/T = 0.0735 is incommensurate with the period, so frames of ~17 cycles tile one period densely
    # (max gap ~1 frame of 96); each target phase uses the nearest real frame, no blending in time, so the interface never smears or jumps
    targets = [(k / M + 0.0) for k in range(M)] if a.mode == "video" else [0.25]
    picks = []
    for tp in targets:
        d = np.abs(((ph - tp + 0.5) % 1.0) - 0.5); d[~ok] = 9
        picks.append(int(np.argmin(d)))
    phases = ts[picks]
    state = {}
    W, H = 1920, 1080
    tmp = Path(__file__).resolve().parent / "frames"; tmp.mkdir(exist_ok=True)
    for k, t in enumerate(phases):
        A = planes(files[picks[k]], n, npl, (0, 1, 2, 3, 4))
        ux, uy, om, f, cs = A
        rgb, ins = field_rgb(a.field, ux, uy, om, f, cs, n, state)
        img = np.concatenate([rgb, ins[..., None]], axis=2)
        th = th_max * math.sin(2 * math.pi * t / T)
        rot = ndimage.rotate(img, math.degrees(th), reshape=False, order=1, mode="constant", cval=0.0)
        rot[..., 3] = np.clip(rot[..., 3], 0, 1)
        w = int(a.wfrac * n); h = int(w * H / W)
        cx = int(n * (0.5 + a.xc)); cy = n // 2
        crop = rot[cy - h // 2: cy + h // 2, cx - w // 2: cx + w // 2]
        fig = plt.figure(figsize=(W / 200, H / 200), dpi=200, facecolor=GROUND)
        ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
        ax.imshow(crop, origin="lower", interpolation="bicubic", aspect="auto")
        fn = tmp / (f"still.png" if a.mode == "still" else f"f_{k:04d}.png")
        fig.savefig(fn, facecolor=GROUND); plt.close(fig)
    if a.mode == "video":
        ff = REPO / ".venv/lib/python3.12/site-packages/imageio_ffmpeg/binaries/ffmpeg-linux-x86_64-v7.0.2"
        subprocess.run([str(ff), "-y", "-framerate", "30", "-i", str(tmp / "f_%04d.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "23",
                        str(OUT / "title_loop.mp4")], check=True)
        print("wrote", OUT / "title_loop.mp4")


if __name__ == "__main__":
    main()
