"""Title-slide background v2: REAL-TIME sequence of the L8 flow (32.5 rpm, theta 7; 100 frames per cycle), cycles 81-90 after the oxygen release.
Colour = oxygen concentration (indigo -> cyan -> lime-white), warm glow = log dissipation rate (EDR). Lab frame (rotation by the rocking angle).
No loop in the data: the oxygen keeps entering, so no period repeats. The file fades in from and out to the ground colour (a deliberate loop).
Usage: uv run python render_title_v2.py still [--sign 1] | video
"""
import argparse
import glob
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
from scipy import ndimage

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO)]
from scripts.postprocess import _t_scales  # noqa: E402

RUN = REPO / "runs/rec_l8_rpm32.5_th7"
OUT = Path(__file__).resolve().parent.parent / "public"
GROUND = np.array([0.043, 0.082, 0.140])
EMBER = np.array([1.0, 0.42, 0.12])
OXY = LinearSegmentedColormap.from_list("oxy", ["#0b1626", "#1d2a78", "#1f6fc0", "#14b5cc", "#79e6b4", "#f4f7a8"])
W, H = 1920, 1080


def load(f, n=256, npl=10):
    mm = np.memmap(f, dtype=np.float32, mode="r", offset=16, shape=(npl, n, n))
    return [np.array(mm[k]) for k in (0, 1, 3, 4, 9)]


def render(f, t, T, th_max, sign, state):
    ux, uy, fv, cs, oxy = load(f)
    n = ux.shape[0]; dx = 1.0 / n
    inside = cs > 0.5; liquid = inside & (fv > 0.5); core = ndimage.binary_erosion(liquid, iterations=2)
    q = np.log10(2 * (np.gradient(ux, dx, axis=1) ** 2 + np.gradient(uy, dx, axis=0) ** 2) + (np.gradient(ux, dx, axis=0) + np.gradient(uy, dx, axis=1)) ** 2 + 1e-12)
    lo, hi = state.setdefault("lim", (np.percentile(q[core], 30), np.percentile(q[core], 99.7)))
    w = np.clip((np.where(core, q, lo) - lo) / (hi - lo), 0, 1) ** 2.2
    w = ndimage.gaussian_filter(w, 0.8)
    v = np.clip(ndimage.gaussian_filter(oxy, (1.2, 4.0)) / 0.12, 0, 1) ** 0.5
    col = np.clip(OXY(v)[..., :3] + 0.85 * w[..., None] * EMBER, 0, 1)
    # colour carried into the air by the nearest liquid cell, so the interface has no dark halo; the interface itself comes from the continuous
    # VOF fraction (upsampled and smoothed), not from a thresholded mask: a thresholded mask is what looked pixelated
    idx = ndimage.distance_transform_edt(~liquid, return_distances=False, return_indices=True)
    col = col[idx[0], idx[1]]
    K = 4
    z = lambda a: ndimage.zoom(a, (K, K) + (1,) * (a.ndim - 2), order=3)  # noqa: E731
    col3 = np.clip(z(col), 0, 1)
    f3 = ndimage.gaussian_filter(z(np.clip(fv, 0, 1) * (cs > 0.5)), 0.35 * K)
    s = np.clip((f3 - 0.5) / 0.16 + 0.5, 0, 1); fsm = s * s * (3 - 2 * s)
    s2 = np.clip((ndimage.gaussian_filter(z(cs.astype(float)), 0.3 * K) - 0.5) / 0.12 + 0.5, 0, 1); ins = s2 * s2 * (3 - 2 * s2)
    img = np.concatenate([col3, fsm[..., None], ins[..., None]], axis=2)
    th = sign * th_max * math.sin(2 * math.pi * t / T)
    rot = ndimage.rotate(img, math.degrees(th), reshape=False, order=3, mode="constant", cval=0.0)
    N = rot.shape[0]; h = int(N * 9 / 16)
    crop = rot[N // 2 - h // 2: N // 2 + h // 2, :]
    out = ndimage.zoom(crop, (H / crop.shape[0], W / crop.shape[1], 1), order=3)
    fl = np.clip(out[..., 3:4], 0, 1); al = np.clip(out[..., 4:5], 0, 1)
    rgb = np.clip(out[..., :3], 0, 1) * fl + GROUND * (1 - fl)
    comp = np.clip(rgb * al + GROUND * (1 - al), 0, 1)[::-1]            # array row 0 = lowest y: flip for display
    return (comp * 255).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["still", "video"])
    ap.add_argument("--sign", type=float, default=-1.0)   # -1 puts the free surface level in the lab frame (checked on stills)
    ap.add_argument("--c0", type=float, default=81.0)
    ap.add_argument("--c1", type=float, default=90.0)
    ap.add_argument("--fps", type=int, default=30)
    a = ap.parse_args()
    p = json.load(open(RUN / "params.json")); _, T = _t_scales(p)
    th_max = math.radians(p["theta_max"][0])
    files = sorted(glob.glob(str(RUN / "fields/snap_*.bin")))
    ts = np.array([float(np.fromfile(f, dtype=np.float64, count=1, offset=8)[0]) for f in files])
    sel = [i for i, t in enumerate(ts) if a.c0 <= t / T <= a.c1]
    state = {}
    if a.mode == "still":
        i = min(sel, key=lambda i: abs(ts[i] / T - 89.25))
        Image.fromarray(render(files[i], ts[i], T, th_max, a.sign, state)).save(Path(__file__).resolve().parent / f"still_v2_{'p' if a.sign > 0 else 'm'}.png")
        return
    ff = REPO / ".venv/lib/python3.12/site-packages/imageio_ffmpeg/binaries/ffmpeg-linux-x86_64-v7.0.2"
    nfade = int(1.2 * a.fps); nfr = len(sel)
    cmd = [str(ff), "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(a.fps), "-i", "-",
           "-vf", f"fade=t=in:st=0:d=1.2,fade=t=out:st={nfr / a.fps - 1.2}:d=1.2", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", str(OUT / "title_v2.mp4")]
    pr = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    for k, i in enumerate(sel):
        pr.stdin.write(render(files[i], ts[i], T, th_max, a.sign, state).tobytes())
        if k % 100 == 0:
            print(k, nfr, flush=True)
    pr.stdin.close(); pr.wait(); print("wrote", OUT / "title_v2.mp4")


if __name__ == "__main__":
    main()
