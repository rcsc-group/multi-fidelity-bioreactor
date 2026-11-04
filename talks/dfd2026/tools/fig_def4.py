"""Definitions slide, four panels (user layout 2026-11-02). One real L8 run: runs/rec_l8_rpm32.5_th7 (32.5 rpm, theta 7, 100 snapshots
per cycle, cycles 78-90; release at cycle 80). Lab frame: each frame is rotated by the rocking angle. Real time.
  videos : tau_loop.mp4 / edr_loop.mp4  cycles 78-80 (before the release, periodic flow, looped on the slide)
           oxy_rel.mp4 / tracer_rel.mp4 cycles 80-90 (from the release; played once, holds the last frame)
           tau = |du/dy + dv/dx| (mu constant), sqrt colour scale; eps = 2 u_x^2 + 2 v_y^2 + (u_y + v_x)^2, log10 scale;
           oxygen (plane 9) and tracer c2 (plane 7) on 0..1.
  plots  : def4_kla.png  C*(t) with the one-period fit at C* = 0.10 (scripts/postprocess.py _kla_period_at_threshold)
           def4_mix.png  chi(t) from the bulk-liquid sums of mix8_new_32p5_7 (diary 2026-10-29), first crossing of 0.5
Usage: uv run --project <repo> python fig_def4.py plots | videos
"""
import functools
import glob
import json
import math
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO)]
from scripts import postprocess as pp  # noqa: E402

RUN = REPO / "runs/rec_l8_rpm32.5_th7"
MIXRUN = Path("/oscar/scratch/eaguerov/replay_gate/runs/mix8_new_32p5_7")
OUT = Path(__file__).resolve().parent.parent / "public"
FF = REPO / ".venv/lib/python3.12/site-packages/imageio_ffmpeg/binaries/ffmpeg-linux-x86_64-v7.0.2"
INK, MUTE, CURVE, ACC = "#14202e", "#7a8797", "#1f6fc0", "#e8590c"
cm = LinearSegmentedColormap.from_list
MAPS = dict(tau=cm("tau", ["#f3f0f6", "#efc3dc", "#e0568f", "#9b2c93", "#3a1d6e"]),
            edr=cm("emb", ["#f6f1ea", "#fdd9a8", "#f7924a", "#c2410c", "#4a1d2a"]),
            oxy=cm("oxy", ["#eef3f8", "#b9dcee", "#5fb7d6", "#1f6fc0", "#1d2a78"]),
            tracer=cm("trc", ["#eef3f8", "#b9dcee", "#5fb7d6", "#1f6fc0", "#1d2a78"]))
p = json.load(open(RUN / "params.json")); T_bio, T = pp._t_scales(p)
TH = math.radians(p["theta_max"][0])
files = sorted(glob.glob(str(RUN / "fields/snap_*.bin")))
ts = np.array([float(np.fromfile(f, dtype=np.float64, count=1, offset=8)[0]) for f in files])
K, WOUT, PAD = 4, 960, 140
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.serif": ["cmr10"], "axes.formatter.use_mathtext": True,
                     "font.size": 15, "axes.labelsize": 17, "axes.linewidth": 1.3, "axes.edgecolor": MUTE,
                     "xtick.color": MUTE, "ytick.color": MUTE, "axes.labelcolor": INK, "axes.spines.top": False, "axes.spines.right": False})


def planes(f):
    mm = np.memmap(f, dtype=np.float32, mode="r", offset=16, shape=(10, 256, 256))
    return [np.array(mm[k]) for k in range(10)]


def quantity(kind, pl):
    ux, uy = pl[0], pl[1]; dx = 1.0 / 256
    uxx, uxy = np.gradient(ux, dx, axis=1), np.gradient(ux, dx, axis=0)
    uyx, uyy = np.gradient(uy, dx, axis=1), np.gradient(uy, dx, axis=0)
    if kind == "tau":
        return np.sqrt(np.abs(uxy + uyx))
    if kind == "edr":
        return np.log10(2 * uxx ** 2 + 2 * uyy ** 2 + (uxy + uyx) ** 2 + 1e-12)
    return pl[9] if kind == "oxy" else pl[7]


LIMSEL = None   # snapshot indices of cycles 78-80: the colour limits of tau/edr come from these, for every video (loop and release match)


@functools.cache
def limits(kind):
    global LIMSEL
    if kind in ("oxy", "tracer"):
        return 0.0, 1.0
    LIMSEL = LIMSEL or [i for i, t in enumerate(ts) if 78.0 <= t / T <= 80.0]
    v = []
    for i in LIMSEL[::20]:
        pl = planes(files[i]); liq = (pl[4] > 0.5) & (pl[3] > 0.5)
        v.append(quantity(kind, pl)[ndimage.binary_erosion(liq, iterations=3)])
    v = np.concatenate(v)
    return (0.0, np.percentile(v, 99.0)) if kind == "tau" else (np.percentile(v, 20), np.percentile(v, 99.7))


def colour(kind, pl, lim):
    """RGB on the 256^2 grid. Velocity-gradient fields are taken only 3 cells inside the liquid and carried outward from there
    (the one-sided gradients across the free surface are not shown); scalars use the whole liquid."""
    liq = (pl[4] > 0.5) & (pl[3] > 0.5)
    gas = (pl[4] > 0.5) & (pl[3] <= 0.5)
    src = (liq & (ndimage.distance_transform_edt(~gas) > 8)) if kind in ("tau", "edr") else liq   # 8 cells off the free surface only
    v = np.clip((quantity(kind, pl) - lim[0]) / (lim[1] - lim[0]), 0, 1)
    if kind == "oxy":
        v = np.clip(v / 0.3, 0, 1) ** 0.5
    if kind in ("tau", "edr"):   # near the surface and walls: smooth fill (normalised Gaussian of the core values), no stripes
        m = src.astype(float)
        fill = ndimage.gaussian_filter(v * m, 4.0) / np.maximum(ndimage.gaussian_filter(m, 4.0), 1e-6)
        v = ndimage.gaussian_filter(np.where(src, v, fill), 1.0)
        src = ndimage.binary_dilation(liq, iterations=1) & (ndimage.gaussian_filter(m, 4.0) > 1e-4)
    col = MAPS[kind](v)[..., :3]
    idx = ndimage.distance_transform_edt(~src, return_distances=False, return_indices=True)
    return col[idx[0], idx[1]]


def frame(i, base, box, scal=None, wb=0.0):
    """base field, cross-faded into the scalar with weight wb; for oxygen the gas turns to the c = 1 colour (the source)."""
    pl = planes(files[i]); fv, cs = pl[3], pl[4]
    col = colour(base, pl, limits(base))
    if scal and wb > 0:
        col = (1 - wb) * col + wb * colour(scal, pl, limits(scal))
    z = lambda a: ndimage.zoom(a, (K, K) + (1,) * (a.ndim - 2), order=3)  # noqa: E731
    f3 = ndimage.gaussian_filter(z(np.clip(fv, 0, 1) * (cs > 0.5)), 0.35 * K)
    s = np.clip((f3 - 0.5) / 0.16 + 0.5, 0, 1); al = s * s * (3 - 2 * s)
    s2 = np.clip((ndimage.gaussian_filter(z(cs.astype(float)), 0.3 * K) - 0.5) / 0.12 + 0.5, 0, 1); wall = s2 * s2 * (3 - 2 * s2)
    img = np.concatenate([np.clip(z(col), 0, 1), al[..., None], wall[..., None]], axis=2)
    th = -TH * math.sin(2 * math.pi * ts[i] / T)
    img = np.pad(img, ((PAD, PAD), (PAD, PAD), (0, 0)))
    img = ndimage.rotate(img, math.degrees(th), reshape=False, order=1, mode="constant", cval=0.0)
    r0, r1, c0, c1 = box
    img = img[r0:r1, c0:c1]
    a = np.clip(img[..., 3:4], 0, 1); w = np.clip(img[..., 4:5], 0, 1)
    ground = np.ones(3); grey = np.array([0.93, 0.94, 0.96])
    g = wb if scal == "oxy" else 0.0
    gas = (1 - g) * grey + g * np.array(MAPS["oxy"](1.0)[:3])
    rgb = img[..., :3] * a + gas * (1 - a)
    rgb = rgb * w + ground * (1 - w)
    edge = np.clip(ndimage.gaussian_gradient_magnitude(img[..., 4], 1.0) * 3.0, 0, 1)[..., None]   # bag outline
    rgb = rgb * (1 - edge) + np.array([0.48, 0.53, 0.59]) * edge
    rgb = rgb[::-1]
    hh = int(round(rgb.shape[0] * WOUT / rgb.shape[1] / 2)) * 2
    out = ndimage.zoom(rgb, (hh / rgb.shape[0], WOUT / rgb.shape[1], 1), order=1)
    return (np.clip(out, 0, 1) * 255).astype(np.uint8)


def crop_box():
    """Fixed crop around the bag that holds it at every rocking angle (K-upsampled pixels)."""
    cs = np.pad(ndimage.zoom((planes(files[0])[4] > 0.5).astype(float), K, order=1), PAD)
    u = np.zeros_like(cs, dtype=bool)
    for a in np.linspace(-TH, TH, 9):          # union of the bag over the rocking angles, rotated exactly as frame() does
        u |= ndimage.rotate(cs, math.degrees(a), reshape=False, order=1) > 0.05
    r = np.where(u.any(axis=1))[0]; c = np.where(u.any(axis=0))[0]
    # gas crop: top of the box = highest liquid row over all frames (each rotated by its own angle), plus a margin
    top = 0
    for i in range(0, len(files), 10):
        pl = planes(files[i])
        liq = np.pad(ndimage.zoom(((pl[4] > 0.5) & (pl[3] > 0.5)).astype(float), K, order=1), PAD)
        th = -TH * math.sin(2 * math.pi * ts[i] / T)
        top = max(top, int(np.where((ndimage.rotate(liq, math.degrees(th), reshape=False, order=1) > 0.05).any(axis=1))[0][-1]))
    n = cs.shape[0]; m = 6
    return max(r[0] - m, 0), min(top + 28, r[-1] + m + 1, n), max(c[0] - m, 0), min(c[-1] + m + 1, n)



FADE = 1.0   # cycles of cross-fade after the release


def video(name, base, scal=None):
    """Loop (scal None): cycle 79-80, 100 frames, periodic. Release: the SAME cycle 79-80 frames, then the release (cycle 80) with a
    cross-fade into the scalar, to cycle 90. The slide starts the release video at the loop's current time: no jump."""
    sel = [i for i, t in enumerate(ts) if 79.0 <= t / T < 80.0] if scal is None else [i for i, t in enumerate(ts) if 79.0 <= t / T <= 90.0]
    box = crop_box()
    fps = 100.0 / (T * T_bio)              # 100 snapshots per cycle: real time
    wb = lambda i: float(np.clip((ts[i] / T - 80.0) / FADE, 0, 1)) if scal else 0.0  # noqa: E731
    first = frame(sel[0], base, box)
    h, w = first.shape[:2]
    cmd = [str(FF), "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{w}x{h}", "-r", f"{fps:.4f}", "-i", "-",
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-g", "25", "-movflags", "+faststart", str(OUT / name)]
    pr = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    for i in sel:
        pr.stdin.write(frame(i, base, box, scal, wb(i)).tobytes())
    pr.stdin.close(); pr.wait()
    print("wrote", name, len(sel), "frames", f"{fps:.2f} fps", f"{w}x{h}")
def axes(figsize):
    fig, ax = plt.subplots(figsize=figsize)
    ax.patch.set_alpha(0)
    return fig, ax


def plot_kla():
    t, c = pp._compute_c_star(RUN)
    s = (t - 80.0 * T) * T_bio
    k = int(np.argmax(c >= 0.10)); half = 0.5 * T
    w = (t >= t[k] - half) & (t <= t[k] + half)
    slope, _ = np.polyfit(t[w], np.log(1 - np.clip(c[w], 0, 1 - 1e-10)), 1)
    kla = -slope * 3600 / T_bio
    m = (s >= 0) & (s <= 30)
    fig, ax = axes((4.6, 2.45))
    ax.plot(s[m], c[m], color=CURVE, lw=2.8)
    sl = -slope * (1 - c[k]) / T_bio
    xs = np.array([s[k] - 9, s[k] + 9]); ax.plot(xs, c[k] + sl * (xs - s[k]), color=ACC, lw=2.0, ls=(0, (4, 3)))
    ax.plot(s[k], c[k], "o", mfc="none", mec=ACC, mew=2.4, ms=12, zorder=5)
    ax.set_xlim(0, 30); ax.set_ylim(0, 0.3)
    ax.set_xlabel(r"$t$ (s)"); ax.set_ylabel(r"$C^*$")
    fig.tight_layout(pad=0.3); fig.savefig(OUT / "def4_kla.png", dpi=220, transparent=True); plt.close(fig)
    print(f"kLa {kla:.1f} 1/h at C* = 0.10, t = {s[k]:.2f} s")


def plot_mix():
    tn, b0, b1, b2 = pp._joined_cols(MIXRUN, "tr_oxy.dat", [pp._COL_B0, pp._COL_B1, pp._COL_B2])
    pn = json.load(open(MIXRUN / "params.json")); Tb_n, _ = pp._t_scales(pn)
    r = np.loadtxt(MIXRUN / "mix_release.dat", skiprows=1)
    s2max = r[4] / r[2] - (r[3] / r[2]) ** 2
    k = tn > r[1]
    s = np.concatenate([[0.0], (tn[k] - r[1]) * Tb_n])
    chi = np.concatenate([[0.0], np.clip(1 - (b2[k] / b0[k] - (b1[k] / b0[k]) ** 2) / s2max, 0, 1)])
    m = s <= 30
    tm = s[int(np.argmax(chi >= 0.5))]
    fig, ax = axes((4.6, 2.45))
    ax.plot(s[m], chi[m], color=CURVE, lw=2.8)
    ax.plot([0, tm], [0.5, 0.5], color=MUTE, lw=1.1, ls=(0, (3, 3))); ax.plot([tm, tm], [0, 0.5], color=ACC, lw=2.4)
    ax.plot(tm, 0.5, "o", color=ACC, ms=8, zorder=5)
    ax.set_xlim(0, 30); ax.set_ylim(0, 1.0)
    ax.set_xlabel(r"$t$ (s)"); ax.set_ylabel(r"$\chi$")
    fig.tight_layout(pad=0.3); fig.savefig(OUT / "def4_mix.png", dpi=220, transparent=True); plt.close(fig)
    print(f"t_mix(0.5) {tm:.2f} s")


if __name__ == "__main__":
    if sys.argv[1] == "plots":
        plot_kla(); plot_mix()
    elif sys.argv[1] == "still":   # still BASE CYCLE [SCAL] -> /oscar/scratch/eaguerov/tmp/still_BASE.png (render check)
        from PIL import Image
        base, cyc = sys.argv[2], float(sys.argv[3]); scal = sys.argv[4] if len(sys.argv) > 4 else None
        i = int(np.argmin(np.abs(ts / T - cyc)))
        wb = float(np.clip((cyc - 80.0) / FADE, 0, 1)) if scal else 0.0
        Image.fromarray(frame(i, base, crop_box(), scal, wb)).save(f"/oscar/scratch/eaguerov/tmp/still_{base}.png")
    else:   # videos [names]: one process per video runs them in parallel
        spec = {"tau_loop.mp4": ("tau", None), "oxy_rel.mp4": ("tau", "oxy"), "edr_loop.mp4": ("edr", None), "tracer_rel.mp4": ("edr", "tracer")}
        for name in sys.argv[2:] or spec:
            video(name, *spec[name])
