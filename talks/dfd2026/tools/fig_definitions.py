"""Definitions slide: four TEXT-FREE panels built from the code's own definitions and one real L8 run (32.5 rpm, theta 7; runs/rec_l8_rpm32.5_th7).
 kLa   : C*(t) after the oxygen release; ln(1 - C*) fitted over ONE rocking period centred on C* = 0.10 (scripts/postprocess.py _kla_period_at_threshold)
 mixing: chi(t) = 1 - sigma^2/sigma^2_max of the tracer, first crossing of 0.5; snapshots of the tracer c2 at 0, t_mix/2 and t_mix
 shear : tau = mu (du/dy + dv/dx) on the liquid, tau_95 = 95th percentile over the liquid (src/BioReactor.c, event normcal)
 EDR   : eps = mu (2 u_x^2 + 2 v_y^2 + (u_y + v_x)^2) on the liquid
Symbols and words are set by the slide (deck font), not in the images. Transparent PNGs ../public/def_<name>_<mode>.png.
Usage: uv run python fig_definitions.py light|dark
"""
import glob
import json
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

MODE = sys.argv[1] if len(sys.argv) > 1 else "light"
RUN = REPO / "runs/rec_l8_rpm32.5_th7"
OUT = Path(__file__).resolve().parent.parent / "public"
cm = LinearSegmentedColormap.from_list
if MODE == "dark":
    MUTE, CURVE, ACC, HIST = "#8fa3bf", "#14b5cc", "#ff8a1f", "#9b2c93"
    OXY = cm("oxy", ["#0b1626", "#1d2a78", "#1f6fc0", "#14b5cc", "#79e6b4", "#f4f7a8"])
    TAU = cm("tau", ["#0b1626", "#3a1d6e", "#9b2c93", "#ff5fa2", "#ffe3ee"])
    EMB = cm("emb", ["#0b1626", "#4a1d2a", "#b8441f", "#ff8a1f", "#fff0c2"])
else:
    MUTE, CURVE, ACC, HIST = "#7a8797", "#1f6fc0", "#e8590c", "#b03a8c"
    OXY = cm("oxy", ["#e9f0f7", "#b9dcee", "#5fb7d6", "#1f6fc0", "#1d2a78"])
    TAU = cm("tau", ["#f3f0f6", "#efc3dc", "#e0568f", "#9b2c93", "#3a1d6e"])
    EMB = cm("emb", ["#f6f1ea", "#fdd9a8", "#f7924a", "#c2410c", "#4a1d2a"])
p = json.load(open(RUN / "params.json")); T_bio, T = pp._t_scales(p)
files = sorted(glob.glob(str(RUN / "fields/snap_*.bin")))
ts = np.array([float(np.fromfile(f, dtype=np.float64, count=1, offset=8)[0]) for f in files])
t_rel = 80.0 * T    # tracer and oxygen released at cycle 80 (n_mix_cycles)
K = 4               # upsampling of the 256^2 snapshots for display


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTE); ax.spines[s].set_linewidth(1.4)
    ax.set_xticks([]); ax.set_yticks([])
    ax.patch.set_alpha(0)


def save(fig, name):
    fig.savefig(OUT / f"def_{name}_{MODE}.png", dpi=200, transparent=True)
    plt.close(fig)
    print("wrote", name, MODE)


def planes(i):
    mm = np.memmap(files[i], dtype=np.float32, mode="r", offset=16, shape=(10, 256, 256))
    return [np.array(mm[k]) for k in range(10)]


def field_img(q, mask, cmap, lo, hi):
    """RGBA of q on the liquid: cubic upsampling, soft edge from the upsampled mask (no staircase)."""
    qf = np.where(mask, q, np.nanmean(q[mask]))
    v = np.clip((ndimage.zoom(qf, K, order=3) - lo) / (hi - lo), 0, 1)
    a = ndimage.gaussian_filter(ndimage.zoom(mask.astype(float), K, order=1), 0.35 * K)
    s = np.clip((a - 0.5) / 0.16 + 0.5, 0, 1)
    rgba = cmap(v)
    rgba[..., 3] = s * s * (3 - 2 * s)
    return rgba


def bbox(mask):
    r = np.where(mask.any(axis=1))[0]; c = np.where(mask.any(axis=0))[0]
    return r[0], r[-1] + 1, c[0], c[-1] + 1


def crop(img, bb):
    r0, r1, c0, c1 = bb
    return img[r0 * K:r1 * K, c0 * K:c1 * K]


def fig_kla():
    t, c = pp._compute_c_star(RUN)
    s = (t - t_rel) * T_bio
    m = (s >= 0) & (s <= 120)
    k = int(np.argmax(c >= 0.10)); half = 0.5 * T
    w = (t >= t[k] - half) & (t <= t[k] + half)
    slope, ic = np.polyfit(t[w], np.log(1 - np.clip(c[w], 0, 1 - 1e-10)), 1)
    fig, ax = plt.subplots(figsize=(6.0, 3.0))
    ax.plot(s[m], c[m], color=CURVE, lw=3.0)
    # slope of the fitted law at the crossing, drawn as a tangent: dC*/dt = kLa (1 - C*) with C* = 0.10
    sl = -slope * (1 - c[k]) / T_bio
    xs = np.array([0.0, 40.0]); ax.plot(xs, c[k] + sl * (xs - s[k]), color=ACC, lw=2.2, ls=(0, (4, 3)))
    ax.plot(s[k], c[k], "o", mfc="none", mec=ACC, mew=2.6, ms=15, zorder=5)
    ax.set_ylim(0, 1); ax.set_xlim(0, 120)
    style(ax); fig.tight_layout(); save(fig, "kla")
    print("kLa fit", -slope * 3600 / T_bio, "1/h at C* = 0.10")


def fig_mix():
    t, c1, c2 = pp._joined_cols(RUN, "tr_oxy.dat", [8, 9])
    _, fl = pp._joined_cols(RUN, "vol_frac_interf.dat", [pp._COL_F_LIQ])
    fm = fl.mean(); mean = c1 / fm; var = c2 / fm - mean ** 2
    i0 = int(np.nonzero(c1 > 1e-10 * fm)[0][0]); chi = np.clip(1 - var / var[i0], 0, 1); chi[:i0] = 0
    s = (t - t[i0]) * T_bio
    m = (s >= 0) & (s <= 40)
    tm = s[i0 + int(np.argmax(chi[i0:] >= 0.5))]
    fig = plt.figure(figsize=(6.0, 2.9))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.35, 1], wspace=0.12, hspace=0.18)
    ax = fig.add_subplot(gs[:, 0])
    ax.plot(s[m], chi[m], color=CURVE, lw=3.0)
    ax.plot([0, tm], [0.5, 0.5], color=MUTE, lw=1.2, ls=(0, (3, 3))); ax.plot([tm, tm], [0, 0.5], color=ACC, lw=2.6)
    ax.plot(tm, 0.5, "o", color=ACC, ms=9, zorder=5)
    ax.set_ylim(0, 1.02); ax.set_xlim(0, 40)
    style(ax)
    for j, sec in enumerate((0.0, 0.5 * tm, tm)):
        a2 = fig.add_subplot(gs[j, 1])
        pl = planes(int(np.argmin(np.abs((ts - t[i0]) * T_bio - sec))))
        f, cs, cc2 = pl[3], pl[4], pl[7]
        mask = (cs > 0.5) & (f > 0.5)
        a2.imshow(crop(field_img(cc2, mask, OXY, 0.0, 1.0), bbox(mask)), origin="lower", interpolation="lanczos"); a2.set_axis_off()
    fig.subplots_adjust(left=0.05, right=0.98, top=0.98, bottom=0.04)
    save(fig, "mix")
    print("t_mix 50 %", tm, "s")


def grads(ux, uy, n=256):
    dx = 1.0 / n
    return np.gradient(ux, dx, axis=1), np.gradient(ux, dx, axis=0), np.gradient(uy, dx, axis=1), np.gradient(uy, dx, axis=0)


def fig_field(kind):
    i = int(np.argmin(np.abs(ts / T - 86.30)))
    pl = planes(i); ux, uy, f, cs = pl[0], pl[1], pl[3], pl[4]
    uxx, uxy, uyx, uyy = grads(ux, uy)
    mask = (cs > 0.5) & (f > 0.5)
    bb = bbox(mask)
    fig, ax = plt.subplots(figsize=(6.0, 2.0))
    if kind == "tau":
        q = np.abs(uxy + uyx)
        img = field_img(np.sqrt(q), mask, TAU, 0.0, np.sqrt(np.percentile(q[mask], 99.0)))
        p95 = np.percentile(q[mask], 95)
    else:
        q = np.log10(2 * uxx ** 2 + 2 * uyy ** 2 + (uxy + uyx) ** 2 + 1e-12)
        img = field_img(q, mask, EMB, np.percentile(q[mask], 20), np.percentile(q[mask], 99.7))
    ax.imshow(crop(img, bb), origin="lower", interpolation="lanczos"); ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1); save(fig, kind)
    if kind == "tau":   # small histogram with the 95th percentile: inset on the slide
        fig, ax = plt.subplots(figsize=(3.0, 1.7))
        v = q[mask]
        ax.hist(v, bins=80, color=HIST, range=(0, np.percentile(v, 99.5)), histtype="stepfilled")
        ax.axvline(p95, color=ACC, lw=2.6)
        style(ax); ax.spines["left"].set_visible(False); fig.tight_layout(); save(fig, "tau_hist")


if __name__ == "__main__":
    fig_kla(); fig_mix(); fig_field("tau"); fig_field("edr")
