"""Slide 6 right side (round 7): one line plot, the same speed profile from L8, L9 and L10 overlaid in the level colours.

Frames: phase-matched snapshots as in tools/fig_s6_chain.py (phase 0.25):
  runs/flow_l8_rpm32.5_th7/fields (2026-10-11 03:18), runs/scal_l9_rpm32.5_th7/fields (2026-10-11 08:18),
  runs/v_l10_rpm32.5_th7/fields (2026-10-12 14:52).
Location: the L8 cell (row 125, col 33), where the block-averaged |u_L8 - u_L10| is largest (tools/fig_s6_zoom.py).
Two candidate profiles through that location:
  V: speed along the vertical line x = col 33, from the bottom wall to row 125 (liquid cells only)
  H: speed along the horizontal line y = row 125, across the domain (liquid cells only)
Separation = RMS(L8 - L10) / RMS(L10) over the liquid points, with L8 and L10 on the same 256-cell grid (L10 block-averaged
4x4, liquid cells only). The profile with the larger separation is plotted (rule fixed before the numbers were read).
Output: figs/s6_profile.png (4.9 x 4.9 in), data/s6_profile.json
"""
import json
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from scripts.fields import load_snapshot, snapshots, liquid_mask  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s6_profile.png"
RUNS = [("L8", "flow_l8_rpm32.5_th7"), ("L9", "scal_l9_rpm32.5_th7"), ("L10", "v_l10_rpm32.5_th7")]
COL = {"L8": "#0072B2", "L9": "#009E73", "L10": "#CC79A7"}
COST = {"L8": "50 core-h", "L9": "550 core-h", "L10": "1,139 core-h"}
ROW0, COL0 = 125, 33              # L8 cell of the crop centre (from fig_s6_zoom.py)
TARGET_PHASE = 0.25

picked = {}
for lab, run in RUNS:
    d = REPO / "runs" / run
    _, T = _t_scales(json.load(open(d / "params.json")))
    files = snapshots(d)
    tt = []
    for f in files:
        with open(f, "rb") as fh:
            fh.read(8)
            (t,) = struct.unpack("<d", fh.read(8))
        tt.append(t)
    ph = (np.array(tt) / T) % 1.0
    dph = np.minimum(np.abs(ph - TARGET_PHASE), 1 - np.abs(ph - TARGET_PHASE))
    picked[lab] = files[int(np.argmin(dph))]

snaps = {lab: load_snapshot(picked[lab]) for lab, _ in RUNS}
spd, msk = {}, {}
for lab in snaps:
    s = snaps[lab]
    msk[lab] = liquid_mask(s)
    spd[lab] = np.where(msk[lab], np.hypot(s["ux"], s["uy"]), np.nan)


def block_mean(a, k):
    n = a.shape[0] // k
    b = a.reshape(n, k, n, k)
    with np.errstate(all="ignore"):
        return np.nanmean(b, axis=(1, 3))


# common 256 grid
s8 = spd["L8"]
s9_256 = block_mean(spd["L9"], 2)
s10_256 = block_mean(spd["L10"], 4)
n256 = 256
print("shapes", s8.shape, s9_256.shape, s10_256.shape)


def profile(arr, n, kind):
    k = arr.shape[0] // n
    if kind == "V":
        col = COL0 * k
        return arr[:, col]
    row = ROW0 * k
    return arr[row, :]


cand = {}
for kind in ("V", "H"):
    x8 = np.arange(n256) if kind == "H" else np.arange(n256)
    p8 = s8[:, COL0] if kind == "V" else s8[ROW0, :]
    p10 = s10_256[:, COL0] if kind == "V" else s10_256[ROW0, :]
    p9 = s9_256[:, COL0] if kind == "V" else s9_256[ROW0, :]
    ok = np.isfinite(p8) & np.isfinite(p10)
    sep10 = np.sqrt(np.mean((p8[ok] - p10[ok]) ** 2)) / np.sqrt(np.mean(p10[ok] ** 2))
    ok9 = np.isfinite(p8) & np.isfinite(p9)
    sep9 = np.sqrt(np.mean((p8[ok9] - p9[ok9]) ** 2)) / np.sqrt(np.mean(p9[ok9] ** 2))
    cand[kind] = dict(n_points=int(ok.sum()), sep_L8_vs_L10=float(sep10), sep_L8_vs_L9=float(sep9))
    print(kind, cand[kind])
choice = max(("V", "H"), key=lambda k: cand[k]["sep_L8_vs_L10"])
print("choice:", choice)

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14, "axes.linewidth": 1.1,
                     "xtick.direction": "in", "ytick.direction": "in"})
fig, ax = plt.subplots(figsize=(4.9, 4.9), dpi=200)
ends = {}
for lab, _ in RUNS:
    n = spd[lab].shape[0]
    arr = spd[lab]
    if choice == "V":
        col = COL0 * (n // n256)
        y = arr[: ROW0 * (n // n256) + 1, col]
        x = (np.arange(len(y)) + 0.5) / n - 0.0
        xx = x
        vals = y
    else:
        row = ROW0 * (n // n256)
        vals = arr[row, :]
        xx = (np.arange(n) + 0.5) / n
    ok = np.isfinite(vals)
    xx, vals = xx[ok], vals[ok]
    ax.plot(xx, vals, color=COL[lab], ls=":", lw=1.6)
    ends[lab] = (xx[-1], vals[-1])
    ax.text(xx[-1] + 0.01, vals[-1], f"{lab}  {COST[lab]}", color=COL[lab], fontsize=14, va="center", ha="left")
ax.set_xlim(0.352, 0.672)
ax.set_xlabel("distance (domain units)" if choice == "V" else "position across the domain (domain units)", fontsize=14)
ax.set_ylabel(r"speed $|u|$ (nondimensional)", fontsize=14)
ax.tick_params(labelsize=13)
ax.grid(which="major", ls=":", alpha=0.4)
fig.subplots_adjust(left=0.16, right=0.98, top=0.97, bottom=0.14)
fig.savefig(OUT, dpi=200)
(TALK / "data" / "s6_profile.json").write_text(json.dumps(dict(candidates=cand, choice=choice, cell=[ROW0, COL0],
                                                               picked={k: picked[k].name for k in picked}), indent=1))
print("wrote", OUT)
