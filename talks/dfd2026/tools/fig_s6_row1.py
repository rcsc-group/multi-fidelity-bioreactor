"""Slide 6, round 9, row 1: three whole-tank panels (L8 -> L9 -> L10), equal aspect, with a zoom inset (1/3 of panel width)
at the region of largest |u_L8 - u_L10| (cell edges drawn), thick arrows between panels, level labels in level colours.

Frames: phase 0.25 snapshots as in tools/fig_s6_chain.py. Crop centre: L8 cell (row 125, col 33) (tools/fig_s6_zoom.py).
Window: 32 L8 cells per side (64 L9, 128 L10 cells).
Output: figs/s6_row1.png (10 x 2.5 in = 960 x 240 px on the slide).
"""
import json
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import Normalize

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from scripts.fields import load_snapshot, snapshots, liquid_mask  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s6_row1.png"
RUNS = [("L8", "flow_l8_rpm32.5_th7"), ("L9", "scal_l9_rpm32.5_th7"), ("L10", "v_l10_rpm32.5_th7")]
COL = {"L8": "#0072B2", "L9": "#009E73", "L10": "#CC79A7"}
CI, CJ, WIN = 125, 33, 32
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
vmax = float(np.nanpercentile(spd["L10"], 99))

# bounding box of the liquid (L8 grid), used for all panels (same physical box)
rr = np.where(msk["L8"].any(axis=1))[0]
cc = np.where(msk["L8"].any(axis=0))[0]
r0, r1, c0, c1 = 60, 220, 0, 256   # tank window: rows 60 to 220 (liquid is rows 91 to 147), full width
aspect = (c1 - c0) / (r1 - r0)                  # width / height in cells
print("liquid box rows", r0, r1, "cols", c0, c1, "aspect", round(aspect, 3))

cmap = plt.get_cmap("cividis").copy()
cmap.set_bad("#e6e6e6")                         # domain outside the liquid, light grey (not white)

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14})
W_IN, H_IN = 10.0, 2.5
fig = plt.figure(figsize=(W_IN, H_IN), dpi=200)
PW = 2.95
PH = PW / aspect
GAP = (W_IN - 3 * PW - 0.1) / 2
Y0 = (H_IN - PH) / 2 - 0.05
XS = [0.05 + i * (PW + GAP) for i in range(3)]
print("panel width in", PW, "height in", round(PH, 3), "gap in", round(GAP, 3))
for i, (lab, _) in enumerate(RUNS):
    n_lv = snaps[lab]["n"]
    k = n_lv // 256
    # whole tank, in L8 cell units mapped to the liquid box
    field = spd[lab][r0 * k: r1 * k, c0 * k: c1 * k]
    ax = fig.add_axes([XS[i] / W_IN, Y0 / H_IN, PW / W_IN, PH / H_IN])
    ax.imshow(field, origin="lower", cmap=cmap, vmin=0, vmax=vmax, interpolation="nearest",
              extent=[c0, c1, r0, r1])
    ax.set_xlim(c0, c1)
    ax.set_ylim(r0, r1)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.9)
    # zoom window marker (L8 cell units)
    ax.add_patch(mpatches.Rectangle((CJ - WIN // 2 - 0.5 + 0.5, CI - WIN // 2), WIN, WIN, fill=False, ec="white", lw=1.2, ls="-"))
    # level label above
    fig.text((XS[i] + PW / 2) / W_IN, (Y0 + PH + 0.1) / H_IN, lab, color=COL[lab], fontsize=18, ha="center", va="bottom")
    # inset: 1/3 of panel width, bottom-right corner
    S = PW / 3.0
    ins = fig.add_axes([(XS[i] + PW - S - 0.06) / W_IN, (Y0 + 0.06) / H_IN, S / W_IN, S / H_IN])
    ra, rb = (CI - WIN // 2) * k, (CI + WIN // 2) * k
    ca, cb = (CJ - WIN // 2) * k, (CJ + WIN // 2) * k
    sub = spd[lab][ra:rb, ca:cb]
    ec_alpha, ec_lw = {1: (0.85, 0.7), 2: (0.5, 0.35), 4: (0.3, 0.15)}[k]
    xs = np.arange(ca, cb + 1)
    ys = np.arange(ra, rb + 1)
    ins.pcolormesh(xs, ys, sub, cmap=cmap, vmin=0, vmax=vmax, shading="flat",
                   edgecolors=(0.1, 0.1, 0.1, ec_alpha), linewidth=ec_lw)
    ins.set_xlim(ca, cb)
    ins.set_ylim(ra, rb)
    ins.set_xticks([])
    ins.set_yticks([])
    ins.set_aspect("equal")
    for sp in ins.spines.values():
        sp.set_linewidth(1.2)
        sp.set_edgecolor(COL[lab])
# thick arrows between panels
ymid = (Y0 + PH / 2) / H_IN
for i in range(2):
    xa = (XS[i] + PW + 0.07) / W_IN
    xb = (XS[i + 1] - 0.07) / W_IN
    fig.add_artist(FancyArrowPatch((xa, ymid), (xb, ymid), transform=fig.transFigure, arrowstyle="simple,head_length=1.0,head_width=0.9,tail_width=0.42",
                                   mutation_scale=28, color="0.25", lw=0))
fig.savefig(OUT, dpi=200)
print("wrote", OUT)
