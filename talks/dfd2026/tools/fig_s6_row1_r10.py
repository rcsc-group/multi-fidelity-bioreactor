"""Slide 6, round 10, row 1: three tanks (L8 -> L9 -> L10), cropped to the liquid plus a thin gas band, equal aspect,
ONE shared colour scale, one shared colour bar with tick labels at the right, thin arrows with 'warm start' above.
No zoom insets, no white boxes (user, round 10).

Frames: phase 0.25 snapshots as in tools/fig_s6_chain.py.
Crop: L8 grid rows 84..156 (liquid rows 91..147 plus a thin gas band), all 256 columns.
Colour limits: vmin, vmax = 5th and 60th percentile of |u| over the liquid cells of the three levels combined.
  Each level contributes its liquid cells on the common 256-cell grid (block mean of liquid cells, cells with >= half
  liquid), so each level has equal weight in the pool.
Saturation: fraction of each level's liquid cells (full resolution) with |u| > vmax (saturated high) and < vmin (low).
Output: figs/s6_row1.png (10 x 1.45 in = 960 x 139 px), data/s6_row1_r10.json
"""
import json
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import Normalize
import matplotlib.cm as cm

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from scripts.fields import load_snapshot, snapshots, liquid_mask  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s6_row1.png"
RUNS = [("L8", "flow_l8_rpm32.5_th7"), ("L9", "scal_l9_rpm32.5_th7"), ("L10", "v_l10_rpm32.5_th7")]
COL = {"L8": "#0072B2", "L9": "#009E73", "L10": "#CC79A7"}
R0, R1 = 84, 156
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


def block_liquid(a, m, k):
    n = a.shape[0] // k
    A = (a * m).reshape(n, k, n, k).sum(axis=(1, 3))
    M = m.reshape(n, k, n, k).sum(axis=(1, 3))
    out = np.where(M >= 0.5 * k * k, A / np.maximum(M, 1), np.nan)
    return out


pool = []
for lab in snaps:
    k = snaps[lab]["n"] // 256
    if k == 1:
        vals = spd[lab]
        grid = block_liquid(np.nan_to_num(spd[lab]), msk[lab], 1)
    else:
        grid = block_liquid(np.nan_to_num(spd[lab]), msk[lab], k)
    pool.append(grid[np.isfinite(grid)])
pool = np.concatenate(pool)
vmin, vmax = (float(np.percentile(pool, 5)), float(np.percentile(pool, 60)))
print(f"colour limits (5th to 60th percentile, pooled, each level equal weight): vmin={vmin:.4f} vmax={vmax:.4f}")
sat = {}
for lab in snaps:
    v = spd[lab][msk[lab]]
    sat[lab] = dict(high=float(np.mean(v > vmax)), low=float(np.mean(v < vmin)), n_liquid=int(v.size))
    print(lab, "liquid cells", v.size, "above vmax", round(100 * sat[lab]["high"], 2), "% below vmin", round(100 * sat[lab]["low"], 2), "%")
sat_total_high = float(sum(sat[l]["high"] * sat[l]["n_liquid"] for l in sat) / sum(sat[l]["n_liquid"] for l in sat))
sat_total_low = float(sum(sat[l]["low"] * sat[l]["n_liquid"] for l in sat) / sum(sat[l]["n_liquid"] for l in sat))
print(f"all levels (cell-weighted): above vmax {100*sat_total_high:.2f} %, below vmin {100*sat_total_low:.2f} %")

cmap = plt.get_cmap("cividis").copy()
cmap.set_bad("#e6e6e6")
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14})
W_IN = 10.0
PW = 2.6
aspect = 256 / (R1 - R0)
PH = PW / aspect
H_IN = PH + 0.42
GAP = 0.4
CB_W = 0.16
X0 = 0.05
XS = [X0 + i * (PW + GAP) for i in range(3)]
print("panel width", PW, "height", round(PH, 3), "figure height", round(H_IN, 3))
fig = plt.figure(figsize=(W_IN, H_IN), dpi=200)
Y0 = 0.02
for i, (lab, _) in enumerate(RUNS):
    k = snaps[lab]["n"] // 256
    field = spd[lab][R0 * k: R1 * k, :]
    ax = fig.add_axes([XS[i] / W_IN, Y0 / H_IN, PW / W_IN, PH / H_IN])
    ax.imshow(field, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest",
              extent=[0, 256, R0, R1])
    ax.set_xlim(0, 256)
    ax.set_ylim(R0, R1)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.9)
    fig.text((XS[i] + PW / 2) / W_IN, (Y0 + PH + 0.06) / H_IN, lab, color=COL[lab], fontsize=18, ha="center", va="bottom")
ymid = (Y0 + PH / 2) / H_IN
for i in range(2):
    xa = (XS[i] + PW + 0.06) / W_IN
    xb = (XS[i + 1] - 0.06) / W_IN
    fig.add_artist(FancyArrowPatch((xa, ymid), (xb, ymid), transform=fig.transFigure, arrowstyle="-|>",
                                   mutation_scale=12, lw=1.4, color="0.25"))
    fig.text((xa + xb) / 2, ymid + 0.05 / H_IN * 10, "warm start", fontsize=12, color="0.25", ha="center", va="bottom")
cx = (XS[2] + PW + 0.3) / W_IN
cax = fig.add_axes([cx, Y0 / H_IN, CB_W / W_IN, PH / H_IN])
sm = cm.ScalarMappable(norm=Normalize(vmin, vmax), cmap=cmap)
cb = fig.colorbar(sm, cax=cax)
cb.ax.tick_params(labelsize=12)
fig.savefig(OUT, dpi=200)
json.dump(dict(vmin=vmin, vmax=vmax, percentiles="5th and 60th of pooled liquid cells, each level equal weight",
               saturation=sat, saturation_all_levels=dict(above_vmax=sat_total_high, below_vmin=sat_total_low),
               crop_rows=[R0, R1], picked={k: v.name for k, v in picked.items()}),
          open(TALK / "data" / "s6_row1_r10.json", "w"), indent=1)
print("wrote", OUT)
