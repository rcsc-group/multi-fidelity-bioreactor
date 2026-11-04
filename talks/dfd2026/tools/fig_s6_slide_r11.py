"""Slide 6, round 11: panels in the right half, one row per level (L8 top, L9, L10), equal aspect, shared limits.
Canvas = slide 1024 x 576 px (figure 10.24 x 5.76 in at 100 dpi). Equations and text are HTML on the slide.
Limits: 10th and 90th percentiles of pooled liquid |u| (data/limits_r11.json, from tools/limits_r11.py 90).
Output: figs/s6_slide_r11.png (1024 x 576, white background except the tanks, labels, arrows, colour bar)."""
import json, struct, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from scripts.fields import load_snapshot, snapshots, liquid_mask
from scripts.postprocess import _t_scales
TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s6_slide_r11.png"
LIM = json.load(open(TALK / "data" / "limits_r11.json"))
VMIN, VMAX = LIM["vmin"], LIM["vmax"]
RUNS = [("L8", "flow_l8_rpm32.5_th7"), ("L9", "scal_l9_rpm32.5_th7"), ("L10", "v_l10_rpm32.5_th7")]
COL = {"L8": "#0072B2", "L9": "#009E73", "L10": "#CC79A7"}
R0, R1 = 84, 156
W, H = 1024, 576
ROW_TOP = [62, 212, 362]
ROW_H = 150
PANEL_X0, PANEL_W = 588, 396
PANEL_H = PANEL_W * (R1 - R0) / 256
LABEL_X = 500
BAR_X0, BAR_W = 938, 10
BAR_Y0, BAR_Y1 = ROW_TOP[0] + (ROW_H - PANEL_H) / 2, ROW_TOP[2] + ROW_H - (ROW_H - PANEL_H) / 2
cmap = plt.get_cmap("cividis").copy(); cmap.set_bad("#e6e6e6")
fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
fig.patch.set_facecolor("white")
def fx(x): return x / W
def fy(y_top, h): return 1 - (y_top + h) / H
for i, (lab, run) in enumerate(RUNS):
    d = REPO / "runs" / run
    _, T = _t_scales(json.load(open(d / "params.json")))
    files = snapshots(d); tt = []
    for f in files:
        with open(f, "rb") as fh:
            fh.read(8); (t,) = struct.unpack("<d", fh.read(8)); tt.append(t)
    ph = (np.array(tt) / T) % 1.0
    dph = np.minimum(np.abs(ph - 0.25), 1 - np.abs(ph - 0.25))
    s = load_snapshot(files[int(np.argmin(dph))])
    m = liquid_mask(s); spd = np.where(m, np.hypot(s["ux"], s["uy"]), np.nan)
    k = s["n"] // 256
    field = spd[R0 * k: R1 * k, :]
    py = ROW_TOP[i] + (ROW_H - PANEL_H) / 2
    ax = fig.add_axes([fx(PANEL_X0), fy(py, PANEL_H), PANEL_W / W, PANEL_H / H])
    ax.imshow(field, origin="lower", cmap=cmap, vmin=VMIN, vmax=VMAX, interpolation="nearest", extent=[0, 256, R0, R1])
    ax.set_xlim(0, 256); ax.set_ylim(R0, R1); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_linewidth(0.9)
    fig.text(fx(PANEL_X0 - 16), fy(py + PANEL_H / 2 - 16, 32), lab, color=COL[lab], fontsize=30, ha="right", va="center")
# thin down arrows between rows, with "warm start" beside them
xa = (PANEL_X0 + PANEL_W * 0.5) / W
for i in range(2):
    y_from = ROW_TOP[i] + (ROW_H - PANEL_H) / 2 + PANEL_H + 6
    y_to = ROW_TOP[i + 1] + (ROW_H - PANEL_H) / 2 - 6
    fig.add_artist(matplotlib.patches.FancyArrowPatch((xa, 1 - y_from / H), (xa, 1 - y_to / H), transform=fig.transFigure,
                                                      arrowstyle="-|>", mutation_scale=12, lw=1.0, color="0.25"))
    fig.text((PANEL_X0 + PANEL_W * 0.5 + 10) / W, 1 - ((y_from + y_to) / 2) / H, "warm start", fontsize=16, color="0.25", ha="left", va="center")
fig.savefig(OUT, dpi=100, facecolor="white")
print("wrote", OUT, "limits", round(VMIN, 4), round(VMAX, 4), "panel h", round(PANEL_H, 1))
