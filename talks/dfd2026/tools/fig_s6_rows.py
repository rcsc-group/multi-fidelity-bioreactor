"""Slide 6 ("the chain"), round 6: three stacked rows (L8, L9, L10), EQUAL ASPECT (no stretch).

Same data and frame choice as tools/fig_s6_chain.py (phase 0.25 of the rocking cycle, same source frames):
  .../multi-fidelity-bioreactor/runs/flow_l8_rpm32.5_th7/fields/snap_*.bin   (2026-10-11 03:18)
  .../multi-fidelity-bioreactor/runs/scal_l9_rpm32.5_th7/fields/snap_*.bin   (2026-10-11 08:18)
  .../multi-fidelity-bioreactor/runs/v_l10_rpm32.5_th7/fields/snap_*.bin     (2026-10-12 14:52)
  data/costs.json (chain-set medians: L8 50, L9 550, L10 1,139 core-h per run)
Panel height is fixed by the slide height. Panel width follows from the domain aspect (equal aspect, no stretch).
Output: figs/s6_rows.png
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

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from scripts.fields import load_snapshot, snapshots, liquid_mask  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s6_rows.png"
TARGET_PHASE = 0.25
RUNS = [("L8", "flow_l8_rpm32.5_th7"), ("L9", "scal_l9_rpm32.5_th7"), ("L10", "v_l10_rpm32.5_th7")]
COL = {"L8": "#0072B2", "L9": "#009E73", "L10": "#CC79A7"}
lv = json.load(open(TALK / "data" / "costs.json"))["levels"]
COST_TXT = {"L8": f"{lv['L8_own']:.0f} core-h", "L9": f"{lv['L9_incl']:.0f} core-h", "L10": f"{lv['L10_incl']:,.0f} core-h"}

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
speeds = {lab: np.where(liquid_mask(snaps[lab]), np.hypot(snaps[lab]["ux"], snaps[lab]["uy"]), np.nan) for lab in snaps}
vmax = float(np.nanpercentile(speeds["L10"], 99))

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14})
W_IN, H_IN = 10.0, 4.92
PH = 1.36                     # panel height (in), set by the slide height
GAP = 0.34                    # vertical gap between rows (holds the down arrow)
fig = plt.figure(figsize=(W_IN, H_IN), dpi=200)
rows_top = H_IN - 0.06
X0 = 1.25                     # panel left edge (in), after the level label
for i, (lab, _) in enumerate(RUNS):
    s = snaps[lab]
    n = s["n"]
    m = liquid_mask(s)
    rr = np.where(m.any(axis=1))[0]
    cc_ = np.where(m.any(axis=0))[0]
    cc = lambda idx: -0.5 + (idx + 0.5) / n
    pad = 0.01
    x_lo, x_hi = cc(cc_.min()) - pad, cc(cc_.max()) + pad
    y_lo, y_hi = cc(rr.min()) - pad, cc(rr.max()) + pad
    aspect = (x_hi - x_lo) / (y_hi - y_lo)          # width / height of the liquid domain, cell units
    PW = PH * aspect
    y0 = rows_top - (i + 1) * PH - i * GAP
    ax = fig.add_axes([X0 / W_IN, y0 / H_IN, PW / W_IN, PH / H_IN])
    ax.imshow(speeds[lab], origin="lower", extent=[-0.5, 0.5, -0.5, 0.5], cmap="cividis", vmin=0, vmax=vmax, interpolation="nearest")
    ax.contour(np.linspace(-0.5, 0.5, n), np.linspace(-0.5, 0.5, n), s["f"], levels=[0.5], colors="k", linewidths=0.7)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect("equal", adjustable="box")      # equal aspect: no stretch
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.9)
    yc = y0 + PH / 2
    fig.text(0.1 / W_IN, yc / H_IN, lab, color=COL[lab], fontsize=20, ha="left", va="center")
    fig.text((X0 + PW + 0.3) / W_IN, yc / H_IN, COST_TXT[lab], fontsize=16, color="0.15", ha="left", va="center")
    if i < 2:
        xa = (X0 + PW * 0.5) / W_IN
        ya_top = y0 - 0.03
        ya_bot = y0 - GAP + 0.03
        fig.add_artist(FancyArrowPatch((xa, ya_top / H_IN), (xa, ya_bot / H_IN), transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=16, lw=1.5, color="0.25"))
        fig.text((X0 + PW * 0.5 + 0.15) / W_IN, ((ya_top + ya_bot) / 2) / H_IN, "warm start", fontsize=14, color="0.25", ha="left", va="center")
fig.savefig(OUT, dpi=200)
print("wrote", OUT, "vmax", round(vmax, 4))
