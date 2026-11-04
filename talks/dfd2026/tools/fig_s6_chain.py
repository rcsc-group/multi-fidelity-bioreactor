"""Slide 6 ("the chain"): one settled flow at three grid levels, phase-matched. Round 2: larger panels, sans fonts.

Data read directly (mtime from `stat -c '%y %n'`):
  .../multi-fidelity-bioreactor/runs/flow_l8_rpm32.5_th7/fields/snap_*.bin   (2026-10-11 03:18, L8, n=256)
  .../multi-fidelity-bioreactor/runs/scal_l9_rpm32.5_th7/fields/snap_*.bin   (2026-10-11 08:18, L9, n=512)
  .../multi-fidelity-bioreactor/runs/v_l10_rpm32.5_th7/fields/snap_*.bin     (2026-10-12 14:52, L10, n=1024)
  .../runs/*/params.json (period, fidelity)
  data/costs.json (from tools/cost_s8.py, sacct per run id; cost_table.tsv 2026-11-04 09:22)
Phase: frac(t / T), T from scripts.postprocess._t_scales. The frame nearest phase 0.25 is shown for each run.
Shared colour scale: speed in liquid, 99th percentile of L10.
Output: figs/s6_chain.png (960 px wide on the slide = 10.0 in at 96 px/in, so 14 pt text is 14 pt on the slide).
"""
import json
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors
from matplotlib.patches import FancyArrowPatch

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from scripts.fields import load_snapshot, snapshots, liquid_mask  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s6_chain.png"
OUT.parent.mkdir(parents=True, exist_ok=True)
TARGET_PHASE = 0.25
RUNS = [("L8", "flow_l8_rpm32.5_th7"), ("L9", "scal_l9_rpm32.5_th7"), ("L10", "v_l10_rpm32.5_th7")]
COL = {"L8": "#0072B2", "L9": "#009E73", "L10": "#CC79A7"}

lv = json.load(open(TALK / "data" / "costs.json"))["levels"]
LABEL = {
    "L8": f"median {lv['L8_own']:.0f} core-h per run\n(n = {lv['L8_n']}, cold, chain set)",
    "L9": f"median {lv['L9_incl']:.0f} core-h per run",
    "L10": f"median {lv['L10_incl']:,.0f} core-h per run\n(warm; run plus its source)",
}

picked = {}
for lab, run in RUNS:
    d = REPO / "runs" / run
    p = json.load(open(d / "params.json"))
    _, T = _t_scales(p)
    files = snapshots(d)
    tt = []
    for f in files:
        with open(f, "rb") as fh:
            fh.read(8)
            (t,) = struct.unpack("<d", fh.read(8))
        tt.append(t)
    tt = np.array(tt)
    ph = (tt / T) % 1.0
    dph = np.minimum(np.abs(ph - TARGET_PHASE), 1 - np.abs(ph - TARGET_PHASE))
    k = int(np.argmin(dph))
    picked[lab] = files[k]
    print(f"{lab}: T={T:.6f} frames={len(files)} chosen {files[k].name} t={tt[k]:.5f} phase={ph[k]:.4f}")

snaps = {lab: load_snapshot(picked[lab]) for lab, _ in RUNS}
speeds = {lab: np.where(liquid_mask(snaps[lab]), np.hypot(snaps[lab]["ux"], snaps[lab]["uy"]), np.nan) for lab in snaps}
vmax = float(np.nanpercentile(speeds["L10"], 99))

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14})
W_IN, H_IN = 10.0, 3.1
fig = plt.figure(figsize=(W_IN, H_IN), dpi=200)
PW, PH, Y0 = 2.9, 2.03, 0.62
XS = [0.08, 3.56, 7.04]
for i, (lab, _) in enumerate(RUNS):
    ax = fig.add_axes([XS[i] / W_IN, Y0 / H_IN, PW / W_IN, PH / H_IN])
    s = snaps[lab]
    n = s["n"]
    ax.imshow(speeds[lab], origin="lower", extent=[-0.5, 0.5, -0.5, 0.5], cmap="cividis", vmin=0, vmax=vmax, interpolation="nearest")
    ax.contour(np.linspace(-0.5, 0.5, n), np.linspace(-0.5, 0.5, n), s["f"], levels=[0.5], colors="k", linewidths=0.7)
    m = liquid_mask(s)
    rows = np.where(m.any(axis=1))[0]
    cols = np.where(m.any(axis=0))[0]
    cc = lambda idx: -0.5 + (idx + 0.5) / n
    pad = 0.01
    ax.set_xlim(cc(cols.min()) - pad, cc(cols.max()) + pad)
    ax.set_ylim(cc(rows.min()) - pad, cc(rows.max()) + pad)
    ax.set_aspect('auto')
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.9)
    fig.text((XS[i] + PW / 2) / W_IN, (Y0 + PH + 0.1) / H_IN, lab, color=COL[lab], fontsize=18, ha="center", va="bottom")
    fig.text((XS[i] + PW / 2) / W_IN, (Y0 - 0.1) / H_IN, LABEL[lab], ha="center", va="top", fontsize=14, color="0.15", linespacing=1.15)

ymid = (Y0 + PH / 2) / H_IN
for i in range(2):
    x0 = (XS[i] + PW + 0.06) / W_IN
    x1 = (XS[i + 1] - 0.06) / W_IN
    fig.add_artist(FancyArrowPatch((x0, ymid), (x1, ymid), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=18, lw=1.6, color="0.25"))
    fig.text((x0 + x1) / 2, ymid + 0.09, "warm\nstart", ha="center", va="bottom", fontsize=14, color="0.25", linespacing=1.0)

fig.savefig(OUT, dpi=200)
print("wrote", OUT, "vmax", round(vmax, 4))
