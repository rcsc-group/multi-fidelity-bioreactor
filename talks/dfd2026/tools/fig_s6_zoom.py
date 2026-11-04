"""Slide 6 (round 6, B): right-hand column, the same zoomed crop at L8, L9 and L10, with cell edges.

Frames: same phase-matched snapshots as tools/fig_s6_chain.py (phase 0.25 of the rocking cycle):
  runs/flow_l8_rpm32.5_th7/fields (2026-10-11 03:18), runs/scal_l9_rpm32.5_th7/fields (2026-10-11 08:18),
  runs/v_l10_rpm32.5_th7/fields (2026-10-12 14:52).
Crop: where |u_L8 - u_L10| (both block-averaged to the L8 grid, 256 x 256) is largest after a 9 x 9 box smooth, and the
window lies fully in the liquid. Window = 32 L8 cells per side; that is 64 L9 cells and 128 L10 cells.
Output: figs/s6_zoom.png (4.9 x 5.0 in, rendered at 96 px/in on the slide, so 14 pt labels are 14 pt).
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
OUT = TALK / "figs" / "s6_zoom.png"
RUNS = [("L8", "flow_l8_rpm32.5_th7"), ("L9", "scal_l9_rpm32.5_th7"), ("L10", "v_l10_rpm32.5_th7")]
COL = {"L8": "#0072B2", "L9": "#009E73", "L10": "#CC79A7"}
WIN = 32            # window in L8 cells per side
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
    print(lab, picked[lab].name)

snaps = {lab: load_snapshot(picked[lab]) for lab, _ in RUNS}
spd, msk = {}, {}
for lab in snaps:
    s = snaps[lab]
    m = liquid_mask(s)
    spd[lab] = np.hypot(s["ux"], s["uy"])
    msk[lab] = m


def block(a, m, k):
    """Block-average a (n x n) to (n/k x n/k), averaging liquid cells only. NaN where liquid fraction < 0.5."""
    n = a.shape[0] // k
    A = (a * m).reshape(n, k, n, k).sum(axis=(1, 3))
    M = m.reshape(n, k, n, k).sum(axis=(1, 3))
    out = np.where(M >= 0.5 * k * k, A / np.maximum(M, 1), np.nan)
    return out, (M >= 0.5 * k * k)


s8, ok8 = block(spd["L8"], msk["L8"], 1)
s10, ok10 = block(spd["L10"], msk["L10"], 4)          # 1024 -> 256
diff = np.where(ok8 & ok10, np.abs(s8 - s10), 0.0)
k = 9
pad = np.pad(diff, k // 2, mode="constant")
sm = np.zeros_like(diff)
for i in range(k):
    for j in range(k):
        sm += pad[i:i + diff.shape[0], j:j + diff.shape[1]]
sm /= k * k
n8 = s8.shape[0]
# candidate centres: window fully inside the liquid of both L8 and L10
full_ok = np.zeros_like(ok8)
h = WIN // 2
for i in range(h, n8 - h):
    for j in range(h, n8 - h):
        full_ok[i, j] = ok8[i - h:i + h, j - h:j + h].all() and ok10[i - h:i + h, j - h:j + h].all()
cand = np.where(full_ok, sm, -np.inf)
ci, cj = np.unravel_index(np.argmax(cand), cand.shape)
print(f"crop centre (L8 cell) row {ci} col {cj}; mean |du| there {sm[ci, cj]:.4f}; max over domain {np.nanmax(np.where(np.isfinite(cand), cand, np.nan)):.4f}")
vmax = float(np.nanpercentile(np.where(msk["L10"], spd["L10"], np.nan), 99))

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14})
W_IN, H_IN = 4.9, 5.0
fig = plt.figure(figsize=(W_IN, H_IN), dpi=200)
PS = 1.42                 # panel side (in)
X0 = 0.12
rows_top = H_IN - 0.12
for idx, (lab, _) in enumerate(RUNS):
    n_lv = snaps[lab]["n"]
    f = n_lv // 256                          # cells of this level per L8 cell
    r0, r1 = (ci - h) * f, (ci + h) * f
    c0, c1 = (cj - h) * f, (cj + h) * f
    field = np.where(msk[lab], spd[lab], np.nan)[r0:r1, c0:c1]
    y0 = rows_top - (idx + 1) * PS - idx * 0.22
    ax = fig.add_axes([X0 / W_IN, y0 / H_IN, PS / W_IN, PS / H_IN])
    xs = np.arange(c0, c1 + 1)
    ys = np.arange(r0, r1 + 1)
    ec_alpha, ec_lw = {1: (0.85, 0.7), 2: (0.5, 0.35), 4: (0.3, 0.15)}[f]
    ax.pcolormesh(xs, ys, field, cmap="cividis", vmin=0, vmax=vmax, shading="flat",
                  edgecolors=(0.1, 0.1, 0.1, ec_alpha), linewidth=ec_lw)
    ax.set_xlim(c0, c1)
    ax.set_ylim(r0, r1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.9)
    yc = y0 + PS / 2
    cells = n_lv
    lab_txt = {"L8": f"L8   2$^8$ cells per side", "L9": f"L9   2$^9$ cells per side", "L10": f"L10   2$^{{10}}$ cells per side"}[lab]
    fig.text((X0 + PS + 0.2) / W_IN, (yc + 0.12) / H_IN, lab_txt, color=COL[lab], fontsize=14, ha="left", va="center")
    cost_txt = {"L8": "50 core-h per run", "L9": "550 core-h per run", "L10": "1,139 core-h per run"}[lab]
    fig.text((X0 + PS + 0.2) / W_IN, (yc - 0.14) / H_IN, cost_txt, color="0.15", fontsize=14, ha="left", va="center")
fig.savefig(OUT, dpi=200)
print("wrote", OUT, "vmax", round(vmax, 4))
