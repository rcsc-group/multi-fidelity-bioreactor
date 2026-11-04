"""Slide 6 right side, round 8: vertical speed profile from the bottom wall, L8 / L9 / L10 overlaid.

Same frames as tools/fig_s6_chain.py (phase 0.25): runs/flow_l8_rpm32.5_th7, runs/scal_l9_rpm32.5_th7,
runs/v_l10_rpm32.5_th7 (fields/snap_*.bin). Line: column COL0 = 33 (L8 grid; L9 and L10 use the same physical column),
from the bottom wall (first liquid cell) up to row ROW0 = 125 (L8 grid).
x axis: distance from the bottom wall, starting at 0, limited to the data range.
End labels use the slide 8 rule (tools/cost_rule_r7.py): each run charged its own cost plus its source when the source is
not in the set; the set is the full chain set (40 L8, 16 L9, 13 L10). Labels show the mean charge per run, as on slide 8.
Output: figs/s6_profile.png (4.9 x 3.75 in = 470 x 360 px on the slide), data/s6_profile_r8.json
"""
import csv
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
import scripts.test_h1_mfbml as H1  # noqa: E402
import scripts.test_h2_mfbml as H2  # noqa: E402

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s6_profile.png"
RUNS = [("L8", "flow_l8_rpm32.5_th7"), ("L9", "scal_l9_rpm32.5_th7"), ("L10", "v_l10_rpm32.5_th7")]
COL = {"L8": "#0072B2", "L9": "#009E73", "L10": "#CC79A7"}
COL0, ROW0 = 33, 125
TARGET_PHASE = 0.25

# ---- cost labels (slide 8 rule, full chain set)
COST = {}
with open("/oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        COST[r["run"]] = float(r["corehc_all"])
l8_set = {H1.lf_run(r, t) for r in H1.RPMS for t in H1.LF_THETAS} & set(COST)
l9_names = [H2.l9_run(r, t) for r, t in H2.l9_points()]
l10_names = [H1.hf_run(r, t) for r, t in H1.hf_points()]
chain_set = l8_set | set(l9_names)
mean_L8 = float(np.mean([COST[n] for n in l8_set]))
mean_L9 = float(np.mean([COST[n] + 0.0 for n in l9_names]))          # sources of warm L9 runs are L8 runs in the set
mean_L10 = float(np.mean([COST[n] for n in l10_names]))               # sources of L10 runs are in the set (x8: L8; xlevel: fig9 L9 in set)
labels = {"L8": f"{mean_L8:,.0f} core-h", "L9": f"{mean_L9:,.0f} core-h", "L10": f"{mean_L10:,.0f} core-h"}
print("labels (mean charge per run, full chain set):", labels)

# ---- fields
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
curves, meta = {}, {}
for lab in snaps:
    s = snaps[lab]
    n = s["n"]
    k = n // 256
    m = liquid_mask(s)
    spd = np.where(m, np.hypot(s["ux"], s["uy"]), np.nan)
    col = COL0 * k
    colv = spd[:, col]
    r_top = ROW0 * k
    ok = np.where(np.isfinite(colv[: r_top + 1]))[0]
    r0 = int(ok[0])
    idx = np.arange(r0, r_top + 1)
    dist = (idx - r0 + 0.5) / n            # distance from the first liquid cell (the wall), domain units
    vals = colv[r0: r_top + 1]
    curves[lab] = (dist, vals)
    meta[lab] = dict(n=n, first_liquid_row=r0, top_row=r_top, col=col, file=picked[lab].name)
    print(lab, meta[lab], "dist max", round(float(dist[-1]), 4))

dmax = max(c[0][-1] for c in curves.values())
x_phys = -0.5 + (COL0 + 0.5) / 256
y_wall = -0.5 + (meta["L8"]["first_liquid_row"] + 0.5) / 256
y_top = -0.5 + (ROW0 + 0.5) / 256
print(f"plotted line: x = {x_phys:.4f}, from the wall at y = {y_wall:.4f} to y = {y_top:.4f} (domain units, origin at centre)")

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14, "axes.linewidth": 1.1,
                     "xtick.direction": "in", "ytick.direction": "in"})
fig, ax = plt.subplots(figsize=(4.9, 3.75), dpi=200)
for lab, _ in RUNS:
    dist, vals = curves[lab]
    ax.plot(dist, vals, color=COL[lab], ls=":", lw=1.7)
    ypos = {"L8": 0.44, "L9": 0.37, "L10": 0.30}[lab]
    ax.text(dmax, ypos, f"{lab}  {labels[lab]}", color=COL[lab], fontsize=13, ha="right", va="center")
ax.set_xlim(0, dmax)
ax.set_ylim(-0.01, 0.58)
ax.set_xlabel("distance from wall (domain units)", fontsize=14)
ax.set_ylabel(r"speed $|u|$ (nondimensional)", fontsize=14)
ax.tick_params(labelsize=13)
ax.grid(which="major", ls=":", alpha=0.4)
fig.subplots_adjust(left=0.17, right=0.97, top=0.97, bottom=0.17)
fig.savefig(OUT, dpi=200)
(TALK / "data" / "s6_profile_r8.json").write_text(json.dumps(dict(labels=labels, mean_charge=dict(L8=mean_L8, L9=mean_L9, L10=mean_L10),
                                                                  meta=meta, x_phys=x_phys, y_wall=y_wall, y_top=y_top, dmax=dmax), indent=1))
print("wrote", OUT)
