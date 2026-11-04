"""Slide 9 ("result 3", 2D campaign template): predicted mean and 2 sigma/mean over (rpm, theta) for the three gated QoIs.

Baseline data, round 3. The fit is the campaign.py fit (fit_qoi_info + max_crit), called read-only. campaign.py `fit` writes
only the maximum and a PNG (state/fit0.json, docs/campaign2d/fit0_*.png); it does not save the maps. So this script
recomputes the same fit and saves the maps it draws (data/s9_maps.npz).

Inputs (mtime from `stat -c '%y %n'`):
  scripts/campaign2d/state/gates.json   (2026-11-04 10:36; read only; FAILED-GATE runs dropped as in `fit`)
  scripts/campaign2d/state/fit0.json    (2026-11-04 10:36; for cross-check of the maxima only)
  experiments/multifidelity/hydro_dataset_v2.json (2026-10-08; hydro QoIs via campaign.load_db)
  experiments/multifidelity/kla_mix_dataset.csv   (kLa QoIs via campaign.load_db; runs/c2d_* rows)
  runs/c2d_* (none finished at the time of writing, so the campaign rows are empty)

Two rules stated on the figure: (1) edr_95 is NOT the 95th percentile. The fit uses edr_mean_t as a stand-in
(campaign.QOI_SOURCE), the same stand-in as `fit`. (2) The left map shows exp(predicted ln mean), i.e. a median, because
the model is fitted on ln y. Level colours follow figstyle: L8 blue, L9 green, L10 pink.
Output: figs/s9_heatmaps.png (slide 9), data/s9_maps.npz.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s9_heatmaps.png"
spec = importlib.util.spec_from_file_location("campaign", REPO / "scripts/campaign2d/campaign.py")
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)

LEVEL_COL = {8: "#0072B2", 9: "#009E73", 10: "#CC79A7"}
QOIS = [("tau_95", r"$\tau_{95}$ (Pa)"), ("edr_95", r"EDR (stand-in: $\langle\varepsilon\rangle_t$)"), ("kLa_1T_10", r"$k_La$ (1/h)")]
gates = json.load(open(REPO / "scripts/campaign2d/state/gates.json"))
fit0 = json.load(open(REPO / "scripts/campaign2d/state/fit0.json"))

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14, "axes.linewidth": 1.0})
fig, axes = plt.subplots(3, 2, figsize=(10.0, 5.6), dpi=200)
save = {}
for i, (q, lab) in enumerate(QOIS):
    d = C.exclude_failed_gate(C.load_db(q), gates)
    m, info = C.fit_qoi_info(d, q)
    res = C.max_crit(m)
    r, t = res["r"], res["t"]
    RR, TT = np.meshgrid(r, t, indexing="ij")
    G = np.column_stack([RR.ravel(), TT.ravel()])
    wm_all = []
    for s in range(0, len(G), 400):
        mu, sc, nu, w = m.predict_dist(G[s:s + 400])
        wm_all.append(w @ mu)
    pred = np.exp(np.concatenate(wm_all)).reshape(RR.shape)
    crit = res["grid"]
    print(f"{q}: n_L8={int((d.level == 8).sum())} n_L9+={int((d.level >= 9).sum())} max 2s/mean {res['max']:.3f} "
          f"at rpm {res['rpm']:.2f} theta {res['theta']:.2f}; fit0.json max {fit0[q]['max']:.3f}; "
          f"pred range {pred.min():.4g}..{pred.max():.4g} (ratio {pred.max() / pred.min():.1f})")
    save[q + "_pred"], save[q + "_crit"], save[q + "_r"], save[q + "_t"] = pred, crit, r, t
    use_log = pred.max() / pred.min() > 10
    norm = LogNorm(pred.min(), pred.max()) if use_log else Normalize(pred.min(), pred.max())
    ax = axes[i, 0]
    im = ax.pcolormesh(r, t, pred.T, shading="auto", cmap="viridis", norm=norm)
    ax2 = axes[i, 1]
    im2 = ax2.pcolormesh(r, t, crit.T, shading="auto", cmap="viridis", vmin=0, vmax=1.0)
    ax2.contour(r, t, crit.T, levels=[0.30], colors="k", linewidths=1.6)
    for a in (ax, ax2):
        s_ = d.copy()
        for lev in (8, 9, 10):
            ss = s_[s_.level == lev]
            a.scatter(ss.rpm, ss.theta, s=7, marker="o", facecolors=LEVEL_COL[lev], edgecolors="white", linewidths=0.5, zorder=3, clip_on=False)
        a.set_xlim(r.min(), r.max())
        a.set_ylim(t.min(), t.max())
        a.tick_params(labelsize=13)
    ax.set_ylabel(lab + "\n" + r"$\theta$ (deg)", fontsize=13)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.045)
    cb.ax.tick_params(labelsize=11)
    if i == 2:
        cax2 = fig.add_axes([0.935, 0.2, 0.012, 0.7])
        cb2 = fig.colorbar(im2, cax=cax2)
        cb2.ax.tick_params(labelsize=11)
        cb2.set_label(r"$2\sigma/\mu$", fontsize=12)
    if i == 2:
        ax.set_xlabel("rpm", fontsize=13)
        ax2.set_xlabel("rpm", fontsize=13)
    if i == 0:
        ax.set_title("predicted mean (median)", fontsize=13, pad=3)
        ax2.set_title(r"$2\sigma/\mu$ (black: 0.30)", fontsize=13, pad=3)
handles = [Line2D([], [], marker="o", ls="none", mfc=LEVEL_COL[lev], mec="white", ms=7) for lev in (8, 9, 10)]
fig.legend(handles, ["L8 runs", "L9 runs", "L10 runs"], loc="lower center", ncol=3, frameon=False, fontsize=13, bbox_to_anchor=(0.5, 0.0))
fig.subplots_adjust(left=0.13, right=0.915, top=0.91, bottom=0.175, hspace=0.45, wspace=0.5)
fig.savefig(OUT, dpi=200)
np.savez(TALK / "data" / "s9_maps.npz", **save)
print("wrote", OUT)
