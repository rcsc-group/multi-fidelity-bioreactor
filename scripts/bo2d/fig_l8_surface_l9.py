"""L8 surrogate surface with L9 scatter (user, 2026-11-11 ~10:05). Per QoI (tau_95, kLa_1T_10, dtmix_0.25): GP
(FullBayes, ordinary) fitted to ALL L8 rows of truth snapshot v2; filled contours of its mean (ln, own colour scale);
L8 probes as small black dots; L9 runs as large circles coloured by their L9 value on their OWN colour scale (same
colormap), so the two patterns compare by shape, not by absolute value (L8 kLa ~3.7x L9). Grid predicted in chunks.
  uv run python scripts/bo2d/fig_l8_surface_l9.py -> figures/bo2d/fig_l8_surface_l9.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts.bo2d import bo  # noqa: E402

C = bo.C
S = ROOT / "scripts/bo2d/state"
QS = (("tau_95", r"$\tau_{95}$ (Pa)"), ("kLa_1T_10", r"$k_La$ (1/h)"), ("dtmix_0.25", r"$t_{\mathrm{mix}}$ (s)"))
NR, NT = 46, 31

plt.rcParams.update(fs.rcparams())
fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.4), sharey=True, gridspec_kw={"wspace": 0.45})
r = np.linspace(*C.RPM, NR); t = np.linspace(*C.THETA, NT)
RR, TT = np.meshgrid(r, t, indexing="ij"); G = np.column_stack([RR.ravel(), TT.ravel()])
for ax, (q, lab) in zip(axes, QS):
    d = pd.read_csv(S / f"truth_snapshot_v2_{q}.csv")
    d = d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]
    d8, d9 = d[d.level == 8], d[d.level == 9]
    m = bo.fit_single(d8, q)
    mu = np.concatenate([np.ravel(m.predict_mean(G[i:i + 200])) for i in range(0, len(G), 200)]).reshape(NR, NT)
    cs = ax.contourf(RR, TT, np.exp(mu), levels=14, cmap="viridis", alpha=0.9)
    cb = fig.colorbar(cs, cax=ax.inset_axes([1.03, 0.0, 0.04, 1.0])); cb.set_label(f"L8 surrogate {lab}")
    ax.plot(d8.rpm, d8.theta, ".", color="k", ms=3, alpha=0.6)
    sc = ax.scatter(d9.rpm, d9.theta, c=d9.y.astype(float), cmap="viridis", s=90, edgecolors="w", linewidths=1.4,
                    norm=matplotlib.colors.LogNorm())
    cb2 = fig.colorbar(sc, cax=ax.inset_axes([0.0, 1.06, 1.0, 0.05]), orientation="horizontal"); cb2.set_label(f"L9 {lab}"); cb2.ax.xaxis.set_label_position("top"); cb2.ax.xaxis.set_ticks_position("top")
    ax.set_xlabel("rocking speed (rpm)"); ax.set_xlim(*C.RPM); ax.set_ylim(*C.THETA)
axes[0].set_ylabel(r"rocking angle $\theta$ (deg)")
fig.subplots_adjust(left=0.05, right=0.94, top=0.8, bottom=0.12)
out = ROOT / "figures/bo2d/fig_l8_surface_l9.png"
fig.savefig(out, dpi=170)
print("wrote", out)
