"""Posterior maps of L9 kLa_1T_10 for both BO arms from the real data (2026-11-08). Columns: arms. Rows: mean, uncertainty
(one shared colour bar per row, log colour for the mean). Marker size encodes level (L7 small -> L9 large), colour per
figstyle. Dashed white = tau_95 cap (model median). Also prints the chain model at its own L9 points (self-check).
  uv run python scripts/bo2d/fig_maps.py OUT_MAPS.png OUT_SURF.png"""
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts import figstyle as fs  # noqa: E402
from scripts.bo2d import bo, loop as L  # noqa: E402

lncap = np.log(json.loads(L.CFG.read_text())["cap"])
G = bo.grid(); n = bo.MAP_N
R, T = G[:, 0].reshape(n, n), G[:, 1].reshape(n, n)
MS = {7: 4, 8: 7, 9: 12}; MK = {7: "v", 8: "s", 9: "o"}
NAME = {"l9only": "L9 only", "chain": "L7 → L8 → L9"}
res = {}
for arm in ("l9only", "chain"):
    Mk, Mt, K, _ = L.fit_arm(arm)
    (mom, _), (mt, _) = Mk(G), Mt(G)
    res[arm] = dict(mu=np.exp(mom[9][0]), sd=np.sqrt(mom[9][1]), mt=mt[9][0], K=K)
    X9 = K[9][["rpm", "theta"]].to_numpy(float)
    m9 = np.exp(Mk(X9)[0][9][0])
    print(arm, "model vs data at own L9 points:", [f"{a:.1f}/{b:.1f}" for a, b in zip(m9, K[9].y)])
plt.rcParams.update(fs.rcparams())
fig, axs = plt.subplots(2, 2, figsize=(9.5, 7.6), sharex=True, sharey=True)
lo = min(r["mu"].min() for r in res.values()); hi = max(r["mu"].max() for r in res.values())
shi = max(r["sd"].max() for r in res.values())
for c, arm in enumerate(("l9only", "chain")):
    r = res[arm]
    for row, (Z, kw) in enumerate(((r["mu"], dict(cmap="viridis", norm=LogNorm(lo, hi))),
                                   (r["sd"], dict(cmap="magma", vmin=0, vmax=shi)))):
        ax = axs[row, c]
        im = ax.pcolormesh(R, T, Z.reshape(n, n), shading="auto", **kw)
        ax.contour(R, T, r["mt"].reshape(n, n), levels=[lncap], colors="w", linewidths=1.6, linestyles="--")
        for lev in (9, 8, 7):          # L7 last: it shares sites with L8
            if lev in r["K"] and len(r["K"][lev]):
                ax.plot(r["K"][lev].rpm, r["K"][lev].theta, MK[lev], ms=MS[lev], mfc=fs.LEVEL_COLOUR[lev],
                        mec="k" if lev == 7 else "w", mew=0.6, ls="none")
        if c == 1:
            fig.colorbar(im, ax=axs[row, :], location="right", shrink=0.9,
                         label="mean L9 kLa (1/h)" if row == 0 else "uncertainty, sd of ln kLa")
    axs[0, c].set_title(NAME[arm]); axs[1, c].set_xlabel("rpm")
for row in (0, 1):
    axs[row, 0].set_ylabel("θ_max (deg)")
h = [Line2D([], [], marker=MK[l], ms=MS[l], mfc=fs.LEVEL_COLOUR[l], mec="0.3", ls="none", label=f"L{l} run")
     for l in (7, 8, 9)] + [Line2D([], [], color="0.3", ls="--", lw=1.6, label="τ95 cap (feasible: left / below)")]
fig.legend(handles=h, frameon=False, loc="lower center", ncol=4, bbox_to_anchor=(0.45, -0.04))
fig.savefig(sys.argv[1], dpi=170, bbox_inches="tight")
# 3D surfaces of the mean, shared z
fig = plt.figure(figsize=(11, 4.8))
for c, arm in enumerate(("l9only", "chain")):
    ax = fig.add_subplot(1, 2, c + 1, projection="3d")
    r = res[arm]
    Z = np.log10(r["mu"]).reshape(n, n)
    ax.plot_surface(R, T, Z, cmap="viridis", vmin=np.log10(lo), vmax=np.log10(hi), lw=0, antialiased=True)
    ax.set_zlim(np.log10(lo), np.log10(hi)); ax.set_title(NAME[arm])
    ax.set_xlabel("rpm"); ax.set_ylabel("θ_max (deg)"); ax.set_zlabel("log10 mean L9 kLa (1/h)")
    ax.view_init(28, -128)
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=170, bbox_inches="tight")
print("wrote", sys.argv[1:])
