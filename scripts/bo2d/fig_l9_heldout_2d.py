"""2D campaign held-out test (user, 2026-11-11 ~10:00: "show me how bad they actually are").
Data: truth snapshot v2 (real runs: L7 / L8 / L9 rows of the 2D campaign; c2d L9 dropped). Leave one L9 SITE out;
train on all L7, all L8 and the other L9 rows; predict ln y at the held-out site.
Models: L9-only GP (fit_single); 2-level L8 -> L9 (FullBayes, basis [1, mean f8, theta]); 3-level L7 -> L8 -> L9
(bo.fit_chain3, BO_THETA basis, = the campaign chain). 95 % band = latent + noise variance of the top level.
Figure 1 (fig_l9_heldout_2d.png): rows = QoI (tau_95, kLa_1T_10, dtmix_0.25), cols = model; x = held-out L9, y =
prediction; filled = theta = 7 site, hollow = off-line. Figure 2 (fig_probes_2d.png): training probes per QoI in
(rpm, theta), colour = level. Prints RMSE (ln) per model and site class, and 95 % coverage.
  uv run python scripts/bo2d/fig_l9_heldout_2d.py
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

bo.CHAIN_THETA = True
C = bo.C
S = ROOT / "scripts/bo2d/state"
QS = (("tau_95", r"$\tau_{95}$ (Pa)"), ("kLa_1T_10", r"$k_La$ (1/h)"), ("dtmix_0.25", r"$t_{\mathrm{mix}}$ (s)"))
MODELS = ("L9 only", "L8 → L9", "L7 → L8 → L9")


def load(q):
    d = pd.read_csv(S / f"truth_snapshot_v2_{q}.csv")
    d = d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)].copy()
    return {l: d[d.level == l].reset_index(drop=True) for l in (7, 8, 9)}


def predict(name, D, d9, X, q):
    if name == "L9 only":
        mom, _ = bo.single_moments(bo.fit_single(d9, q), X)
    elif name == "L8 → L9":
        m8 = bo.fit_single(D[8], q)
        Xt, yt = bo._xy(d9)
        m9 = bo._ThetaBasis(C.DOMAIN, lf_model=bo._Mean(m8), grid=C.hf_grid(C.noise_floor(q), float(np.std(yt)))).train(Xt, yt)
        mu, V, n = bo.latent(m9, X)
        mom = {9: (mu, bo.rho_of(m9) ** 2 * bo.latent(m8, X)[1] + V, n)}
    else:
        mom, _ = bo.chain_moments(bo.fit_chain3(D[7], D[8], d9, q), X)
    mu, V, n = mom[9]
    return float(np.ravel(mu)[0]), float(np.sqrt(np.ravel(V)[0] + np.ravel(n)[0]))


def main():
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(3, 3, figsize=(10.5, 10.0))
    fig2, ax2 = plt.subplots(1, 3, figsize=(12.0, 3.8), sharey=True)
    c9 = fs.level_colour(9)
    for row, (q, lab) in enumerate(QS):
        D = load(q)
        d9 = D[9]
        key = d9.rpm.round(2).astype(str) + "_" + d9.theta.round(2).astype(str)
        sites = key.unique()
        for L, ms in ((7, 3), (8, 5), (9, 8)):
            ax2[row].plot(D[L].rpm, D[L].theta, "o", ms=ms, color=fs.level_colour(L), mfc=fs.level_colour(L), alpha=0.8, ls="none",
                          label=f"L{L} ({len(D[L])})")
        ax2[row].set_title(lab); ax2[row].set_xlabel("rocking speed (rpm)"); ax2[row].grid(**fs.GRID_KW)
        ax2[row].legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=False)
        for col, name in enumerate(MODELS):
            yt, mu, sd, on = [], [], [], []
            for s in sites:
                te = key == s
                X = d9[te][["rpm", "theta"]].to_numpy(float)[:1]
                m, sdv = predict(name, D, d9[~te], X, q)
                yt.append(float(np.log(d9[te].y.astype(float)).mean())); mu.append(m); sd.append(sdv)
                on.append(bool(np.isclose(X[0, 1], 7.0)))
            yt, mu, sd, on = map(np.array, (yt, mu, sd, on))
            ax = axes[row, col]
            for sel, mfc in ((on, c9), (~on, "w")):
                ax.errorbar(np.exp(yt[sel]), np.exp(mu[sel]), yerr=[np.exp(mu[sel]) - np.exp(mu[sel] - 1.96 * sd[sel]),
                            np.exp(mu[sel] + 1.96 * sd[sel]) - np.exp(mu[sel])], fmt="o", ms=5, color=c9, mfc=mfc,
                            ecolor=c9, elinewidth=0.8, alpha=0.85)
            lo, hi = np.exp(min(yt.min(), (mu - 1.96 * sd).min())) * 0.9, np.exp(max(yt.max(), (mu + 1.96 * sd).max())) * 1.1
            ax.plot([lo, hi], [lo, hi], color="0.5", lw=0.8, ls=":")
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.grid(**fs.GRID_KW)
            if row == 0:
                ax.set_title(name)
            if col == 0:
                ax.set_ylabel(f"predicted {lab}")
            if row == 2:
                ax.set_xlabel("held-out L9")
            e = mu - yt; cov = np.mean(np.abs(e) <= 1.96 * sd)
            print(f"{q:11s} {name:14s} n={len(yt)} RMSE ln all {np.sqrt(np.mean(e ** 2)):.3f} | theta7 "
                  f"{np.sqrt(np.mean(e[on] ** 2)):.3f} (n {on.sum()}) | off-line {np.sqrt(np.mean(e[~on] ** 2)):.3f} "
                  f"(n {(~on).sum()}) | sd ln y9 {np.std(yt):.3f} | 95 % coverage {cov:.2f}", flush=True)
    ax2[0].set_ylabel(r"rocking angle $\theta$ (deg)")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="none", color=c9, mfc=c9), Line2D([], [], marker="o", ls="none", color=c9, mfc="w")]
    fig.legend(h, [r"$\theta$ = 7 site", "off-line site"], loc="upper center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(ROOT / "figures/bo2d/fig_l9_heldout_2d.png", dpi=170)
    fig2.tight_layout(); fig2.savefig(ROOT / "figures/bo2d/fig_probes_2d.png", dpi=170)
    print("wrote figures/bo2d/fig_l9_heldout_2d.png, fig_probes_2d.png")


if __name__ == "__main__":
    main()
