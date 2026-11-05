"""Backup validation slide, top row: ours (L10, theta 7) vs Kim, with a per-point aleatoric band.

Aleatoric sigma per point, from that ONE run only (heteroscedastic):
  - per rocking cycle c (after dropping the first cycle after the restart) store max/min of the signed
    spatial-mean shear (col 10) and max of the spatial-mean EDR (COL_EDISS);
  - the plotted statistic is taken over N_CYC = 3 cycles (hydro_dataset_v2.qois):
        tau_amp = 0.5 (max_c max10 - min_c min10),  edr_max = max_c edr;
  - bootstrap: draw N_CYC cycles with replacement from the run's cycles, recompute the statistic,
    B = 4000 draws; sigma = sd of ln(statistic).
  - trend test: OLS slope of ln(per-cycle value) vs cycle index; |t| > 2 flags a drift (transient,
    not noise). Reported in data/s7a_band.json and the speaker notes, not drawn (hollow = mean in the deck grammar).
Band = exp(ln y +- 2 sigma) around our points (shaded, linear interpolation between rpm points).
No text in the figure except axes, ticks and legend.
Writes figs/s7a_band.png and data/s7a_band.json.
"""
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
sys.argv = [sys.argv[0]]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
warnings.filterwarnings("ignore")
import scripts.hydro_dataset_v2 as H  # noqa: E402

PINK, KIM = "#CC79A7", "black"
RPMS = [17.5, 20, 22.5, 25, 27.5, 30, 32.5]
B = 4000
rng = np.random.default_rng(0)
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14, "axes.linewidth": 1.1,
                     "xtick.direction": "in", "ytick.direction": "in"})


def cycles(run):
    d = REPO / "runs" / run
    p = json.load(open(d / "params.json"))
    arr = H._load_dat(d / "shear_stress.dat")
    T_bio, T = H._t_scales(p)
    a = p.get("geometry", {}).get("a", 0.25)
    U = a / T_bio
    scale, escale = H.RHO * U**2, H.RHO * U**3 / a
    cyc = arr[:, 1] / T
    c0 = math.ceil(p.get("t_checkpoint", 0.0) / T - 1e-6) + 1
    c1 = math.floor(cyc.max() + 1e-6)
    rows = []
    for s in range(c0, c1):
        m = (cyc >= s) & (cyc < s + 1)
        if m.sum() >= 10:
            rows.append((arr[m, 10].max() * scale, arr[m, 10].min() * scale, arr[m, H.COL_EDISS].max() * escale))
    return np.array(rows)


def boot(c):
    n = len(c)
    idx = rng.integers(0, n, size=(B, H.N_CYC))
    tau = 0.5 * (c[idx, 0].max(axis=1) - c[idx, 1].min(axis=1))
    edr = c[idx, 2].max(axis=1)
    return float(np.log(tau).std(ddof=1)), float(np.log(edr).std(ddof=1))


def trend_t(v):
    x = np.arange(len(v), dtype=float)
    y = np.log(v)
    X = np.c_[np.ones_like(x), x]
    beta, res, *_ = np.linalg.lstsq(X, y, rcond=None)
    dof = len(v) - 2
    if dof < 1:
        return float("nan")
    s2 = float(((y - X @ beta) ** 2).sum()) / dof
    se = math.sqrt(s2 / ((x - x.mean()) ** 2).sum())
    return float(beta[1] / se) if se > 0 else float("inf")


out = {}
for r in RPMS:
    run = f"l10_fig13a_xlevel_rpm{r:g}"
    q = H.qois(run)
    c = cycles(run)
    s_tau, s_edr = boot(c)
    t_tau = trend_t(0.5 * (c[:, 0] - c[:, 1]))
    t_edr = trend_t(c[:, 2])
    out[r] = dict(run=run, n_cyc=int(len(c)), tau_amp_t=q["tau_amp_t"], edr_max_t=q["edr_max_t"],
                  sigma_ln_tau=s_tau, sigma_ln_edr=s_edr, trend_t_tau=t_tau, trend_t_edr=t_edr)
    print(f"{r:5g} rpm  cycles {len(c)}  sigma ln tau_amp {s_tau:.4f} (trend t {t_tau:+.2f})  "
          f"sigma ln edr_max {s_edr:.4f} (trend t {t_edr:+.2f})")
json.dump({str(k): v for k, v in out.items()}, open(TALK / "data/s7a_band.json", "w"), indent=1)

kim = pd.read_csv(REPO / "experiments/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv", skiprows=[1])
kim = kim[(kim["RPM"] >= 17.5) & (kim["RPM"] <= 32.5)].sort_values("RPM")
rk = kim["RPM"].to_numpy(float)

fig, axs = plt.subplots(1, 2, figsize=(9.6, 2.55), dpi=200)
for ax, key, skey, tkey, kv, ylab, ttl, mk in [
    (axs[0], "tau_amp_t", "sigma_ln_tau", "trend_t_tau", kim["tau_liq_mean"].to_numpy(float), "Pa", r"$\tau_{\mathrm{amp}}$ (Pa)", "o"),
    (axs[1], "edr_max_t", "sigma_ln_edr", "trend_t_edr", kim["Ediss_liq_mean"].to_numpy(float), "W m$^{-3}$", r"$\varepsilon_{\max}$ (W m$^{-3}$)", "s"),
]:
    y = np.array([out[r][key] for r in RPMS])
    s = np.array([out[r][skey] for r in RPMS])
    drift = np.array([abs(out[r][tkey]) > 2 for r in RPMS])
    ax.fill_between(RPMS, y * np.exp(-2 * s), y * np.exp(2 * s), color=PINK, alpha=0.25, lw=0)
    ax.plot(rk, kv, color=KIM, marker=mk, ms=6.5, mfc=KIM, mec=KIM, lw=1.3, ls="-")
    ax.plot(RPMS, y, color=PINK, lw=1.1, ls=":")
    ax.plot(RPMS, y, ls="none", marker=mk, ms=5.5, mfc=PINK, mec=PINK)
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 2.0, 5.0)))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("rpm", fontsize=14)
    ax.set_title(ttl, fontsize=14, pad=8)
    ax.grid(which="major", ls=":", alpha=0.4)
    ax.tick_params(labelsize=13)
    ax.set_xticks([20, 25, 30])
handles = [Line2D([], [], color=KIM, marker="o", ms=6.5, mfc=KIM, lw=1.3),
           Line2D([], [], color=PINK, marker="o", ms=5.5, mfc=PINK, lw=1.1, ls=":"),
           Patch(color=PINK, alpha=0.25, lw=0)]
fig.legend(handles, ["Kim et al. 2025", "ours, L10", r"$\pm 2\sigma$, cycle to cycle"], loc="upper center", ncol=3,
           frameon=False, fontsize=13, bbox_to_anchor=(0.5, 1.01), handlelength=2.0, columnspacing=1.2)
fig.subplots_adjust(left=0.135, right=0.985, top=0.62, bottom=0.2, wspace=0.42)
fig.savefig(TALK / "figs/s7a_band.png", dpi=200)
print("wrote", TALK / "figs/s7a_band.png")
