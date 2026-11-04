"""Slide 7 ("result 1"), round 2. Two figures, each two panels of half-slide width.

  figs/s7a_kim.png   (a) tau_amp_t and edr_max_t vs rpm at theta 7: ours (L10 training run, plus the second L10 run at 25 rpm) vs Kim
  figs/s7b_predict.png (b) predicted vs observed with 95 % bands: LOO (13 runs) and theta 5.5 test set (6 runs)

Data read directly (mtime from `stat -c '%y %n'`):
  .../multi-fidelity-bioreactor/experiments/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv   (2026-07-21 11:11; digitised Kim 2025 Fig 13,
      theta 7. tau_liq_mean -> our tau_amp_t (Pa), Ediss_liq_mean -> our edr_max_t (W m^-3). Mapping: diary 2026-10-22.)
  .../experiments/multifidelity/h3_band_check_kim.json  (2026-10-23 09:34; LOO points, variant B = L9 variance propagated)
  .../experiments/multifidelity/test_h4_testset_kim.json (2026-10-23 09:34; means checked against the recomputed test predictions)
  .../runs/l10_fig13a_xlevel_rpm<r>, runs/x8_l10_rpm25_th7, runs/x8_l10_rpm<r>_th5.5 (shear_stress.dat via scripts/hydro_dataset_v2.qois)
Kim et al. (2025) = arXiv:2504.05421.
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
sys.argv = [sys.argv[0]]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
warnings.filterwarnings("ignore")
import scripts.test_h4_testset as T  # noqa: E402
from scripts.hydro_dataset_v2 import qois  # noqa: E402

OUT = TALK / "figs"
OUT.mkdir(parents=True, exist_ok=True)
PINK = "#CC79A7"
KIM = "black"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14, "axes.linewidth": 1.1,
                     "xtick.direction": "in", "ytick.direction": "in"})
GRID = {"which": "major", "ls": ":", "alpha": 0.4}

# ------------------------------------------------------------------ (a)
kim = pd.read_csv(REPO / "experiments/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv", skiprows=[1])
kim = kim[(kim["RPM"] >= 17.5) & (kim["RPM"] <= 32.5)].sort_values("RPM")
rk = kim["RPM"].to_numpy(float)
kim_tau = kim["tau_liq_mean"].to_numpy(float)
kim_edr = kim["Ediss_liq_mean"].to_numpy(float)

hd = json.load(open(REPO / "experiments/multifidelity/h3_band_check_kim.json"))
rpm_th7 = sorted({p["rpm"] for p in hd["tau_95_t"]["points_B"] if p["theta"] == 7.0})
ours = {r: qois(T.H1.hf_run(r, 7.0)) for r in rpm_th7} if hasattr(T, "H1") else {}
if not ours:
    import scripts.test_h1_mfbml as H1  # noqa: E402
    ours = {r: qois(H1.hf_run(r, 7.0)) for r in rpm_th7}
second = qois("x8_l10_rpm25_th7")
for r in rpm_th7:
    print(f"(a) training run {r:g}: tau_amp={ours[r]['tau_amp_t']:.6g} Pa, edr_max={ours[r]['edr_max_t']:.6g} W/m3")
print(f"(a) second run 25 rpm (x8_l10_rpm25_th7): tau_amp={second['tau_amp_t']:.6g}, edr_max={second['edr_max_t']:.6g}")
for r, t, e in zip(rk, kim_tau, kim_edr):
    print(f"(a) Kim {r:g} rpm: tau {t:.6g} Pa, edr {e:.6g} W/m3")
for r in rpm_th7:
    k = int(np.where(rk == r)[0][0])
    print(f"(a) {r:g}: ours/Kim-1 tau {100*(ours[r]['tau_amp_t']/kim_tau[k]-1):+.1f} %, edr {100*(ours[r]['edr_max_t']/kim_edr[k]-1):+.1f} %")
k25 = int(np.where(rk == 25.0)[0][0])
print(f"(a) second 25 rpm run/Kim-1 tau {100*(second['tau_amp_t']/kim_tau[k25]-1):+.1f} %, edr {100*(second['edr_max_t']/kim_edr[k25]-1):+.1f} %")

figa, axa = plt.subplots(1, 2, figsize=(9.6, 2.55), dpi=200)
for ax, key, kv, ylab, ttl in [
    (axa[0], "tau_amp_t", kim_tau, "Pa", r"$\tau_{\mathrm{amp}}$ (Pa)"),
    (axa[1], "edr_max_t", kim_edr, "W m$^{-3}$", r"$\varepsilon_{\max}$ (W m$^{-3}$)"),
]:
    mk = "o" if key.startswith("tau") else "s"
    ax.plot(rk, kv, color=KIM, marker=mk, ms=6.5, mfc=KIM, mec=KIM, lw=1.3, ls="-")
    ax.plot(rpm_th7, [ours[r][key] for r in rpm_th7], color=PINK, marker=mk, ms=5.5, mfc=PINK, mec=PINK, lw=1.1, ls=":")
    ax.plot([25.0], [second[key]], color=PINK, marker=mk, ms=6.5, mfc="white", mec=PINK, mew=1.4, ls="none")
    ax.set_yscale("log")
    ax.set_xlabel("rpm", fontsize=14)
    ax.set_ylabel(ylab, fontsize=14)
    ax.set_title(ttl, fontsize=14, pad=3)
    ax.grid(**GRID)
    ax.tick_params(labelsize=13)
    ax.set_xticks([20, 25, 30])
handles = [Line2D([], [], color=KIM, marker="o", ms=6.5, mfc=KIM, lw=1.3, ls="-"),
           Line2D([], [], color=PINK, marker="o", ms=5.5, mfc=PINK, lw=1.1, ls=":"),
           Line2D([], [], color=PINK, marker="o", ms=6.5, mfc="white", mew=1.4, ls="none")]
figa.legend(handles, ["Kim et al. 2025 (digitised)", "ours, L10 (training run)", "ours, L10 (second run, 25 rpm)"],
            loc="upper center", ncol=3, frameon=False, fontsize=13, bbox_to_anchor=(0.5, 1.01), handlelength=2.2, columnspacing=1.2)
figa.subplots_adjust(left=0.135, right=0.985, top=0.66, bottom=0.2, wspace=0.42)
figa.savefig(OUT / "s7a_kim.png", dpi=200)
print("wrote", OUT / "s7a_kim.png")

# ------------------------------------------------------------------ (b)
import scripts.test_h1_mfbml as H1  # noqa: E402
import scripts.test_h2_mfbml as H2  # noqa: E402
lf_pts, _ = T.load_points([(r, t) for r in T.RPMS for t in T.LF_THETAS], T.lf_run)
l9_pts, _ = T.load_points(T.l9_points(), T.l9_run)
hf_pts, _ = T.load_points(T.hf_points(), T.hf_run)
test = {r: (T.qois(f"x8_l10_rpm{r:g}_th5.5"), T.qois(f"x8_l9_rpm{r:g}_th5.5")) for r in T.TEST_RPM}
hk = sorted(hf_pts, key=lambda k: (k[1], k[0]))
X = np.array(hk)
Xt = np.array([[r, T.TH] for r in sorted(test)])
test_pts = {}
for q in ("tau_95_t", "edr_max_t"):
    k8 = sorted(lf_pts, key=lambda k: (k[1], k[0]))
    lf = T.train_lf2d(np.array(k8), np.log([lf_pts[k][q] for k in k8]))
    d9 = {**{k: v[q] for k, v in l9_pts.items()}, **{(r, T.TH): test[r][1][q] for r in test}}
    k9 = sorted(d9, key=lambda k: (k[1], k[0]))
    m1 = T.FullBayesKRRLRGP(T.BOUNDS, lf_model=lf, basis="linear").train(np.array(k9), np.log([d9[k] for k in k9]))
    ad = T.Offset(m1)
    yh = np.log([hf_pts[k][q] for k in hk])
    m2 = T.FullBayesKRRLRGP(T.BOUNDS, lf_model=ad, basis="linear").train(X, yh)
    mean = np.ravel(m2.predict_mean(Xt))
    qq = m2.quantile(Xt, [0.025, 0.975])
    ad.off = T.DELTA
    rho = (np.ravel(m2.predict_mean(Xt)) - mean) / T.DELTA
    ad.off = 0.0
    hw = np.sqrt(((qq[:, 1] - qq[:, 0]) / 2) ** 2 + (1.96 * rho * T.sd_log(m1, Xt)) ** 2)
    obs = np.array([test[r][0][q] for r in sorted(test)])
    test_pts[q] = dict(obs=obs, pred=np.exp(mean), lo=np.exp(mean - hw), hi=np.exp(mean + hw))
    ref = json.load(open(REPO / "experiments/multifidelity/test_h4_testset_kim.json"))[q]["A"]["pred"]
    print(f"(b) {q}: max |recomputed - JSON pred| = {np.max(np.abs(np.exp(mean) - np.array(ref))):.2e}")

figb, axb = plt.subplots(1, 2, figsize=(9.6, 2.55), dpi=200)
panels = [("leave-one-out, 13 L10 runs", "loo"), ("test set, theta 5.5, 6 new runs", "test")]
for ax, (ttl, which) in zip(axb, panels):
    for q, mk, filled in (("tau_95_t", "o", False), ("edr_max_t", "s", True)):
        if which == "loo":
            pts = hd[q]["points_B"]
            obs = np.exp([p["obs"] for p in pts])
            pred = np.exp([p["mean"] for p in pts])
            lo = np.exp([p["lo"] for p in pts])
            hi = np.exp([p["hi"] for p in pts])
        else:
            tp = test_pts[q]
            obs, pred, lo, hi = tp["obs"], tp["pred"], tp["lo"], tp["hi"]
        ax.errorbar(obs, pred, yerr=[pred - lo, hi - pred], fmt=mk, ms=6.5, mfc=(PINK if filled else "white"), mec=PINK,
                    color=PINK, ecolor=PINK, elinewidth=1.0, capsize=2.5, ls="none")
    lo_, hi_ = ax.get_xlim()[0] * 0.8, ax.get_xlim()[1] * 1.25
    ax.plot([lo_, hi_], [lo_, hi_], color="0.7", lw=1.0, zorder=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(ttl, fontsize=14, pad=3)
    ax.set_xlabel("observed", fontsize=14)
    ax.set_ylabel("predicted", fontsize=14)
    ax.grid(**GRID)
    ax.tick_params(labelsize=13)
handles = [Line2D([], [], color=PINK, marker="o", mfc="white", mec=PINK, ms=6.5, ls="none"),
           Line2D([], [], color=PINK, marker="s", mfc=PINK, mec=PINK, ms=6.5, ls="none"),
           Line2D([], [], color=PINK, lw=1.0, ls="-")]
figb.legend(handles, [r"$\tau_{95}$ (Pa)", r"$\varepsilon_{\max}$ (EDR, W m$^{-3}$)", "95 % band"],
            loc="upper center", ncol=3, frameon=False, fontsize=13, bbox_to_anchor=(0.5, 1.01), handlelength=1.6, columnspacing=1.2)
figb.subplots_adjust(left=0.1, right=0.985, top=0.66, bottom=0.2, wspace=0.42)
figb.savefig(OUT / "s7b_predict.png", dpi=200)
print("wrote", OUT / "s7b_predict.png")
