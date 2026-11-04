"""Slide 7 ("result 1"): (a) validation against Kim at theta 7; (b) held-out prediction with 95 % bands.

Data read directly (mtime at writing, `stat -c '%y %n'`):
  /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/experiments/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv
      (2026-07-21 11:11; digitised Kim 2024 Fig 13, theta 7. Columns tau_liq_mean -> our tau_amp_t, Ediss_liq_mean -> our edr_max_t;
       the mapping and the 2026-10-22 comparison are in diary.md, entry 2026-10-22 "CORRECTION")
  .../experiments/multifidelity/h3_band_check_kim.json  (2026-10-23 09:34; LOO points, variant B = L9 variance propagated)
  .../experiments/multifidelity/test_h4_testset_kim.json (2026-10-23 09:34; used only for the verdict, not the bands)
  .../runs/<x8_l10_rpm*_th7 / x8_l10_rpm*_th5.5>/shear_stress.dat, params.json, results (through scripts/hydro_dataset_v2.qois)
  .../experiments/multifidelity/test_h1_mfbml, test_h2_mfbml, h3_band_check, test_h4_testset (code, same recipe as the JSON)
Per-point 95 % bands of the test set are recomputed with the same code path (variant A) because the JSON keeps only the mean
and the mean half-width. The recomputed means are checked against test_h4_testset_kim.json and must agree.

Outputs (figs/):  s7a_kim.png  (tau_amp_t and edr_max_t vs rpm, theta 7)
                  s7b_predict.png  (LOO and theta-5.5 test set, predicted vs observed, 95 % bands)
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
sys.argv = [sys.argv[0]]              # test_h4_testset parses sys.argv at import
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
warnings.filterwarnings("ignore")
import scripts.test_h4_testset as T   # noqa: E402
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from scripts.figstyle import level_colour, KIM, MK, series_kw, rcparams, GUIDE  # noqa: E402

OUT = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026/figs")
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update(rcparams())
plt.rcParams.update({"font.size": 14})
PINK = level_colour(10)
# quantity -> (figstyle marker, statistic for fill). tau_amp_t and edr_max_t are time extrema -> filled.
# tau_95 is a percentile of the spatial distribution (not an extremum) -> hollow.
MARK = {"tau": MK["tau"], "edr": MK["ediss"]}

# ----------------------------------------------------------------- (a) Kim validation
kim = pd.read_csv(REPO / "experiments/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv", skiprows=[1])
kim = kim[(kim["RPM"] >= 17.5) & (kim["RPM"] <= 32.5)].sort_values("RPM")
rpms_kim = kim["RPM"].to_numpy()
tau_kim = kim["tau_liq_mean"].to_numpy()
edr_kim = kim["Ediss_liq_mean"].to_numpy()

hd = json.load(open(REPO / "experiments/multifidelity/h3_band_check_kim.json"))
loo = {q: hd[q]["points_B"] for q in ("tau_95_t", "edr_max_t")}
rpm_th7 = sorted({p["rpm"] for p in loo["tau_95_t"] if p["theta"] == 7.0})
ours = {}
for r in rpm_th7:
    run = T.hf_run(r, 7.0)
    ours[r] = qois(run)
    print(f"(a) run {run}: tau_amp_t={ours[r]['tau_amp_t']:.6g} edr_max_t={ours[r]['edr_max_t']:.6g}")
print("(a) Kim theta 7 rows used:", {r: (round(t, 6), round(e, 6)) for r, t, e in zip(rpms_kim, tau_kim, edr_kim)})
rel_tau, rel_edr = [], []
for r in rpm_th7:
    if r in set(rpms_kim):
        k = int(np.where(rpms_kim == r)[0][0])
        rel_tau.append(ours[r]["tau_amp_t"] / tau_kim[k] - 1)
        rel_edr.append(ours[r]["edr_max_t"] / edr_kim[k] - 1)
print("(a) ours/Kim - 1, tau_amp:", np.round(100 * np.array(rel_tau), 1), "% ; edr_max:", np.round(100 * np.array(rel_edr), 1), "%")

fig, axes = plt.subplots(2, 1, figsize=(4.8, 3.9), dpi=200, sharex=True)
for ax, (key, kv, lab) in zip(axes, [("tau_amp_t", tau_kim, r"$\tau_{\mathrm{amp}}$"), ("edr_max_t", edr_kim, r"$\varepsilon_{\max}$")]):
    qk = "tau" if key.startswith("tau") else "edr"
    ax.plot(rpms_kim, kv, **series_kw(KIM, MARK[qk], stat="max", ours=False))
    ours_v = [ours[r][key] for r in rpm_th7]
    ax.plot(rpm_th7, ours_v, **series_kw(PINK, MARK[qk], stat="max", ours=True))
    ax.set_yscale("log")
    ax.set_ylabel(lab, fontsize=14)
    ax.grid(**{"which": "major", "ls": ":", "alpha": 0.4})
axes[1].set_xlabel("rpm", fontsize=14)
axes[0].tick_params(labelsize=13)
axes[1].tick_params(labelsize=13)
handles = [Line2D([], [], **series_kw(KIM, MARK["tau"], stat="max")),
           Line2D([], [], **series_kw(PINK, MARK["tau"], stat="max", ours=True))]
fig.legend(handles, ["Kim et al. (digitised)", "ours, L10"], loc="lower center", ncol=2, frameon=False, fontsize=14, bbox_to_anchor=(0.5, -0.04))
fig.subplots_adjust(left=0.22, right=0.97, top=0.97, bottom=0.25, hspace=0.12)
fig.savefig(OUT / "s7a_kim.png", dpi=200)
print("wrote", OUT / "s7a_kim.png")

# ----------------------------------------------------------------- (b) predictions with 95 % bands
lf_pts, _ = T.load_points([(r, t) for r in T.RPMS for t in T.LF_THETAS], T.lf_run)
l9_pts, _ = T.load_points(T.l9_points(), T.l9_run)
hf_pts, _ = T.load_points(T.hf_points(), T.hf_run)
test = {}
for r in T.TEST_RPM:
    h = T.qois(f"x8_l10_rpm{r:g}_th{T.TH:g}")
    q9 = T.qois(f"x8_l9_rpm{r:g}_th{T.TH:g}")
    test[r] = (h, q9)
print("(b) test runs:", sorted(test))
hk = sorted(hf_pts, key=lambda k: (k[1], k[0]))
X = np.array(hk)
Xt = np.array([[r, T.TH] for r in sorted(test)])
test_pts = {}
for q in ("tau_95_t", "edr_max_t"):
    k8 = sorted(lf_pts, key=lambda k: (k[1], k[0]))
    lf = T.train_lf2d(np.array(k8), np.log([lf_pts[k][q] for k in k8]))
    l9_new = {(r, T.TH): test[r][1][q] for r in test}
    d9 = {**{k: v[q] for k, v in l9_pts.items()}, **l9_new}
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
    test_pts[q] = dict(rpm=sorted(test), obs=obs, pred=np.exp(mean), lo=np.exp(mean - hw), hi=np.exp(mean + hw))
    ref = json.load(open(REPO / "experiments/multifidelity/test_h4_testset_kim.json"))[q]["A"]["pred"]
    print(f"(b) {q}: recomputed pred vs JSON pred max |diff| = {np.max(np.abs(np.exp(mean) - np.array(ref))):.2e}")
    print(f"(b) {q}: test obs {np.round(obs, 5)}")

fig, axes = plt.subplots(1, 2, figsize=(4.8, 2.9), dpi=200)
axes = np.atleast_1d(axes)
panels = [("LOO, 13 runs", "loo"), ("test, theta 5.5", "test")]
for ax, (title, which) in zip(axes, panels):
    for q, qk in (("tau_95_t", "tau"), ("edr_max_t", "edr")):
        if which == "loo":
            pts = loo[q]
            obs = np.exp([p["obs"] for p in pts])
            pred = np.exp([p["mean"] for p in pts])
            lo = np.exp([p["lo"] for p in pts])
            hi = np.exp([p["hi"] for p in pts])
        else:
            tp = test_pts[q]
            obs, pred, lo, hi = tp["obs"], tp["pred"], tp["lo"], tp["hi"]
        ec = PINK
        ax.errorbar(obs, pred, yerr=[pred - lo, hi - pred], fmt=MARK[qk], ms=6, mfc=(ec if qk == "edr" else "w"), mec=ec, color=ec, ecolor=ec, elinewidth=0.9, capsize=2.5, ls="none")
    allv = np.concatenate([np.array(ax.get_xlim()), np.array(ax.get_ylim())])
    lo_, hi_ = 0.9 * min(ax.get_xlim()[0], ax.get_ylim()[0]), 1.1 * max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([lo_, hi_], [lo_, hi_], color=GUIDE, lw=1.0, zorder=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(title, fontsize=14, pad=4)
    ax.set_xlabel("observed", fontsize=14)
    ax.tick_params(labelsize=12)
axes[0].set_ylabel("predicted", fontsize=14)
handles = [Line2D([], [], **{**series_kw(PINK, MARK["tau"]), "ls": "none"}), Line2D([], [], **{**series_kw(PINK, MARK["edr"], stat="max"), "ls": "none"}), Line2D([], [], color=PINK, lw=0.9, marker="|", ms=8, ls="none")]
fig.legend(handles, [r"$\tau_{95}$", r"$\varepsilon_{\max}$", "bars: 95 % band"], loc="lower center", ncol=3, frameon=False, fontsize=13, bbox_to_anchor=(0.5, -0.04), handletextpad=0.4, columnspacing=1.0)
fig.subplots_adjust(left=0.16, right=0.98, top=0.9, bottom=0.33, wspace=0.45)
fig.savefig(OUT / "s7b_predict.png", dpi=200)
print("wrote", OUT / "s7b_predict.png")
