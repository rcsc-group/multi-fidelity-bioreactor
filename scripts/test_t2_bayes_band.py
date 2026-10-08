"""T2 (pre-registered, diary 2026-10-08): Bayesian residual amplitude in KRR-LR-GPR, re-run of test A.

Pre-registration (diary.md, quoted):
  "T2 (zero compute): prior on the residual amplitude in KRR-LR-GPR (fully Bayesian over the GP
   amplitude, replacing the plug-in sigma2 that collapses with 1 residual dof). PASS = pooled 95%
   coverage >= 85% on test A (tau, all 4 cases, n_HF = 3) AND median rel RMSE <= 1.2x the test-A
   value in every case."
Task spec: transfer coefficients flat prior; GP amplitude IG(a0=2, b0=0.1); predictive Student-t,
nu = n - p + 2 a0, location = reference universal-kriging mean. Derivation in
scripts/mfbml_local/krr_lr_gpr_bayes.py. Tests: tests/test_krr_lr_gpr_bayes.py.

Protocol = test A (scripts/test_mf_l10_holdout.py): tau_mean_t, tau_95_t; LF L8, L9; n_HF = 3
(endpoints 17.5, 37.5 + one interior rpm, 7 designs); held-out L10 rpm scored. Data, LooKRR, fit(),
metric functions are imported from test A.

Choices where the spec was silent (closest to the spec):
  * Coverage band = location +- t_{0.975, nu} * scale, scale includes the noise term nr^2 inside
    the amplitude-scaled variance (eq. 3 of the module), i.e. total predictive for a new observation.
    The reference band instead adds the noise unscaled; this is a consequence of putting noise in K.
  * "test-A value" of the median rel RMSE = recomputed MF value in this run (also compared with the
    stored test_mf_l10_holdout.json; identical locations are asserted to 1e-9).
  * Pooled coverage = hits / held-out points over all designs of all 4 cases.
  * Failed fits (Cholesky) are counted and skipped.
Reuse helper for T3/T4: as_bayes(model) re-labels a trained reference model (same trained state; the
T2 run asserts that an independently trained subclass gives identical hyperparameters and mean).

Usage: uv run python scripts/test_t2_bayes_band.py
"""
from __future__ import annotations

import itertools
import json
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import figstyle as fs                                              # noqa: E402
from scripts.plot_mf_fig13 import LooKRR                                        # noqa: E402
from scripts.test_mf_l10_holdout import (load, col, fit, pred, rel_rmse, ENDS, DS,  # noqa: E402
                                          QOIS, LFS, OUT)
from mfbml_local.krr_lr_gpr_bayes import KernelRidgeLinearGaussianProcessBayes as Bayes  # noqa: E402

A0, B0 = 2.0, 0.1


def fit_bayes(xh, yh, xl, yl):
    """Independently trained subclass (same settings as test A's fit(), basis linear)."""
    lf = LooKRR(design_space=DS, params_optimize=True, noise_data=True,
                optimizer_restart=10, seed=42)
    m = Bayes(design_space=DS, optimizer_restart=10, lf_poly_order="linear", seed=42,
              lf_model=lf, noise_prior=None, a0=A0, b0=B0)
    m.train(X=[col(xh), col(xl)], Y=[np.asarray(yh, float), np.asarray(yl, float)])
    return m


def as_bayes(m, a0=A0, b0=B0):
    """View a trained reference model as the Bayes subclass (shares trained state)."""
    m.__class__ = Bayes
    m.a0, m.b0 = a0, b0
    return m


def main():
    t0 = time.time()
    warnings.filterwarnings("ignore")
    D = load()
    old = json.load(open(OUT / "test_mf_l10_holdout.json"))["results"]
    hrpm = np.array(sorted(D[10]))
    interior = [r for r in hrpm if r not in ENDS]
    res, hits_tot, n_tot, models = {}, 0, 0, {}
    ref_hits_tot = 0
    for q in QOIS:
        yh_all = np.array([D[10][r][q] for r in hrpm])
        for lfl in LFS:
            xl = np.array(sorted(D[lfl]))
            yl = np.array([D[lfl][r][q] for r in xl])
            rr, rb, hits, ref_hits, nt, hw, fails, maxdiff = [], [], 0, 0, 0, [], 0, 0.0
            zs, krn, nrs = [], [], []
            for mid in itertools.combinations(interior, 1):
                tr = list(ENDS) + list(mid)
                msk = np.isin(hrpm, tr)
                xh, yh, xt, yt = hrpm[msk], yh_all[msk], hrpm[~msk], yh_all[~msk]
                try:
                    mf = fit("linear", xh, yh, xl, yl)
                    bm = fit_bayes(xh, yh, xl, yl)
                except Exception as e:  # noqa: BLE001
                    fails += 1
                    print("  FAIL", q, lfl, tr, e)
                    continue
                mum, sdm = pred(mf, xt)
                mub, scale, nu = bm.predict_dist(col(xt))
                maxdiff = max(maxdiff, float(np.max(np.abs(mum - mub))))
                assert np.allclose(mum, mub, atol=1e-9), "location differs from the reference"
                zs += list(np.abs(mub - yt) / scale)
                krn.append(float(np.ravel(bm.kernel.param)[0])); nrs.append(float(bm.noise / bm.yh_std))
                rr.append(rel_rmse(mum, yt))
                rb.append(rel_rmse(mub, yt))
                half = scale * stats.t.ppf(0.975, nu)
                hits += int(np.sum(np.abs(mub - yt) <= half))
                ref_hits += int(np.sum(np.abs(mum - yt) <= 1.96 * sdm))
                nt += len(yt)
                hw.append(float(np.median(2 * half / mub)))
                models[(q, lfl, tuple(tr))] = (mf, bm)
            hits_tot += hits
            n_tot += nt
            ref_hits_tot += ref_hits
            old_med = old[f"{q}|L{lfl}|n3"]["summary"]["MF"]["median"]
            res[f"{q}|L{lfl}"] = dict(
                med_ref=float(np.median(rr)), med_bayes=float(np.median(rb)), med_testA_stored=old_med,
                ratio=float(np.median(rb) / np.median(rr)), cov_bayes=hits / nt, cov_ref=ref_hits / nt,
                median_rel_band_width=float(np.median(hw)), nu=float(nu), n_fail=fails,
                max_abs_mean_diff=maxdiff, diag_median_abs_err_over_scale=float(np.median(zs)),
                diag_kernel_param_median=float(np.median(krn)), diag_noise_ratio_median=float(np.median(nrs)), n_points=nt)
    pooled, pooled_ref = hits_tot / n_tot, ref_hits_tot / n_tot
    cov_ok = pooled >= 0.85
    rmse_ok = all(v["ratio"] <= 1.2 for v in res.values())
    passed = bool(cov_ok and rmse_ok)

    lines = [f"{'case':18s} {'RMSE ref':>9s} {'RMSE bayes':>11s} {'ratio':>6s} {'cov ref':>8s} {'cov bayes':>10s} {'rel width':>10s}"]
    for k, v in res.items():
        lines.append(f"{k:18s} {v['med_ref']:9.2%} {v['med_bayes']:11.2%} {v['ratio']:6.3f} "
                     f"{v['cov_ref']:8.0%} {v['cov_bayes']:10.0%} {v['median_rel_band_width']:10.1%}")
    lines.append(f"pooled coverage: reference {pooled_ref:.1%}, Bayes {pooled:.1%} ({hits_tot}/{n_tot}); nu = {nu:.0f}")
    table = "\n".join(lines)
    print(table)
    for k, v in res.items():
        print(f"  diag {k}: median |err|/scale = {v['diag_median_abs_err_over_scale']:.3g}, "
              f"kernel param (10^theta) = {v['diag_kernel_param_median']:.3g}, noise ratio = {v['diag_noise_ratio_median']:.3g}")
    print("PASS" if passed else "FAIL", f"(coverage>=85%: {cov_ok}; rmse ratio<=1.2 in all: {rmse_ok})")

    # figure: representative (median-MF-error) n_HF=3 design per case
    plt.rcParams.update(fs.rcparams())
    f, axes = plt.subplots(2, 2, figsize=(8.0, 6.4), sharex=True)
    xg = np.linspace(15, 37.5, 200)
    for i, q in enumerate(QOIS):
        for j, lfl in enumerate(LFS):
            keys = [k for k in models if k[0] == q and k[1] == lfl]
            errs = []
            for k in keys:
                msk = np.isin(hrpm, k[2])
                yh_all = np.array([D[10][r][q] for r in hrpm])
                errs.append(rel_rmse(pred(models[k][0], hrpm[~msk])[0], yh_all[~msk]))
            k = keys[int(np.argsort(errs)[len(errs) // 2])]
            mf, bm = models[k]
            xl = np.array(sorted(D[lfl]))
            yl = np.array([D[lfl][x][q] for x in xl])
            yh = np.array([D[10][x][q] for x in hrpm])
            trm = np.isin(hrpm, k[2])
            mu, scale, nu_ = bm.predict_dist(col(xg))
            _, sd = pred(mf, xg)
            half = scale * stats.t.ppf(0.975, nu_)
            ax = axes[i, j]
            cl, ch, mk = fs.level_colour(lfl), fs.level_colour(10), fs.MK["tau"]
            stat = "mean" if q == "tau_mean_t" else "max"
            ax.fill_between(xg, mu - half, mu + half, color=ch, alpha=0.18, lw=0, label="MF-Bayes 95%")
            ax.plot(xg, mu - 1.96 * sd, color=fs.NEUTRAL, lw=0.9, ls="-.", label="MF 95%")
            ax.plot(xg, mu + 1.96 * sd, color=fs.NEUTRAL, lw=0.9, ls="-.")
            ax.plot(xg, mu, color=ch, lw=1.4, ls="--", label="MF")
            ax.plot(xl, yl, **fs.series_kw(cl, mk, stat=stat, ours=True, ls="none"), label=f"L{lfl}")
            ax.plot(hrpm[trm], yh[trm], **fs.series_kw(ch, mk, stat=stat, ours=True, ls="none"),
                    label="L10 training")
            ax.plot(hrpm[~trm], yh[~trm],
                    **fs.series_kw(ch, mk, stat=stat, ours=True, ls="none", ms=3.0, alpha=0.6),
                    label="L10 held out")
            ax.set_title(f"LF L{lfl}", fontsize=10)
            ax.grid(**fs.GRID_KW)
            if j == 0:
                ax.set_ylabel("mean wall shear stress (Pa)" if i == 0 else "95th-percentile wall shear stress (Pa)")
            if i == 1:
                ax.set_xlabel("rocking speed (rpm)")
            if i == 0 and j == 1:
                ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    f.tight_layout()
    f.savefig(OUT / "test_t2_bayes_band.png", dpi=200, bbox_inches="tight")
    rt = time.time() - t0
    json.dump(dict(pass_=passed, coverage_ok=cov_ok, rmse_ok=rmse_ok, pooled_cov_bayes=pooled,
                   pooled_cov_ref=pooled_ref, a0=A0, b0=B0, cases=res, table=table, runtime_s=rt),
              open(OUT / "test_t2_bayes_band.json", "w"), indent=1)
    print(f"runtime {rt:.1f} s")


if __name__ == "__main__":
    main()
