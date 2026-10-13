"""H1 evaluation (pre-registered, diary.md): mfbml (T2b full-Bayes residual GP) on (rpm, theta).

Rule (quoted): "mfbml as in T2b (full-Bayes residual GP) on input (rpm, theta); LF = the L8 grid without
theta 2 (40 points); HF = these 6 + the 7 theta-7 L10 points at 17.5-32.5 rpm (xlevel set) = 13.
Leave-one-out over the 13 HF points. PASS if, for tau_95 AND EDR: LOO relative RMSE <= 10%, AND mfbml
RMSE below both baselines (L8 x global ratio; HF-only GP, integrated), AND the 95% intervals cover
>= 80% of the 13 (>= 11). tau_mean reported."
PENDING while fewer than 13 HF runs have runs/<id>/results.json.

Usage: uv run python scripts/test_h1_mfbml.py
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

OUT = ROOT / "experiments" / "multifidelity" / "test_h1_mfbml.json"
BOUNDS = np.array([[15.0, 37.5], [4.0, 9.0]])
RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
LF_THETAS = (4, 5.5, 7, 9)
QOIS = ("tau_95_t", "edr_mean_t", "tau_mean_t")
CRITERIA_QOIS = ("tau_95_t", "edr_mean_t")
N_REQUIRED = 13
RMSE_MAX = 0.10
COVER_MIN = 11   # of 13 (>= 80%)


def lf_run(rpm, theta):
    return f"kmix_l8_rpm{rpm:g}" if theta == 7 else f"t6_l8_rpm{rpm:g}_th{theta:g}"


def hf_run(rpm, theta):
    if theta == 7:
        return f"l10_fig13a_xlevel_rpm{rpm:g}"
    return f"x8_l10_rpm{rpm:g}_th{theta:g}"


def hf_points():
    pts = [(r, 7.0) for r in (17.5, 20, 22.5, 25, 27.5, 30, 32.5)]
    pts += [(r, float(t)) for r in (17.5, 25, 32.5) for t in (4, 9)]
    return sorted(pts, key=lambda p: (p[1], p[0]))


# ---------------------------------------------------------------- pure logic
def ratio_baseline(yh_train, lf_train, lf_test):
    """L8 x global ratio: LF value at the held-out point times mean(HF/LF) over the training points."""
    return float(lf_test) * float(np.mean(np.asarray(yh_train, float) / np.asarray(lf_train, float)))


def rel_rmse_ratio(pred, obs):
    pred, obs = np.asarray(pred, float), np.asarray(obs, float)
    return float(np.sqrt(np.mean((pred / obs - 1.0) ** 2)))


def h1_verdict(stats: dict, n_required: int = N_REQUIRED) -> dict:
    """stats[qoi] = dict(mf_rmse, hf_rmse, ratio_rmse, mf_cover (hits), n) for tau_95_t and edr_mean_t."""
    n = min(stats[q]["n"] for q in CRITERIA_QOIS)
    if n < n_required:
        return dict(verdict=f"PENDING ({n}/{n_required} HF points)", failed=[])
    cover_min = int(np.ceil(0.8 * n_required - 1e-9))   # 11 of 13
    failed = []
    for q in CRITERIA_QOIS:
        s = stats[q]
        if not s["mf_rmse"] <= RMSE_MAX:
            failed.append((q, "rmse<=10%"))
        if not s["mf_rmse"] < s["hf_rmse"]:
            failed.append((q, "beats HF-only"))
        if not s["mf_rmse"] < s["ratio_rmse"]:
            failed.append((q, f"beats ratio"))
        if not s["mf_cover"] >= cover_min:
            failed.append((q, f"coverage>={cover_min}/{n_required}"))
    return dict(verdict="FAIL" if failed else "PASS", failed=failed)


# ---------------------------------------------------------------- data
def has_results(run):
    return (ROOT / "runs" / run / "results.json").exists()


def load_points(points, run_fn):
    """-> dict (rpm, theta) -> qois dict, only for runs with results.json; plus list of missing ids."""
    from scripts.hydro_dataset_v2 import qois
    got, missing = {}, []
    for rpm, th in points:
        run = run_fn(rpm, th)
        if has_results(run):
            got[(float(rpm), float(th))] = qois(run)
        else:
            missing.append(run)
    return got, missing


def train_lf2d(X, y):
    """Same LF KRR as T2b (train_lf: LooKRR seed 42, restart 10, portion_test 0.2), 2-D inputs."""
    from scripts.test_t3_learning_curve import CachedLooKRR
    lf = CachedLooKRR(design_space=BOUNDS, params_optimize=True, noise_data=True,
                      optimizer_restart=10, seed=42)
    lf.train(np.asarray(X, float), np.asarray(y, float), portion_test=0.2)
    return lf


def loo(q, lf_pts, hf_pts):
    from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP
    keys = sorted(lf_pts, key=lambda k: (k[1], k[0]))
    XL = np.array(keys)
    lf = train_lf2d(XL, [lf_pts[k][q] for k in keys])
    hk = sorted(hf_pts, key=lambda k: (k[1], k[0]))
    X = np.array(hk)
    yh = np.array([hf_pts[k][q] for k in hk])
    ylf = np.array([lf_pts[k][q] for k in hk])        # L8 run at the same (rpm, theta)
    rows = []
    for i, k in enumerate(hk):
        tr = np.arange(len(hk)) != i
        xt = X[i:i + 1]
        r = dict(rpm=k[0], theta=k[1], obs=float(yh[i]))
        for name, basis in (("mf", "linear"), ("hf", "ordinary")):
            m = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis=basis).train(X[tr], yh[tr])
            r[f"{name}_pred"] = float(m.predict_mean(xt)[0])
            qq = m.quantile(xt, [0.025, 0.975])[0]
            r[f"{name}_lo"], r[f"{name}_hi"] = float(qq[0]), float(qq[1])
        r["ratio_pred"] = ratio_baseline(yh[tr], ylf[tr], ylf[i])
        rows.append(r)
    obs = np.array([r["obs"] for r in rows])
    hit = lambda n: int(sum(r[f"{n}_lo"] <= r["obs"] <= r[f"{n}_hi"] for r in rows))  # noqa: E731
    s = dict(n=len(rows), mf_rmse=rel_rmse_ratio([r["mf_pred"] for r in rows], obs),
             hf_rmse=rel_rmse_ratio([r["hf_pred"] for r in rows], obs),
             ratio_rmse=rel_rmse_ratio([r["ratio_pred"] for r in rows], obs),
             mf_cover=hit("mf"), hf_cover=hit("hf"))
    return s, rows


def main():
    warnings.filterwarnings("ignore")
    lf_pts, lf_missing = load_points([(r, t) for r in RPMS for t in LF_THETAS], lf_run)
    hf_all = hf_points()
    hf_pts, hf_missing = load_points(hf_all, hf_run)
    n = len(hf_pts)
    print(f"LF points {len(lf_pts)}/40 (missing: {lf_missing}); HF points {n}/{N_REQUIRED}; "
          f"HF missing: {hf_missing}")
    res, stats = {}, {}
    # every HF point needs its L8 counterpart in the LF set
    hf_use = {k: v for k, v in hf_pts.items() if k in lf_pts}
    if len(hf_use) < n:
        print("WARNING: HF points without L8 counterpart dropped:", sorted(set(hf_pts) - set(hf_use)))
    if len(hf_use) >= 3 and len(lf_pts) >= 8:
        for q in QOIS:
            s, rows = loo(q, lf_pts, hf_use)
            stats[q] = s
            res[q] = dict(stats=s, points=rows)
            print(f"\n{q}: n={s['n']}  rel RMSE  MF {s['mf_rmse']:.1%}  HF-only {s['hf_rmse']:.1%}  "
                  f"L8xratio {s['ratio_rmse']:.1%}  coverage MF {s['mf_cover']}/{s['n']}  "
                  f"HF-only {s['hf_cover']}/{s['n']}")
            print(f"{'rpm':>6}{'theta':>6}{'obs':>11}{'MF':>11}{'MF 95% interval':>26}{'HF-only':>11}{'ratio':>11}")
            for r in rows:
                print(f"{r['rpm']:6g}{r['theta']:6g}{r['obs']:11.4g}{r['mf_pred']:11.4g}"
                      f"{'[%.4g, %.4g]' % (r['mf_lo'], r['mf_hi']):>26}{r['hf_pred']:11.4g}{r['ratio_pred']:11.4g}")
        verdict = h1_verdict(stats)
    else:
        verdict = dict(verdict=f"PENDING ({len(hf_use)}/{N_REQUIRED} HF points)", failed=[])
    print(f"\nVERDICT: {verdict['verdict']}", verdict["failed"] or "")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(dict(verdict=verdict, n_hf=len(hf_use), hf_missing=hf_missing, lf_missing=lf_missing,
                   results=res), open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
