"""Is the LF (L8) model the bottleneck in 2D? (2026-11-10)
Hypothesis: MF gains are small in 2D because L8 is sparse (80 rows, campaign-clustered): at held-out L9 sites with no
L8 run, the LF value is an interpolation with its own error. Falsified if B2-lin error at held-out sites WITH an L8 run
(paired) is no lower, relative to L9-only, than at sites WITHOUT one.
Same protocol as mf_ablation.py (20 site draws, n9 in 4, 8, 16). Variants: B1 (L9 only), B2-lin.
Also printed: LOO RMSE of the L8 KRR at the paired sites (LF interpolation error), vs the L9 sd.
Usage: uv run python scripts/bo2d/diag_lf_bottleneck.py"""
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts")); sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))
warnings.filterwarnings("ignore")
from mfbml.methods.kernel_ridge_regression import KernelRidgeRegression as KRR  # noqa: E402
from scripts.bo2d import bo  # noqa: E402
from scripts.bo2d import loop as LP  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402
import campaign as C  # noqa: E402
D = C.DOMAIN


def krr(X, y):
    m = KRR(design_space=D, params_optimize=True, noise_data=True, optimizer_restart=10, seed=42)
    m.train(X, y.reshape(-1, 1), portion_test=0.2)
    return m


def fb(Xh, yh, lf, floor):
    m = FullBayesKRRLRGP(D, lf_model=lf, basis="linear" if lf is not None else "ordinary",
                         grid=C.hf_grid(floor, float(np.std(yh))))
    return m.train(Xh, yh)


for q in ("kLa_1T_10", "dtmix_0.25", "tau_95"):
    d = pd.read_csv(ROOT / "scripts/bo2d/state" / f"truth_snapshot_{q}.csv")
    d = d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]; d = LP.as_objective(d, q)
    d["ly"] = np.log(d.y.astype(float)); d["key"] = d.rpm.round(2).astype(str) + "_" + d.theta.round(2).astype(str)
    floor = C.noise_floor(q)
    xy = lambda s: (s[["rpm", "theta"]].to_numpy(float), s.ly.to_numpy(float))
    d8 = d[d.level == 8]; d9 = d[d.level == 9]; X8, y8 = xy(d8); lf8 = krr(X8, y8)
    paired = set(d8.key) & set(d9.key); sites = d9.key.unique()
    # LF interpolation error at paired sites: drop that site's L8 rows, predict it
    e = []
    for k in list(paired):
        m = d8.key != k; lf = krr(*xy(d8[m])); xs, ys = xy(d8[~m])
        e.append(np.mean(np.ravel(lf.predict(xs))) - ys.mean())
    print(f"{q}: L8 KRR leave-site-out RMSE at the {len(paired)} paired sites = {np.sqrt(np.mean(np.square(e))):.3f} ln; "
          f"L9 sd = {d9.ly.std():.3f}; L8 sd = {d8.ly.std():.3f}", flush=True)
    out = []
    for n9 in (4, 8, 16):
        for rep in range(20):
            rng = np.random.default_rng(100 * rep + n9)
            tr = rng.choice(sites, n9, replace=False); A = d9[d9.key.isin(tr)]; B = d9[~d9.key.isin(tr)]
            Xa, ya = xy(A); Xb, yb = xy(B); pb = B.key.isin(paired).to_numpy()
            for name, lf in (("B1", None), ("B2-lin", lf8)):
                mu, _ = bo.gauss(fb(Xa, ya, lf, floor), Xb); se = (mu - yb) ** 2
                for grp, msk in (("paired", pb), ("unpaired", ~pb)):
                    if msk.any():
                        out.append(dict(n9=n9, rep=rep, v=name, grp=grp, rmse=np.sqrt(se[msk].mean())))
    o = pd.DataFrame(out).groupby(["grp", "v", "n9"]).rmse.median().unstack("n9")
    rel = o.xs("B2-lin", level="v") / o.xs("B1", level="v")
    print(o.round(3).to_string()); print("B2-lin / B1:\n" + rel.round(2).to_string() + "\n", flush=True)
