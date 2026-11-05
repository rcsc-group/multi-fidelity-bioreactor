"""Fitted noise of the L10 stage of the chain (same fit as tools/fig_s7_validation_r2.py, panel b).

Noise sd in ln units, per grid node g: yh_std * sqrt(b_n/a_n)_g * eta_g  (krr_lr_gpr_fullbayes.py:30, :126, :160).
Reports the posterior median and 90 % interval (weights w_g), and the share of posterior mass at the
eta cap (top grid node, eta = 0.3). Writes data/aleatoric_gp_noise_r1.json.
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
sys.argv = [sys.argv[0]]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
warnings.filterwarnings("ignore")
import scripts.test_h4_testset as T  # noqa: E402


def wq(v, w, qs):
    o = np.argsort(v)
    c = np.cumsum(w[o])
    return [float(v[o][np.searchsorted(c, q)]) for q in qs]


lf_pts, _ = T.load_points([(r, t) for r in T.RPMS for t in T.LF_THETAS], T.lf_run)
l9_pts, _ = T.load_points(T.l9_points(), T.l9_run)
hf_pts, _ = T.load_points(T.hf_points(), T.hf_run)
test = {r: (T.qois(f"x8_l10_rpm{r:g}_th5.5"), T.qois(f"x8_l9_rpm{r:g}_th5.5")) for r in T.TEST_RPM}
hk = sorted(hf_pts, key=lambda k: (k[1], k[0]))
out = {}
for q in ("tau_95_t", "edr_max_t", "tau_amp_t"):
    if q not in next(iter(hf_pts.values())):
        print("no", q, "in dataset; skip")
        continue
    k8 = sorted(lf_pts, key=lambda k: (k[1], k[0]))
    lf = T.train_lf2d(np.array(k8), np.log([lf_pts[k][q] for k in k8]))
    d9 = {**{k: v[q] for k, v in l9_pts.items()}, **{(r, T.TH): test[r][1][q] for r in test}}
    k9 = sorted(d9, key=lambda k: (k[1], k[0]))
    m1 = T.FullBayesKRRLRGP(T.BOUNDS, lf_model=lf, basis="linear").train(np.array(k9), np.log([d9[k] for k in k9]))
    m2 = T.FullBayesKRRLRGP(T.BOUNDS, lf_model=T.Offset(m1), basis="linear").train(np.array(hk), np.log([hf_pts[k][q] for k in hk]))
    for name, m in (("L9 stage", m1), ("L10 stage", m2)):
        sd = m.yh_std * np.sqrt(m.b_n / m.a_n) * m.eta
        med, lo, hi = wq(sd, m.w, [0.5, 0.05, 0.95])
        cap = float(m.w[m.eta >= m.eta.max() - 1e-12].sum())
        out[f"{q} {name}"] = dict(sd_ln_median=med, sd_ln_5=lo, sd_ln_95=hi, mass_at_eta_cap=cap, yh_std=m.yh_std, n=m.n)
        print(f"{q:10s} {name}: noise sd (ln) median {med:.4f}  90% [{lo:.4f}, {hi:.4f}]  mass at eta cap {cap:.2f}  (HF std {m.yh_std:.3f}, n {m.n})")
json.dump(out, open(TALK / "data/aleatoric_gp_noise_r1.json", "w"), indent=1)
