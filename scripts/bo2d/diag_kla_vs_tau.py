"""Hypothesis (2026-11-11 08:00, after the user rejected 'flat kLa'): kLa and tau_95 are both set by dissipation, so at
L9 kLa ~ f(tau). Then under a tau cap the feasible kLa maximum is ~f(cap) along the whole cap contour (a ridge of
near-ties), and 5 % regret is reached at many sites. Test on the frozen truth snapshots (L9, same run ids):
Spearman (kLa, tau) over all L9 sites; LOO error of ln kLa ~ a + b ln tau (+ c theta) vs the mean; and the kLa
spread among sites whose tau lies within +-15 % of the cap used by the benchmark (0.014 Pa).
  uv run python scripts/bo2d/diag_kla_vs_tau.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

S = Path(__file__).parent / "state"
k = pd.read_csv(S / "truth_snapshot_kLa_1T_10.csv"); t = pd.read_csv(S / "truth_snapshot_tau_95.csv")
k = k[k.level == 9][["run", "rpm", "theta", "y"]].rename(columns={"y": "kla"})
t = t[t.level == 9][["run", "y"]].rename(columns={"y": "tau"})
d = k.merge(t, on="run").dropna()
d = d[(d.kla > 0) & (d.tau > 0)]
lk, lt, th = np.log(d.kla.to_numpy()), np.log(d.tau.to_numpy()), d.theta.to_numpy()
print(f"n={len(d)}  Spearman(kLa, tau) {spearmanr(lk, lt).correlation:+.2f}   (rpm {spearmanr(lk, d.rpm).correlation:+.2f}, "
      f"theta {spearmanr(lk, th).correlation:+.2f})")


def loo(A, y):
    e = []
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None); e.append(y[i] - A[i] @ c)
    return float(np.sqrt(np.mean(np.square(e))))


one = np.ones_like(lk)
r0 = loo(one[:, None], lk)
print(f"LOO ln kLa: mean {r0:.3f} | a + b ln tau {loo(np.column_stack([one, lt]), lk):.3f} | "
      f"+ c theta {loo(np.column_stack([one, lt, th]), lk):.3f} | rpm + theta (no tau) {loo(np.column_stack([one, d.rpm, th]), lk):.3f}")
cap = 0.014
near = np.abs(lt - np.log(cap)) < np.log(1.15)
print(f"sites with tau within 15 % of the cap: n={near.sum()}  kLa {np.round(d.kla.to_numpy()[near], 1)}  "
      f"at {np.round(d[['rpm', 'theta']].to_numpy()[near], 1).tolist()}")
print(d.assign(tau_mPa=1000 * d.tau).sort_values("tau")[["rpm", "theta", "kla", "tau_mPa"]].round(2).to_string(index=False))
