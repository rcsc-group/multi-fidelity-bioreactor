"""Is the poor pooled L8 -> L9 kLa transfer a PROTOCOL artefact? (2026-11-10, user: "the 1D figures look like a good
transfer - this smells"). The theta = 7 L9 points are fig9 runs (cold start, older binary); the off-line L9 points are
c2d warm starts (xlevel binary). Hypothesis: within each set the map is good; pooled, a protocol offset breaks it.
Falsified if R(within theta=7) and R(off-line) are as poor as the pooled R, and a protocol indicator does not help.
LOO R = mean LOO residual^2 / var(y9), maps on ln values. Usage: uv run python scripts/bo2d/diag_kla_protocol.py"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import loop as LP  # noqa: E402


def loo(A, y):
    e = []
    for i in range(len(y)):
        m = np.arange(len(y)) != i; c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None); e.append(y[i] - A[i] @ c)
    return np.mean(np.square(e)) / y.var()


for q in ("kLa_1T_10", "dtmix_0.25", "tau_95"):
    d = pd.read_csv(Path("scripts/bo2d/state") / f"truth_snapshot_{q}.csv")
    d = d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]; d = LP.as_objective(d, q)
    d["ly"] = np.log(d.y.astype(float)); d["key"] = d.rpm.round(2).astype(str) + "_" + d.theta.round(2).astype(str)
    p = d[d.level == 8].groupby("key").agg(y8=("ly", "mean"), th=("theta", "first"), rpm=("rpm", "first")).join(
        d[d.level == 9].groupby("key").agg(y9=("ly", "mean"), run9=("run", "first")), how="inner")
    p["fig9"] = p.run9.str.startswith("fig9").astype(float)
    one = np.ones(len(p)); f, y = p.y8.to_numpy(), p.y9.to_numpy()
    on = p.fig9.to_numpy() == 1
    res = {"pooled a+b f8": loo(np.column_stack([one, f]), y),
           "pooled + protocol offset": loo(np.column_stack([one, f, p.fig9]), y),
           "pooled + protocol offset + theta": loo(np.column_stack([one, f, p.fig9, p.th]), y)}
    if on.sum() > 3:
        res[f"theta=7 only (n={on.sum()})"] = loo(np.column_stack([np.ones(on.sum()), f[on]]), y[on])
    if (~on).sum() > 3:
        res[f"off-line only (n={(~on).sum()})"] = loo(np.column_stack([np.ones((~on).sum()), f[~on]]), y[~on])
    print(q + ": " + " | ".join(f"{k}: {v:.2f}" for k, v in res.items()))
    if q == "kLa_1T_10":
        for g, m in (("theta=7 (fig9)", on), ("off-line (c2d)", ~on)):
            c = np.polyfit(f[m], y[m], 1)
            print(f"   {g}: y9 = {c[1]:+.2f} + {c[0]:.2f} y8 ; mean(y9 - y8) = {np.mean(y[m] - f[m]):+.2f} ln")
