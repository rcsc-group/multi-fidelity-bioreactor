"""Generalised to any level pair: uv run python scripts/bo2d/diag_transfer_forms_lv.py LO HI. Is the L8 -> L9 transfer of kLa nonlinear or theta-dependent? (2026-11-10, user: "1D kLa did well -> nonlinear
transfer"). Paired real sites only (no GP). Leave-one-out R = mean LOO residual^2 / var(y9) for maps of increasing
flexibility. If an interaction / quadratic form drops R well below the linear one, the transfer is not a global line.
Usage: uv run python scripts/bo2d/diag_transfer_forms.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
import pandas as pd
from scripts.bo2d import loop as LP

LO, HI = int(sys.argv[1]), int(sys.argv[2])
FORMS = {"a + b f8": lambda f, t: [f],
         "+ c theta": lambda f, t: [f, t],
         "+ d theta*f8": lambda f, t: [f, t, t * f],
         "+ e f8^2": lambda f, t: [f, t, f * f],
         "theta only (no f8)": lambda f, t: [t]}
for q in ("kLa_1T_10", "dtmix_0.25", "tau_95"):
    d = pd.read_csv(Path("scripts/bo2d/state") / f"truth_snapshot_{q}.csv")
    d = d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]; d = LP.as_objective(d, q)
    d["ly"] = np.log(d.y.astype(float)); d["key"] = d.rpm.round(2).astype(str) + "_" + d.theta.round(2).astype(str)
    p = d[d.level == LO].groupby("key").agg(y8=("ly", "mean"), th=("theta", "first"), rpm=("rpm", "first")).join(
        d[d.level == HI].groupby("key").ly.mean().rename("y9"), how="inner")
    f, t, y = p.y8.to_numpy(), p.th.to_numpy(), p.y9.to_numpy()
    out = []
    for name, cols in FORMS.items():
        A = np.column_stack([np.ones_like(f)] + cols(f, t)); e = []
        for i in range(len(y)):
            m = np.arange(len(y)) != i
            c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None); e.append(y[i] - A[i] @ c)
        out.append(f"{name}: {np.mean(np.square(e)) / y.var():.2f}")
    print(f"{q} (n={len(y)}, theta {t.min():.1f}-{t.max():.1f}): LOO R  " + " | ".join(out))
