"""Is the frozen BO truth fair to the chain? (2026-11-10, user: "mfbml was great in 1D, trash in 2D - smells".)
Hypothesis: the truth = independent per-level GPs on sparse data, so between data sites truth L8 and truth L9 are less
correlated than the real flow; the chain is judged against an artefact.
Falsified if: corr(truth L8, truth L9) on the grid ~= corr(observed L8, observed L9) at paired real sites.
Also prints the truth vs observation residual at the L9 sites (how well the truth honours the data).
Usage: uv run python scripts/bo2d/diag_truth_corr.py"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import synth_bench as S  # noqa: E402
from scripts.bo2d import bo  # noqa: E402
from scripts.bo2d import loop as LP  # noqa: E402

G = bo.grid()
for obj, q in (("kla", "kLa_1T_10"), ("mix", "dtmix_0.25")):
    os.environ["BO_OBJ"] = obj; S.setup("real")
    for name, qq, l8, l9 in ((q, q, S.F[8][0], S.F[9][0]), ("tau_95", "tau_95", S.F[8][1], S.F[9][1])):
        a, b = l8(G), l9(G)
        d = pd.read_csv(Path(__file__).parent / "state" / f"truth_snapshot_{qq}.csv")
        d = d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]
        d = LP.as_objective(d, qq)
        d["ly"] = np.log(d.y.astype(float)); d["key"] = d.rpm.round(2).astype(str) + "_" + d.theta.round(2).astype(str)
        p = d[d.level == 8].groupby("key").ly.mean().to_frame("y8").join(d[d.level == 9].groupby("key").ly.mean().to_frame("y9"), how="inner")
        d9 = d[d.level == 9]
        res9 = d9.ly.to_numpy() - l9(d9[["rpm", "theta"]].to_numpy(float))
        print(f"{obj} {name}: grid corr(L8,L9) truth = {np.corrcoef(a, b)[0,1]:.2f} | paired real sites n={len(p)} "
              f"corr obs = {np.corrcoef(p.y8, p.y9)[0,1] if len(p) > 2 else float('nan'):.2f} | "
              f"truth L9 residual at L9 data: rms {np.sqrt(np.mean(res9**2)):.3f} ln (sd of data {d9.ly.std():.3f})", flush=True)
