"""Fair model ablation (2026-11-10). Replaces the scoring protocol of mf_ablation.py, which drew sites at random from a
pool where 10 of 30 L9 sites sit on the theta = 7 1D line (zero area in the domain, and L8 is densest there). The user
called out this theta = 7 draw bias three times.

PRE-REGISTERED (written before any result):
- L9 pool: frozen snapshot plus every finished live BO L9 run (bo1_l9_, bo3_l9_, bo4_l9_; N_WARM = 30), deduplicated
  by run. DROPPED: the 8 c2d_l9 warm starts (released after 12 settle cycles; the L9 wave needs about 27).
- Test sites: OFF-LINE L9 sites only (|theta - 7| > 0.05). The theta = 7 L9 sites never score.
- Training draw: n9 sites, of which ON7 (0 or 1, a factor) come from the theta = 7 line; the rest are off-line.
- L8 rows: a factor, L8LINE = keep (all L8 rows) or drop (L8 rows on theta = 7 removed except at a training L9 site).
- L7 rows: all of them (incl. the dense screen scr_l7_ runs that have finished).
- Variants: B1 (L9 only), U2-lin (paper), B2-lin, B2-quad, B3-lin, B3-quad (same code as mf_ablation.py).
- Score: RMSE (ln) over the held-out off-line sites; 20 draws per n9 in {4, 8}; also 90 % coverage.
Usage: uv run python scripts/bo2d/mf_ablation_fair.py QOI N9 ON7 L8LINE  -> CSV rows"""
import sys
import time
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts")); sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))
warnings.filterwarnings("ignore")
Q, N9, ON7, L8LINE = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]

from scripts.bo2d import mf_ablation as A  # noqa: E402  (reuse its model builders verbatim)
from scripts.bo2d import bo  # noqa: E402
from scripts.bo2d import loop as LP  # noqa: E402

ON = lambda s: np.isclose(s.theta.astype(float), 7.0, atol=0.05)


def load(q):
    d = A.load(q)
    live = LP.rows(q, ["bo1_l9_", "bo3_l9_", "bo4_l9_", "scr_l7_"], {7, 9})
    live = live[np.isfinite(live.y.astype(float)) & (live.y.astype(float) > 0)]
    if len(live):
        live = LP.as_objective(live, q); live["ly"] = np.log(live.y.astype(float))
        live["key"] = live.rpm.round(2).astype(str) + "_" + live.theta.round(2).astype(str)
        d = pd.concat([d, live], ignore_index=True).drop_duplicates("run")
    return d[~d.run.str.startswith("c2d_l9")]


d = load(Q); floor = A.C.noise_floor(Q); xy = A.xy
d9 = d[d.level == 9]; off = d9[~ON(d9)].key.unique(); on = d9[ON(d9)].key.unique()
print(f"# {Q}: L7 {int((d.level == 7).sum())}, L8 {int((d.level == 8).sum())} "
      f"({int(ON(d[d.level == 8]).sum())} on theta=7), L9 off-line sites {len(off)}, on-line sites {len(on)}", flush=True)
X7, y7 = xy(d[d.level == 7]); lf7 = A.krr(X7, y7)
print("qoi,n9,on7,l8line,rep,variant,rmse,cover90,seconds,ntest", flush=True)
CACHE = {}
for rep in range(20):
    rng = np.random.default_rng(1000 * rep + 10 * N9 + ON7)
    tr = list(rng.choice(off, N9 - ON7, replace=False)) + list(rng.choice(on, ON7, replace=False))
    Atr = d9[d9.key.isin(tr)]; B = d9[~d9.key.isin(tr) & ~ON(d9)]
    d8 = d[d.level == 8]
    if L8LINE == "drop":
        d8 = d8[~ON(d8) | d8.key.isin(tr)]
    k8 = tuple(d8.index)                      # L8 fits depend only on which L8 rows are kept: fit once per set
    if k8 not in CACHE:
        X8, y8 = xy(d8); CACHE[k8] = (X8, y8, A.krr(X8, y8), A.fb(X8, y8, lf7, False, floor))
    X8, y8, lf8, fb8 = CACHE[k8]
    Xa, ya = xy(Atr); Xb, yb = xy(B)
    V = {
        "B1": lambda: bo.gauss(A.fb(Xa, ya, None, False, floor), Xb),
        "U2-lin": lambda: A.up_pred(A.upstream(Xa, ya, X8, y8, "linear", lf=A._Fixed(lf8.predict)), Xb),
        "B2-lin": lambda: bo.gauss(A.fb(Xa, ya, lf8, False, floor), Xb),
        "B2-quad": lambda: bo.gauss(A.fb(Xa, ya, lf8, True, floor), Xb),
        "B3-lin": lambda: bo.gauss(A.fb(Xa, ya, bo._Mean(fb8), False, floor), Xb),
        "B3-quad": lambda: bo.gauss(A.fb(Xa, ya, bo._Mean(fb8), True, floor), Xb),
    }
    for name, f in V.items():
        t0 = time.time()
        try:
            mu, sd = f()
            rmse = float(np.sqrt(np.mean((mu - yb) ** 2))); cov = float(np.mean(np.abs(mu - yb) <= A.Z90 * sd))
        except Exception as e:  # report, never hide
            rmse, cov = np.nan, np.nan; print(f"# {name} rep {rep}: {type(e).__name__}: {e}", flush=True)
        print(f"{Q},{N9},{ON7},{L8LINE},{rep},{name},{rmse:.4f},{cov:.3f},{time.time() - t0:.1f},{len(yb)}", flush=True)
