"""Screen-then-L9 (user design, 2026-11-10 22:30: "screen with L8 very finely, then a few L9 jobs to map out the
landscape"): frozen-truth 2D BO, arm "screen".

Pre-registered BEFORE any run (2026-11-10 ~22:35):
- Screen, paid up front: L7 + L8 (L8 warm from L7, 26 core-h per site) at N maximin-LHS sites (seed 20261110),
  N from env SCREEN_N (runs at N = 40 and N = 120), plus L7 + L8 at the 4 seed L9 sites.
- Then the same loop as synth_taumf / l9only: aEI, L9 queries only (376 core-h, adds L7 + L8 + L9 at the site),
  BUT both models (objective AND tau) are the L7 -> L8 -> L9 chain (bo.fit_chain3, BO_THETA=1).
- Same truth snapshots, seeds 21-23, cap 0.014, BO_CAPN 120, 5 % regret target, recommend rule.
- Score: core-h to first feasible <= 5 % regret vs frozen l9only (kla 20776 / 10388 / 14840; mix 8904 / 13356 / 7420).
  PASS if screen < l9only in >= 2 of 3 seeds per objective. Saving = geometric mean of l9only / screen.
  uv run python scripts/bo2d/synth_screen.py kla|mix SEEDS      (env SCREEN_N)
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

os.environ["BO_OBJ"] = sys.argv[1]
os.environ.setdefault("BO_CAPN", "120")
os.environ.setdefault("BO_THETA", "1")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import bo  # noqa: E402
from scripts.bo2d import synth_bench as S  # noqa: E402

N, DOE_SEED = int(os.environ.get("SCREEN_N", "40")), 20261110
COST = {9: 376.0}
TAG = f"screen{N}"


def lhs(n):
    rng = np.random.default_rng(DOE_SEED)
    lo, hi = S.C.DOMAIN[:, 0], S.C.DOMAIN[:, 1]
    best, bd = None, -1.0
    for _ in range(500):
        u = (np.argsort(rng.random((2, n)), axis=1).T + rng.random((n, 2))) / n
        d = np.min(np.linalg.norm(u[:, None] - u[None], axis=2) + np.eye(n) * 9)
        if d > bd:
            best, bd = u, d
    return lo + best * (hi - lo)


def models(K, T):
    mk = bo.fit_chain3(K[7], K[8], K[9], S.QK)
    mt = bo.fit_chain3(T[7], T[8], T[9], "tau_95")
    return (lambda X: bo.chain_moments(mk, X)), (lambda X: bo.chain_moments(mt, X))


def run(cap, seed, max_rounds=200):
    rng = np.random.default_rng(seed)
    K, T = {}, {}
    seed9 = np.array(S.SEED9)
    S.add(K, T, 9, seed9, rng)
    X = np.vstack([lhs(N), seed9])
    S.add(K, T, 7, X, rng); S.add(K, T, 8, X, rng)
    spent = 26.0 * len(X)
    print(f"{TAG} seed {seed} screen core-h {spent:.0f}", flush=True)
    G = bo.grid(); lncap = np.log(cap)
    true_opt = S.k9(G)[S.t9(G) <= lncap].max()
    for rd in range(max_rounds):
        if spent >= S.CAP_CH:
            break
        Mk, Mt = models(K, T)
        Xe = pd.concat([K[l][["rpm", "theta"]] for l in K]).drop_duplicates().to_numpy(float)
        mk_e, _ = Mk(Xe); mt_e, _ = Mt(Xe)
        best = bo.effective_best(mk_e, mt_e, lncap)
        if best is None:
            best = float(np.log(K[9].y).min())
        mk, ck = Mk(G); mt, _ = Mt(G)
        picks = bo.pick_budget(bo.aei(mk, ck, mt, best, lncap, COST), G, S.BUDGET, COST)
        for x, _l in picks:
            for lev in (7, 8, 9):
                S.add(K, T, lev, x[None], rng)
            spent += COST[9]
        Mk, Mt = models(K, T)
        xb, _ = bo.recommend_m(Mk, Mt, K[9], lncap)
        if xb is None:
            print(f"{TAG} seed {seed} round {rd} core-h {spent:6.0f} no feasible recommendation", flush=True); continue
        feas = bool(S.t9(xb[None])[0] <= lncap); regret = 1 - np.exp(S.k9(xb[None])[0] - true_opt)
        print(f"{TAG} seed {seed} round {rd} core-h {spent:6.0f} regret {100 * regret:5.1f}% feasible {feas} "
              f"at {xb.round(2)} picks L9 {len(picks)}", flush=True)
        if feas and regret <= S.TARGET:
            print(f"HIT {TAG} seed {seed} core-h {spent:.0f}", flush=True); return spent
    print(f"MISS {TAG} seed {seed} core-h {spent:.0f}", flush=True); return None


if __name__ == "__main__":
    S.setup("real")
    for s in (int(a) for a in sys.argv[2].split(",")):
        run(0.014, s)
