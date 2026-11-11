"""Multi-fidelity on the CONSTRAINT only (user go 2026-11-10 18:55): frozen-truth 2D BO, arm "taumf".

Pre-registered BEFORE any run (2026-11-10 ~19:00):
- Same truth snapshots, acquisition (aEI, levels {9}), recommend rule, seeds 21-23, cap 0.014, BO_CAPN 120, 5 % regret
  target, as the l9only arm of synth_bench.py ('real' mode). Objective model: L9-only GP (= l9only).
- Constraint (tau_95) model: L7 -> L8 -> L9 chain (bo.fit_chain3, BO_THETA=1), i.e. the fair ablation's winner.
- Cheap design, paid up front (charged in core-h): L7 at 30 maximin-LHS sites (seed 20261110) + L8 at 10 of them
  (greedy maximin subset) = the user's 30 / 10 split; plus L7 + L8 at the 4 seed L9 sites. Each L9 query costs 376
  (L7 + L8 + L9 at the site, measured) and adds all three levels' data.
- Score per objective: core-h to the first feasible <= 5 % regret recommendation, vs the frozen l9only arm
  (kla 20776 / 10388 / 14840; mix 8904 / 13356 / 7420) and the perfect-tau ceiling (kla 5936 / 4452 / 4452;
  mix 1484 / 1484 / 4452). PASS if taumf < l9only in >= 2 of 3 seeds. Saving = l9only / taumf (geometric mean).
  uv run python scripts/bo2d/synth_taumf.py kla|mix SEEDS
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

N7, N8, DOE_SEED = 30, 10, 20261110
COST = {9: 376.0}


def doe():
    rng = np.random.default_rng(DOE_SEED)
    lo, hi = S.C.DOMAIN[:, 0], S.C.DOMAIN[:, 1]
    best, bd = None, -1.0
    for _ in range(2000):   # maximin LHS
        u = (np.argsort(rng.random((2, N7)), axis=1).T + rng.random((N7, 2))) / N7
        d = np.min(np.linalg.norm(u[:, None] - u[None], axis=2) + np.eye(N7) * 9)
        if d > bd:
            best, bd = u, d
    sel = [0]
    while len(sel) < N8:   # greedy maximin subset for L8
        dist = np.min(np.linalg.norm(best[:, None] - best[sel][None], axis=2), axis=1)
        dist[sel] = -1
        sel.append(int(np.argmax(dist)))
    X7 = lo + best * (hi - lo)
    return X7, X7[sel]


def models(K, T):
    mk = bo.fit_single(K[9], S.QK)
    mt = bo.fit_chain3(T[7], T[8], T[9], "tau_95")
    return (lambda X: bo.single_moments(mk, X)), (lambda X: bo.chain_moments(mt, X))


def run(cap, seed, max_rounds=200):
    rng = np.random.default_rng(seed)
    K, T = {}, {}
    S.add(K, T, 9, np.array(S.SEED9), rng)
    X7, X8 = doe()
    seed9 = np.array(S.SEED9)
    S.add(K, T, 7, np.vstack([X7, seed9]), rng)
    S.add(K, T, 8, np.vstack([X8, seed9]), rng)
    spent = 5.0 * N7 + 21.0 * N8 + 26.0 * len(seed9)   # L7 5; L8 warm step 21 (26 incl. its L7)
    print(f"taumf seed {seed} cheap design core-h {spent:.0f}", flush=True)
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
            print(f"taumf seed {seed} round {rd} core-h {spent:6.0f} no feasible recommendation", flush=True); continue
        feas = bool(S.t9(xb[None])[0] <= lncap); regret = 1 - np.exp(S.k9(xb[None])[0] - true_opt)
        print(f"taumf seed {seed} round {rd} core-h {spent:6.0f} regret {100 * regret:5.1f}% feasible {feas} "
              f"at {xb.round(2)} picks L9 {len(picks)}", flush=True)
        if feas and regret <= S.TARGET:
            print(f"HIT taumf seed {seed} core-h {spent:.0f}", flush=True); return spent
    print(f"MISS taumf seed {seed} core-h {spent:.0f}", flush=True); return None


if __name__ == "__main__":
    S.setup("real")
    for s in (int(a) for a in sys.argv[2].split(",")):
        run(0.014, s)
