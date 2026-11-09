"""Offline check of the BO driver: augmented EI (Huang et al. 2006) on synthetic multi-fidelity functions (2026-11-07).

Truth (ln units) shaped like our data: kLa9 rises with rpm and theta plus a resonance bump; tau9 rises with both.
L8 = L9 + a theta-dependent offset (L8 over-predicts kLa at low theta, as the data show); L7 = L8 + offset + wiggle.
Noise sd 0.08 per evaluation. Both arms: same 4 L9 seed sites (SEED9), the same core-h budget per round
(4 L9-equivalents), the same aEI code. l9only: levels {9}, cost 97 (L8 spin-up 20 + L9 warm 77).
chain: levels {7, 8, 9}; costs with prerequisites: L7 5, L8 22 (L7 + L8 warm), L9 99 (L7 + L8 + L9);
the prerequisite runs' data are added. Chain seed also has L7 and L8 at the 36 real L8 sites.
Criterion (stated before running): chain reaches <= 5 % regret with fewer core-h than l9only in >= 2 of 3 seeds.
Stop rule (2026-11-07 fix): each arm runs until it first reaches <= 5 % regret or spends CAP_CH core-h,
not a fixed number of rounds (the 10-round rule gave the chain ~800 core-h vs 3880 for l9only).
Mode "held": a second truth (different bump, slopes and L8/L7 discrepancy) and new seeds 11-13,
written before the dev run with the cost cap was seen; it is the test of the locked method.
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import bo  # noqa: E402
C = bo.C
SD = 0.08
QK = "kLa_1T_10"
SD_K = SD_T = SD   # mode 'real' sets the measured floors per QoI
SEED9 = [(16.35, 2.33), (19.14, 2.89), (23.14, 2.13), (36.35, 5.24)]
COST = {"l9only": {9: 97.0}, "chain": {7: 5.0, 8: 22.0, 9: 99.0}}
STEP = {7: 5.0, 8: 17.0, 9: 77.0}
BUDGET = 4 * 97.0


def k9(X):
    r, t = X[:, 0], X[:, 1]
    return np.log(15) + 0.035 * (r - 15) + 0.06 * (t - 2) + 0.9 * np.exp(-((r - 37.5) / 2.5) ** 2) * (t / 7)
def t9(X):
    r, t = X[:, 0], X[:, 1]
    return np.log(0.004) + 0.12 * (r - 15) + 0.15 * (t - 2)
def k8(X): return k9(X) + 0.3 + 0.9 * (7 - X[:, 1]) / 5
def t8(X): return t9(X) - 0.2
def k7(X): return k8(X) + 0.2 + 0.1 * np.sin(X[:, 0] / 3)
def t7(X): return t8(X) - 0.2
F = {7: (k7, t7), 8: (k8, t8), 9: (k9, t9)}
CAP_CH = 3880.0
METHOD = os.environ.get("BO_METHOD", "aei")   # "mf": MF-GP-UCB level rule for the chain (bo.mf_pick)
bo.CHAIN_THETA = os.environ.get("BO_THETA", "0") == "1"   # basis [1, f_lower, theta] for L8/L9 layers
TARGET = 0.05


def _real_truth(q, lev):
    """Truth from the REAL data (2026-11-09): FullBayes GP posterior mean (ln) fitted to every run of level lev for
    QoI q, independent per level (the chain's MF structure is NOT imposed on the truth). Mixing as the rate 1/t."""
    from scripts.bo2d import loop as LP
    d = pd.concat([C.load_db(q), LP.rows(q, ["bo1_", "bo4_", "bo_l7_"], {7, 8, 9})]).drop_duplicates("run")
    d = d[(d.level == lev) & np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]
    d = LP.as_objective(d, q)
    X, y = d[["rpm", "theta"]].to_numpy(float), np.log(d.y.to_numpy(float))
    m = bo.FullBayesKRRLRGP(C.DOMAIN, basis="ordinary", grid=C.hf_grid(C.noise_floor(q), float(np.std(y)))).train(X, y)
    print(f"truth {q} L{lev}: n={len(y)} sd={np.std(y):.3f}", flush=True)
    return lambda Z: m.predict_mean(np.atleast_2d(Z))


def setup(mode):
    """mode 'held': unseen truth for the test of the locked method. mode 'real': truth from the real data,
    objective from env BO_OBJ (kla | mix), measured costs."""
    global k9, t9, k8, t8, k7, t7, F, SD_K, SD_T, COST, BUDGET, CAP_CH, QK
    if mode == "real":
        qk = QK = {"kla": "kLa_1T_10", "mix": "dtmix_0.25"}[os.environ.get("BO_OBJ", "kla")]
        k9, k8, k7 = (_real_truth(qk, l) for l in (9, 8, 7))
        t9, t8, t7 = (_real_truth("tau_95", l) for l in (9, 8, 7))
        F = {7: (k7, t7), 8: (k8, t8), 9: (k9, t9)}
        SD_K, SD_T = C.noise_floor(qk), C.noise_floor("tau_95")
        COST = {"l9only": {9: 371.0}, "chain": {7: 5.0, 8: 26.0, 9: 376.0}}   # measured (sacct, 2026-11-08)
        BUDGET = 4 * 371.0; CAP_CH = 40 * 371.0
        return (21, 22, 23)
    if mode != "held":
        return (1, 2, 3)
    def k9(X):
        r, t = X[:, 0], X[:, 1]
        return np.log(12) + 0.028 * (r - 15) + 0.09 * (t - 2) + 0.6 * np.exp(-((r - 29) / 3.5) ** 2 - ((t - 4) / 2) ** 2)
    def t9(X):
        r, t = X[:, 0], X[:, 1]
        return np.log(0.005) + 0.10 * (r - 15) + 0.22 * (t - 2)
    def k8(X): return k9(X) + 0.5 + 0.6 * (7 - X[:, 1]) / 5 - 0.01 * (X[:, 0] - 25)
    def t8(X): return t9(X) - 0.1
    def k7(X): return k8(X) + 0.4 + 0.15 * np.cos(X[:, 1])
    def t7(X): return t8(X) - 0.15
    F = {7: (k7, t7), 8: (k8, t8), 9: (k9, t9)}
    return (11, 12, 13)


def evaluate(level, X, rng):
    fk, ft = F[level]; n = len(X)
    runs = [f"r{level}_{rng.integers(1e12)}" for _ in range(n)]
    dk = pd.DataFrame(dict(level=level, rpm=X[:, 0], theta=X[:, 1], run=runs, y=np.exp(fk(X) + SD_K * rng.standard_normal(n))))
    dt = dk.copy(); dt["y"] = np.exp(ft(X) + SD_T * rng.standard_normal(n))
    return dk, dt


def add(K, T, lev, X, rng):
    dk, dt = evaluate(lev, X, rng)
    K[lev] = pd.concat([K.get(lev), dk], ignore_index=True); T[lev] = pd.concat([T.get(lev), dt], ignore_index=True)


def models(arm, K, T):
    if arm == "l9only":
        mk, mt = bo.fit_single(K[9], QK), bo.fit_single(T[9], "tau_95")
        return (lambda X: bo.single_moments(mk, X)), (lambda X: bo.single_moments(mt, X))
    mk = bo.fit_chain3(K[7], K[8], K[9], QK); mt = bo.fit_chain3(T[7], T[8], T[9], "tau_95")
    return (lambda X: bo.chain_moments(mk, X)), (lambda X: bo.chain_moments(mt, X))


def run_arm(arm, cap, seed, max_rounds=80):
    rng = np.random.default_rng(seed)
    spent0 = [0.0]
    real = C.load_db("kLa_1T_10"); s8 = real[real.level == 8][["rpm", "theta"]].drop_duplicates().to_numpy(float)[:36]
    K, T = {}, {}
    add(K, T, 9, np.array(SEED9), rng)
    if arm == "chain":
        if os.environ.get("BO_FAIRSEED", "0") == "1":   # 2026-11-08 fair seed: = real chain arm (loop.arm_data)
            from scripts.bo2d import loop as LP
            s7 = LP.arm_data("chain", "tau_95")[7][["rpm", "theta"]].drop_duplicates().to_numpy(float)
            add(K, T, 7, s7, rng); add(K, T, 8, np.array(SEED9), rng)
            spent0[0] = 5.0 * (len(s7) - len(SEED9)) + 26.0 * len(SEED9)   # reused L7 + L7/L8 at the seed sites
        else:
            add(K, T, 7, s8, rng); add(K, T, 8, s8, rng)
    G = bo.grid(); lncap = np.log(cap); spent = spent0[0]
    true_opt = k9(G)[t9(G) <= lncap].max()
    for rd in range(max_rounds):
        if spent >= CAP_CH:
            break
        Mk, Mt = models(arm, K, T)
        Xe = pd.concat([K[l][["rpm", "theta"]] for l in K]).drop_duplicates().to_numpy(float)
        mk_e, _ = Mk(Xe); mt_e, _ = Mt(Xe)
        best = bo.effective_best(mk_e, mt_e, lncap)
        if best is None:
            best = float(np.log(K[9].y).min())
        mk, ck = Mk(G); mt, _ = Mt(G)
        if arm == "chain" and METHOD == "mf":
            picks = bo.mf_pick(mk, mt, best, lncap, G)
        else:
            picks = bo.pick_budget(bo.aei(mk, ck, mt, best, lncap, COST[arm]), G, BUDGET, COST[arm])
        nl = {7: 0, 8: 0, 9: 0}
        for x, l in picks:
            for lev in ([l] if arm == "l9only" else range(7, l + 1)):
                add(K, T, lev, x[None], rng)
            spent += COST[arm][l]; nl[l] += 1
        Mk, Mt = models(arm, K, T)
        xb, _ = bo.recommend_m(Mk, Mt, K[9], lncap)
        if xb is None:
            print(f"{arm:6s} seed {seed} round {rd} core-h {spent:6.0f} no feasible recommendation", flush=True); continue
        feas = bool(t9(xb[None])[0] <= lncap); regret = 1 - np.exp(k9(xb[None])[0] - true_opt)
        print(f"{arm:6s} seed {seed} round {rd} core-h {spent:6.0f} regret {100 * regret:5.1f}% feasible {feas} "
              f"at {xb.round(2)} picks L7/L8/L9 {nl[7]}/{nl[8]}/{nl[9]}", flush=True)
        if feas and regret <= TARGET:
            print(f"HIT {arm} seed {seed} core-h {spent:.0f}", flush=True); return spent
    print(f"MISS {arm} seed {seed} core-h {spent:.0f}", flush=True); return None


if __name__ == "__main__":
    cap = float(sys.argv[1]) if len(sys.argv) > 1 else 0.03
    mode = sys.argv[2] if len(sys.argv) > 2 else "dev"
    seeds = setup(mode)
    if len(sys.argv) > 3:
        seeds = tuple(int(a) for a in sys.argv[3].split(","))
    arms = sys.argv[4].split(",") if len(sys.argv) > 4 else ["chain", "l9only"]
    G = bo.grid()
    print("true constrained optimum:", G[np.argmax(np.where(t9(G) <= np.log(cap), k9(G), -np.inf))], flush=True)
    win = 0
    for seed in seeds:
        h = {arm: run_arm(arm, cap, seed) for arm in arms}
        if len(arms) < 2:
            continue
        ok = h["chain"] is not None and (h["l9only"] is None or h["chain"] < h["l9only"])
        win += ok
        print(f"SEED {seed} chain {h['chain']} l9only {h['l9only']} chain_wins {ok}", flush=True)
    print(f"VERDICT mode {mode}: chain wins {win}/3 -> {'PASS' if win >= 2 else 'FAIL'}", flush=True)
