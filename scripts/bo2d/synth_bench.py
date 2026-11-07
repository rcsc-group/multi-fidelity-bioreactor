"""Offline check of the BO driver on synthetic multi-fidelity functions shaped like our data (2026-11-07).

Truth (ln units): kLa9 rises with rpm and theta plus a resonance bump near 37 rpm; tau9 rises with rpm and theta.
L8 = L9 + a theta-dependent offset (L8 over-predicts kLa at low theta, as the data show, diary 2026-11-06);
L7 = L8 + a further offset and a wiggle. Noise sd 0.08 on every evaluation (~8 %, the warm/cold L9 kLa spread).
Seeds: the 4 gated campaign L9 sites (SEED9) for both arms; the chain also gets L8 and L7 at the 36 real L8 sites.
Pass criteria (stated before running): (1) both arms end within 5 % in kLa of the true constrained optimum;
(2) report core-h of each arm at first reaching 5 % regret. Costs per evaluation: L7 3, L8 17, L9 100 core-h.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import bo  # noqa: E402
C = bo.C
COST = {7: 3.0, 8: 17.0, 9: 100.0}
SD = 0.08


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
# Common L9 seed of both arms: the 4 campaign L9 sites that pass the settled gate (2026-11-07). The theta-7 line
# rows are left out: they hold the constrained optimum, so with them both arms start solved.
SEED9 = [(16.35, 2.33), (19.14, 2.89), (23.14, 2.13), (36.35, 5.24)]


def evaluate(level, X, rng, tag):
    fk, ft = F[level]
    n = len(X)
    runs = [f"{tag}{level}_{rng.integers(1e9)}" for _ in range(n)]
    dk = pd.DataFrame(dict(level=level, rpm=X[:, 0], theta=X[:, 1], run=runs, y=np.exp(fk(X) + SD * rng.standard_normal(n))))
    dt = dk.copy(); dt["y"] = np.exp(ft(X) + SD * rng.standard_normal(n))
    return dk, dt


def run_arm(arm, cap, rounds, q, seed):
    rng = np.random.default_rng(seed)
    real_k = C.load_db("kLa_1T_10")
    s9 = np.array(SEED9)
    s8 = real_k[real_k.level == 8][["rpm", "theta"]].drop_duplicates().to_numpy(float)[:36]
    K, T = {}, {}
    K[9], T[9] = evaluate(9, s9, rng, "s")
    if arm == "chain":
        for lev in (7, 8):
            K[lev], T[lev] = evaluate(lev, s8, rng, "s")
    G = bo.grid(); lncap = np.log(cap); spent = 0.0; hist = []
    true_opt = k9(G)[t9(G) <= lncap].max()
    for rd in range(rounds):
        best, run = bo.best_feasible(K[9], T[9], cap)
        if arm == "l9only":
            mk, mt = bo.fit_single(K[9], "kLa_1T_10"), bo.fit_single(T[9], "tau_95")
        else:
            mk = bo.fit_chain(K[7], K[8], K[9], "kLa_1T_10")[0]; mt = bo.fit_chain(T[7], T[8], T[9], "tau_95")[0]
        a = bo.cei(*bo.gauss(mk, G), *bo.gauss(mt, G), best, lncap)
        B = bo.pick_batch(a, G, q)
        for x in B:
            lev = 9 if arm == "l9only" else bo.next_level(x[None], K)
            if lev is None:
                lev = 9
            dk, dt = evaluate(lev, x[None], rng, "b")
            K[lev] = pd.concat([K.get(lev), dk], ignore_index=True); T[lev] = pd.concat([T.get(lev), dt], ignore_index=True)
            spent += COST[lev]
        if arm == "l9only":
            mk, mt = bo.fit_single(K[9], "kLa_1T_10"), bo.fit_single(T[9], "tau_95")
        else:
            mk = bo.fit_chain(K[7], K[8], K[9], "kLa_1T_10")[0]; mt = bo.fit_chain(T[7], T[8], T[9], "tau_95")[0]
        xb, _ = bo.recommend(mk, mt, K[9], lncap)
        feas = bool(t9(xb[None])[0] <= lncap)
        regret = 1 - np.exp(k9(xb[None])[0] - true_opt)
        hist.append((rd, spent, regret, feas))
        print(f"{arm:6s} seed {seed} round {rd} core-h {spent:6.0f} regret {100 * regret:5.1f}% feasible {feas} at {xb.round(2)}", flush=True)
    return hist


if __name__ == "__main__":
    G = bo.grid()
    cap = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
    print("true constrained optimum:", G[np.argmax(np.where(t9(G) <= np.log(cap), k9(G), -np.inf))], flush=True)
    for arm in ("l9only", "chain"):
        run_arm(arm, cap, rounds=int(sys.argv[2]) if len(sys.argv) > 2 else 6, q=4, seed=1)
