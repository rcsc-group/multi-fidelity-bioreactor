"""Why are 12-62 % of the chain's recommendations infeasible when truth corr(L8, L9) of tau_95 is 0.97? (2026-11-10)
Hypothesis: the chain's tau_95 predictive at L9 is biased or overconfident (bug or model defect), not a data limit.
Test: same data as diag_mf_recommend (L7/L8 at 36 real sites + L9 sites; n9 L9 sites). On the grid, compare the
L9 predictive of tau_95 with the truth: bias, RMSE, 90 % coverage; also near the cap only (|truth - ln cap| < 0.3).
Usage: uv run python scripts/bo2d/diag_chain_tau.py"""
import os
import sys
from pathlib import Path
import numpy as np

os.environ["BO_OBJ"] = "mix"
sys.path.insert(0, str(Path(__file__).parent))
import synth_bench as S  # noqa: E402
from scripts.bo2d import bo  # noqa: E402

S.setup("real"); C = S.C
G = bo.grid(); t9 = S.F[9][1](G); lc = np.log(0.014); near = np.abs(t9 - lc) < 0.3
real = C.load_db("kLa_1T_10"); s8 = real[real.level == 8][["rpm", "theta"]].drop_duplicates().to_numpy(float)[:36]
lo, hi = C.DOMAIN[:, 0], C.DOMAIN[:, 1]
print("n9 rep arm | bias rmse cover90 | near-cap: bias rmse cover90 | rho8 rho9", flush=True)
for n9 in (8, 16):
    for rep in range(3):
        rng = np.random.default_rng(1000 * rep + n9)
        X9 = lo + (hi - lo) * rng.random((n9, 2)); K, T = {}, {}
        S.add(K, T, 9, X9, rng); S.add(K, T, 7, np.vstack([s8, X9]), rng); S.add(K, T, 8, np.vstack([s8, X9]), rng)
        for arm in ("l9only", "chain"):
            if arm == "chain":
                ms = bo.fit_chain3(T[7], T[8], T[9], "tau_95"); mom, _ = bo.chain_moments(ms, G)
                rr = f"{bo.rho_of(ms[1]):.2f} {bo.rho_of(ms[2]):.2f}"
            else:
                m = bo.fit_single(T[9], "tau_95"); mom, _ = bo.single_moments(m, G); rr = ""
            mu, V = mom[9][0], mom[9][1]; e = mu - t9; z = np.abs(e) / np.sqrt(V)
            f = lambda s: f"{e[s].mean():+.3f} {np.sqrt((e[s]**2).mean()):.3f} {np.mean(z[s] < 1.645):.2f}"
            print(f"{n9} {rep} {arm} | {f(slice(None))} | {f(near)} | {rr}", flush=True)
