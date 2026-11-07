"""Falsifiable check (2026-11-07): does the chain's L9 layer inherit L8's low-theta bias, and does basis [1, f8, theta] fix it?
Data as in synth_bench dev seed: L7+L8 at 36 real L8 sites, SEED9 + n9 extra L9 points (random). Report the
recommended point and its true regret for both bases over 5 draws. Prediction: [1,f8] -> theta low, larger regret."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import bo, synth_bench as S  # noqa: E402

cap = 0.03; lncap = np.log(cap); G = bo.grid()
true_opt = S.k9(G)[S.t9(G) <= lncap].max()
real = bo.C.load_db("kLa_1T_10"); s8 = real[real.level == 8][["rpm", "theta"]].drop_duplicates().to_numpy(float)[:36]
for n9 in (4, 12):
    for draw in range(5):
        rng = np.random.default_rng(100 + draw)
        K, T = {}, {}
        X9 = np.vstack([np.array(S.SEED9), bo.C.DOMAIN[:, 0] + rng.random((n9, 2)) * np.ptp(bo.C.DOMAIN, axis=1)])
        S.add(K, T, 9, X9, rng); S.add(K, T, 7, s8, rng); S.add(K, T, 8, s8, rng)
        out = []
        for theta in (False, True):
            bo.CHAIN_THETA = theta
            mk = bo.fit_chain3(K[7], K[8], K[9], "kLa_1T_10"); mt = bo.fit_chain3(T[7], T[8], T[9], "tau_95")
            # what the model believes over the whole grid (recommend_m only scans evaluated L9 points: insensitive)
            mu = bo.chain_moments(mk, G)[0][9][0]; tm = bo.chain_moments(mt, G)[0][9][0]
            ok = tm <= lncap
            xb = G[ok][np.argmax(mu[ok])]
            reg = 1 - np.exp(S.k9(xb[None])[0] - true_opt)
            rmse = np.sqrt(np.mean((mu - S.k9(G)) ** 2))
            out.append(f"{'theta' if theta else 'plain'} argmax={xb.round(2)} regret {100*reg:5.1f}% rmse {rmse:.3f}")
        print(f"n9={4+n9:2d} draw {draw}: " + " | ".join(out), flush=True)
