"""Budget-split study (user's idea, 2026-11-10): for a fixed core-h budget B, which split between L7, L8 and L9 lets mfbml
(recursive KRR-LR-GPR, bo.fit_chain3 with the theta basis) recommend the best constrained optimum? Zero new simulations:
frozen real truth (synth_bench mode real), noisy evaluations at the measured noise floors.

PRE-REGISTERED (written before any result):
- Design: nested random sites, uniform in the domain. n9 L9 sites; L8 at the n9 sites + n8x extra sites; L7 at all L8
  sites + n7x extra sites. Costs incl. prerequisites: L9 site 376, L8-only site 26, L7-only site 5 (measured).
  Split = cost fractions (f7, f8, f9): n9 = floor(f9 B / 376), n8x = floor(f8 B / 26), n7x = floor(f7 B / 5);
  n7x, n8x capped at 300 (cost of the cap reported as spent). Need n9 >= 2.
- Baseline: l9only GP on floor(B / 371) L9 sites (the L9-only arm's cost per run).
- Recommendation (unchanged method, bo.Z_REC): argmax of model-mean objective on the 101 x 101 grid s.t.
  mean_t + Z_REC sd_t <= ln cap. cap = 0.014.
- Score: L9 truth regret of the recommendation; success = regret <= target AND truth ln tau_95 <= ln cap + 0.017 (the
  tau_95 noise floor; a boundary optimum cannot be judged finer than the measurement).
- Reps 0-15 choose the best split per budget (dev); reps 16-31 confirm it once (no re-tuning).
- Main objective: mix (t_mix,25 as the rate; truth L8-L9 corr 0.80 vs observed 0.83). kla secondary (truth unfair: -0.03
  vs 0.38).
Usage: uv run python scripts/bo2d/split_study.py OBJ BUDGET REP  -> CSV rows on stdout"""
import os
import sys
from pathlib import Path
import numpy as np

OBJ, B, REP = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
os.environ["BO_OBJ"] = OBJ
os.environ["BO_THETA"] = "1"
sys.path.insert(0, str(Path(__file__).parent))
import synth_bench as S  # noqa: E402
from scripts.bo2d import bo  # noqa: E402

bo.CHAIN_THETA = True
S.setup("real"); C = S.C
CAP = 0.014; lc = np.log(CAP); TOL = 0.017; NCAP = 300
G = bo.grid(); k9t, t9t = S.F[9][0](G), S.F[9][1](G); opt = k9t[t9t <= lc].max()
lo, hi = C.DOMAIN[:, 0], C.DOMAIN[:, 1]
SPLITS = [(0.0, 0.0, 1.0)] + [(round((1 - f9) * (1 - s8), 3), round((1 - f9) * s8, 3), f9)
                               for f9 in (0.4, 0.55, 0.7, 0.85) for s8 in (0.0, 0.5, 1.0)]


def rec(fk, ft):
    (mk, _) = fk(G); (mt, _) = ft(G)
    ok = mt[9][0] + bo.Z_REC * np.sqrt(mt[9][1]) <= lc
    return int(np.argmax(np.where(ok, mk[9][0], -np.inf))) if ok.any() else None


def row(arm, f, n, spent, i):
    if i is None:
        reg, exc, x = np.nan, np.nan, (np.nan, np.nan)
    else:
        reg, exc, x = 1 - np.exp(k9t[i] - opt), t9t[i] - lc, G[i]
    print(f"{OBJ},{B:.0f},{REP},{arm},{f[0]},{f[1]},{f[2]},{n[0]},{n[1]},{n[2]},{spent:.0f},{reg:.4f},{exc:.4f},"
          f"{bool(i is not None and exc <= TOL)},{x[0]:.3f},{x[1]:.3f}", flush=True)


print("obj,budget,rep,arm,f7,f8,f9,n7x,n8x,n9,spent,regret,excess,feasible_tol,rpm,theta", flush=True)
rng0 = np.random.default_rng(10_000 * REP + int(B))
pool = lo + (hi - lo) * rng0.random((2000, 2))          # one site pool per (budget, rep): splits are paired
n9b = int(B // 371)
rng = np.random.default_rng(rng0.integers(1 << 30)); K, T = {}, {}
S.add(K, T, 9, pool[:n9b], rng)
fk, ft = S.models("l9only", K, T)
row("l9only", (0, 0, 1), (0, 0, n9b), n9b * 371.0, rec(fk, ft))
for f in SPLITS:
    n9 = int(f[2] * B // 376); n8x = min(int(f[1] * B // 26), NCAP); n7x = min(int(f[0] * B // 5), NCAP)
    if n9 < 2:
        continue
    X9 = pool[:n9]; X8 = pool[:n9 + n8x]; X7 = pool[:n9 + n8x + n7x]
    rng = np.random.default_rng(rng0.integers(1 << 30)); K, T = {}, {}
    S.add(K, T, 9, X9, rng); S.add(K, T, 8, X8, rng); S.add(K, T, 7, X7, rng)
    fk, ft = S.models("chain", K, T)
    row("chain", f, (n7x, n8x, n9), n9 * 376.0 + n8x * 26.0 + n7x * 5.0, rec(fk, ft))
