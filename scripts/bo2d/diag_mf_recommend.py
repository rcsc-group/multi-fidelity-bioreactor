"""Does mfbml (recursive KRR-LR-GPR, L7 -> L8 -> L9) find the L9 optimum with fewer L9 runs than an L9-only GP?
(2026-11-10, after the no-go of "search cheap, confirm at the raw cheap optimum": diag_level_optima.py.)

Frozen real truth (synth_bench mode real). Per replicate and per n9:
  - n9 L9 sites drawn uniformly in the domain (same sites for both arms; noisy evaluations at the noise floor);
  - chain also gets L7 + L8 at the 36 real L8 design sites plus at the n9 L9 sites;
  - each arm recommends argmax of the model mean of the objective on the 101x101 grid, subject to
    mean_t + Z_REC * sd_t <= ln cap (P(feasible) >= 0.8, as bo.recommend);
  - score = L9 regret of the recommendation (truth), and L9 feasibility (truth).
Cost: l9only n9 * 371; chain n9 * 376 + (36 + n9) * (5 + 26) minus the L9-site double count (L7/L8 at L9 sites are
  inside the 376 already) -> chain = n9 * 376 + 36 * 31.
No BO loop: this isolates the surrogate (mfbml) from the picker.
Usage: uv run python scripts/bo2d/diag_mf_recommend.py OBJ REP [cap]  -> prints CSV rows"""
import os
import sys
from pathlib import Path
import numpy as np

OBJ, REP = sys.argv[1], int(sys.argv[2])
CAP = float(sys.argv[3]) if len(sys.argv) > 3 else 0.014
os.environ["BO_OBJ"] = OBJ
sys.path.insert(0, str(Path(__file__).parent))
import synth_bench as S  # noqa: E402
from scripts.bo2d import bo  # noqa: E402

S.setup("real")
C = S.C
G = bo.grid(); lc = np.log(CAP)
k9t, t9t = S.F[9][0](G), S.F[9][1](G)
opt = k9t[t9t <= lc].max()
real = C.load_db("kLa_1T_10")
s8 = real[real.level == 8][["rpm", "theta"]].drop_duplicates().to_numpy(float)[:36]
lo, hi = C.DOMAIN[:, 0], C.DOMAIN[:, 1]


def rec(fk, ft):
    (mk, _) = fk(G); (mt, _) = ft(G)
    muk = mk[9][0]; mut, Vt = mt[9][0], mt[9][1]
    ok = mut + bo.Z_REC * np.sqrt(Vt) <= lc
    if not ok.any():
        return None
    return int(np.argmax(np.where(ok, muk, -np.inf)))


print("obj,rep,n9,arm,core_h,regret,feasible,rpm,theta", flush=True)
for n9 in (4, 6, 8, 12, 16, 24):
    rng = np.random.default_rng(1000 * REP + n9)
    X9 = lo + (hi - lo) * rng.random((n9, 2))
    K, T = {}, {}
    S.add(K, T, 9, X9, rng)
    S.add(K, T, 7, np.vstack([s8, X9]), rng); S.add(K, T, 8, np.vstack([s8, X9]), rng)
    for arm, cost in (("l9only", n9 * 371.0), ("chain", n9 * 376.0 + 36 * 31.0)):
        fk, ft = S.models(arm, K, T)
        i = rec(fk, ft)
        if i is None:
            print(f"{OBJ},{REP},{n9},{arm},{cost:.0f},nan,False,nan,nan", flush=True); continue
        reg = 1 - np.exp(k9t[i] - opt)
        print(f"{OBJ},{REP},{n9},{arm},{cost:.0f},{reg:.4f},{bool(t9t[i] <= lc)},{G[i,0]:.3f},{G[i,1]:.3f}", flush=True)
