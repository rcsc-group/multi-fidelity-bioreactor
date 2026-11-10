"""Ceiling test (2026-11-10): how much compute can ANY better tau model save?

Arm "oracle" = the l9only arm of synth_bench.py (same truth snapshot, same acquisition, same costs, same seeds), except
that the tau_95 model is the TRUTH t9 itself (zero variance, zero noise). Its core-h to <= 5 % regret is the best that a
perfect constraint model can give the L9-only search. If oracle / l9only is near 1, multi-fidelity on tau cannot reach
~10x by itself and the method must also save on the objective search.
  uv run python scripts/bo2d/synth_ceiling.py kla|mix SEEDS (e.g. 21,22,23)
"""
import os
import sys
from pathlib import Path

import numpy as np

os.environ["BO_OBJ"] = sys.argv[1]
os.environ.setdefault("BO_CAPN", "40")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import bo  # noqa: E402
from scripts.bo2d import synth_bench as S  # noqa: E402

S.setup("real")
S.COST["oracle"] = S.COST["l9only"]
_models = S.models


def models(arm, K, T):
    if arm != "oracle":
        return _models(arm, K, T)
    mk = bo.fit_single(K[9], S.QK)
    z = lambda X: np.zeros(len(np.atleast_2d(X)))  # noqa: E731
    return (lambda X: bo.single_moments(mk, X)), (lambda X: ({9: (S.t9(np.atleast_2d(X)), z(X) + 1e-12, z(X))}, {9: z(X) + 1}))


S.models = models
for seed in (int(s) for s in sys.argv[2].split(",")):
    S.run_arm("oracle", 0.014, seed)
