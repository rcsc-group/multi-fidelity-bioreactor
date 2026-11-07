"""BO driver pieces (scripts/bo2d/bo.py), synthetic data only."""
import numpy as np
import pandas as pd
from scripts.bo2d import bo


def test_cei_zero_where_infeasible_and_positive_where_feasible():
    mk = np.array([1.0, 1.0]); sk = np.array([0.1, 0.1])
    mt = np.array([-5.0, 5.0]); st = np.array([0.01, 0.01])
    a = bo.cei(mk, sk, mt, st, best=0.5, lncap=0.0)
    assert a[0] > 0.4 and a[1] < 1e-12


def test_pick_batch_respects_exclusion_disc():
    G = bo.grid(); s = -np.linalg.norm(bo.scaled(G) - 0.5, axis=1) + 10
    B = bo.pick_batch(s, G, 3)
    d = np.linalg.norm(bo.scaled(B)[:, None] - bo.scaled(B)[None], axis=2)[np.triu_indices(3, 1)]
    assert len(B) == 3 and d.min() >= bo.R_EXCL


def test_next_level_cascade():
    x = np.array([[25.0, 4.0]])
    row = pd.DataFrame(dict(rpm=[25.0], theta=[4.0]))
    assert bo.next_level(x, {}) == 7
    assert bo.next_level(x, {7: row}) == 8
    assert bo.next_level(x, {7: row, 8: row}) == 9
    assert bo.next_level(x, {7: row, 8: row, 9: row}) is None


def test_best_feasible_uses_same_run_tau():
    dk = pd.DataFrame(dict(run=["a", "b"], y=[30.0, 20.0]))
    dt = pd.DataFrame(dict(run=["a", "b"], y=[0.2, 0.05]))
    v, r = bo.best_feasible(dk, dt, cap=0.1)
    assert r == "b" and abs(v - np.log(20.0)) < 1e-12
