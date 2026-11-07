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


def test_aei_prefers_cheap_level_when_correlated_and_l9_when_not():
    one = np.ones(3)
    mom = {7: (one, one, 0.01 * one), 9: (one, one, 0.01 * one)}
    mt = {9: (-5 * one, 0.01 * one, 0.0 * one)}
    cost = {7: 5.0, 9: 100.0}
    hi = bo.aei(mom, {7: 0.95 * one, 9: one}, mt, 0.5, 0.0, cost)
    lo = bo.aei(mom, {7: 0.01 * one, 9: one}, mt, 0.5, 0.0, cost)
    assert (hi[7] > hi[9]).all() and (lo[7] < lo[9]).all()


def test_pick_budget_stays_within_budget():
    G = bo.grid(); s = {9: np.ones(len(G)), 7: 0.5 * np.ones(len(G))}
    picks = bo.pick_budget(s, G, budget=230.0, cost={9: 100.0, 7: 5.0})
    spent = sum({9: 100.0, 7: 5.0}[l] for _, l in picks)
    assert spent <= 230.0 and sum(l == 9 for _, l in picks) == 2


def test_mf_pick_runs_lowest_uncertain_level():
    """At the EI argmax, L7 is learned (sd < noise) and L8 is not: the rule must run L8."""
    import numpy as np
    from scripts.bo2d import bo
    G = bo.grid()
    n = len(G)
    mu9 = np.zeros(n); mu9[n // 2] = 2.0
    big, small, noise = np.full(n, 0.5), np.full(n, 1e-4), np.full(n, 0.01)
    mom_k = {7: (mu9, small, noise), 8: (mu9, big, noise), 9: (mu9, big, noise)}
    mom_t = {9: (np.full(n, -10.0), small, noise)}
    picks = bo.mf_pick(mom_k, mom_t, best=0.0, lncap=0.0, G=G, q=1)
    assert len(picks) == 1
    assert np.allclose(picks[0][0], G[n // 2]) and picks[0][1] == 8
