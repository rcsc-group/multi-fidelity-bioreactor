"""BO loop driver for the DFD comparison (2026-11-07). State: scripts/bo2d/state/tasks.json.

Arms (run-id tag):
  l9only (bo1): per point, a cold L8 spin-up (start state only, never used as data), then L9 warm from it.
  chain  (bo3): per point, the lowest level not yet run there: L7 cold, then L8 warm from it, then L9 warm from it.
                A point is promoted to the next level only if it is still in the arm's top-q proposals.
Shared seed: SEED9 L9 rows (both arms); the chain also has the real L8 rows and the bo_l7 seed rows.

  uv run python scripts/bo2d/loop.py status
  uv run python scripts/bo2d/loop.py propose ARM [--go]     fit the arm, pick q points, create/advance tasks
  uv run python scripts/bo2d/loop.py advance [--go]         submit the next step of every task whose step finished
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.bo2d import bo, submit as S  # noqa: E402
C = bo.C

TASKS = S.STATE / "tasks.json"
CFG = S.STATE / "config.json"          # frozen: cap, q
TAGS = {"l9only": "bo1", "chain": "bo4"}   # bo3 = chain with unfair archive seed (abandoned 2026-11-08)
SEED9 = ["c2d_l9_rpm16.35_th2.33", "c2d_l9_rpm19.14_th2.89", "c2d_l9_rpm23.14_th2.13", "c2d_l9_rpm36.35_th5.24"]
QK, QT = "kLa_1T_10", "tau_95"
Q_LF = 8      # chain: extra L7+L8 points per round at the next-best candidates (synthetic test 2026-11-07)


def load_tasks():
    return json.loads(TASKS.read_text()) if TASKS.exists() else []


def save_tasks(t):
    S.STATE.mkdir(parents=True, exist_ok=True)
    TASKS.write_text(json.dumps(t, indent=1))


def rows(qoi, prefixes, levels):
    """Fit rows from runs/<prefix>* (finished, gate not FAILED) plus SEED9."""
    src, col = C.QOI_SOURCE[qoi]
    out = []
    for pre in prefixes:
        for p in sorted((C.ROOT / "runs").glob(pre + "*")):
            r = C.campaign_row(p, col, src)
            if r is None or r["level"] not in levels:
                continue
            pp = json.loads((p / "params.json").read_text())
            if float(pp.get("t_checkpoint", 0) or 0) > 0 and C.warm_run_gate(p)[0] == "FAILED-GATE":
                continue
            out.append(r)
    d = pd.DataFrame(out, columns=["level", "rpm", "theta", "run", "y"])
    bad = ~np.isfinite(d.y.astype(float)) | (d.y.astype(float) <= 0)
    if bad.any():   # 9/36 bo_l7 seeds have non-finite kLa (2026-11-08); a NaN row makes the whole chain fit NaN
        print(f"  dropping {int(bad.sum())} non-finite rows: {list(d.run[bad])}", file=sys.stderr)
    return d[~bad].reset_index(drop=True)


def arm_data(arm, qoi):
    seed = rows(qoi, SEED9, {9})
    own = rows(qoi, [TAGS[arm] + "_"], {9})
    d9 = pd.concat([seed, own], ignore_index=True)
    if arm == "l9only":
        return {9: d9}
    # fair seed (user, 2026-11-08): the chain sees L7/L8 only at the 4 shared SEED9 sites, plus its own picks
    real8 = C.load_db(qoi); real8 = _at_seed(real8[real8.level == 8], seed)
    d8 = pd.concat([real8, rows(qoi, [TAGS[arm] + "_l8_"], {8})], ignore_index=True)
    # plus the L7 seed runs off the regular grid lines (theta 4 / 5.5 / 7 rows look hand-placed); charged L7_SEED_COST
    s7 = rows(qoi, ["bo_l7_"], {7})
    d7 = pd.concat([_at_seed(s7, seed), _off_grid(s7), rows(qoi, [TAGS[arm] + "_l7_"], {7})],
                   ignore_index=True).drop_duplicates(subset=["rpm", "theta", "run"])
    return {7: d7, 8: d8, 9: d9}


GRID_THETA = (4.0, 5.5, 7.0)


def _off_grid(d, tol=0.02):
    on = np.zeros(len(d), bool)
    for t_ in GRID_THETA:
        on |= ((d.theta - t_).abs() < tol).to_numpy()
    return d[~on].reset_index(drop=True)


def l7_seed_cost():
    """core-h charged to the chain for the reused off-grid L7 seed runs (from sacct, 5 core-h median)."""
    return 5.0 * len(_off_grid(rows(QT, ["bo_l7_"], {7})))


def _at_seed(d, seed, tol=0.01):
    keep = np.zeros(len(d), bool)
    for r_, t_ in zip(seed.rpm, seed.theta):
        keep |= ((d.rpm - r_).abs() < tol).to_numpy() & ((d.theta - t_).abs() < tol).to_numpy()
    return d[keep].drop_duplicates(subset=["rpm", "theta"]).reset_index(drop=True)


# measured 2026-11-08 (measure_costs.py): L7 5, L8 spin 21, L9 warm 350 core-h; L8 warm taken = L8 spin (70 vs 80 cycles)
COST = {"l9only": {9: 371.0}, "chain": {7: 5.0, 8: 26.0, 9: 376.0}}


def fit_arm(arm):
    """Locked method (tag bo-lock-267421f): l9only = single-level aEI path; chain = theta-basis KOH + MF level rule."""
    K, T = arm_data(arm, QK), arm_data(arm, QT)
    if arm == "l9only":
        mk, mt = bo.fit_single(K[9], QK), bo.fit_single(T[9], QT)
        return (lambda X: bo.single_moments(mk, X)), (lambda X: bo.single_moments(mt, X)), K, T
    bo.CHAIN_THETA = True
    mk, mt = bo.fit_chain3(K[7], K[8], K[9], QK), bo.fit_chain3(T[7], T[8], T[9], QT)
    return (lambda X: bo.chain_moments(mk, X)), (lambda X: bo.chain_moments(mt, X)), K, T


def _has(d, x, tol=0.01):
    return bool(((d.rpm - x[0]).abs() < tol).mul((d.theta - x[1]).abs() < tol).any()) if len(d) else False


def job_state(j):
    out = subprocess.run(["sacct", "-X", "-n", "-P", "-j", str(j), "-o", "State"], capture_output=True, text=True).stdout
    return (out.split() or ["?"])[0]


def propose(arm, go):
    """Asynchronous (2026-11-09, user): keep q tasks in flight per arm; fill the free slots as soon as tasks end.
    Same acquisition; picks within R_EXCL of an in-flight point are skipped (standard async exclusion)."""
    cfg = json.loads(CFG.read_text())
    lncap, q = np.log(cfg["cap"]), int(cfg["q"])
    tasks = load_tasks()
    busy = [t for t in tasks if t["arm"] == arm and t["state"] in ("running", "new")]
    free = q - len(busy)
    if free <= 0:
        print(f"{arm}: {len(busy)} task(s) in flight; no free slot"); return
    pend = np.array([t["x"] for t in busy], float).reshape(-1, 2)
    Mk, Mt, K, T = fit_arm(arm)
    G = bo.grid()
    Xe = pd.concat([K[l][["rpm", "theta"]] for l in K]).drop_duplicates().to_numpy(float)
    best = bo.effective_best(Mk(Xe)[0], Mt(Xe)[0], lncap)
    if best is None:
        best = float(np.log(K[9].y).min())
    mk, ck = Mk(G); mt, _ = Mt(G)
    qq = q + len(pend)
    if arm == "chain":
        picks = bo.mf_pick(mk, mt, best, lncap, G, qq)
    else:
        picks = bo.pick_budget(bo.aei(mk, ck, mt, best, lncap, COST[arm]), G, qq * COST[arm][9], COST[arm], qq)
    if len(pend):
        far = lambda x: np.min(np.linalg.norm(bo.scaled(pend) - bo.scaled(np.array([x])), axis=1)) >= bo.R_EXCL
        picks = [(x, l) for x, l in picks if far(x)]
    picks = picks[:free]
    rec = bo.recommend_m(Mk, Mt, K[9], lncap)
    print(f"{arm}: n9={len(K[9])} best={best:.3f} recommend={rec}")
    for x, lev in picks:
        x = [round(float(x[0]), 3), round(float(x[1]), 3)]
        if arm == "l9only":
            steps = ["spin8", "l9"]
        else:   # cascade: every level below lev must exist at x (warm source, same forcing); skip those present
            # a run that exists but has no finite kLa (L7 at rpm <= ~21: 10 % crossing within half a period of
            # release, 2026-11-08) is still a valid warm source and cannot give data at that level: step up.
            while lev < 9 and _nearest(lev, x) is not None:
                lev += 1
            steps = [f"l{l}" for l in range(7, lev + 1) if _nearest(l, x) is None and not _has(K[l], x)]
            if not steps:
                continue
        t = dict(arm=arm, x=x, steps=steps, i=0, jobs=[], runs=[], state="new")
        print(f"  propose {x} level {lev} steps {steps}")
        tasks.append(t)
    if go:
        save_tasks(tasks)
        advance(go)


def _source_for(t, step):
    """Run to warm-start this step from (same rpm, theta)."""
    x = t["x"]
    if step == "l9":
        if t["arm"] == "l9only":
            return C.run_id(8, *x, tag=TAGS["l9only"])
        cand = [C.run_id(8, *x, tag=TAGS["chain"])]
        return next((c for c in cand if (C.ROOT / "runs" / c / "shear_stress.dat").exists()), None) or _nearest(8, x)
    if step == "l8":
        return C.run_id(7, *x, tag=TAGS["chain"]) if (C.ROOT / "runs" / C.run_id(7, *x, tag=TAGS["chain"])).exists() \
            else _nearest(7, x)
    return None


def _nearest(level, x):
    """A finished run at the same (rpm, theta) (within 0.01) and this level, any name, with a checkpoint, else None."""
    d = pd.concat([C.load_db(QT), rows(QT, ["bo_", "bo4_"], {7, 8})], ignore_index=True)
    d = d[(d.level == level) & ((d.rpm - x[0]).abs() < 0.01) & ((d.theta - x[1]).abs() < 0.01)]
    for r in d.run:
        if (C.ROOT / "runs" / r / "checkpoint.dump").exists():
            return r
    return None


def advance(go):
    tasks = load_tasks()
    for t in tasks:
        if t["state"] in ("done", "failed"):
            continue
        if t["state"] == "running":
            st = job_state(t["jobs"][-1])
            if st in ("PENDING", "RUNNING", "?"):
                continue
            if st != "COMPLETED":
                t["state"] = "failed"; print(f"FAILED {t['runs'][-1]} ({st})")
                d = C.ROOT / "runs" / t["runs"][-1]       # move aside so a retry at the same x can submit
                if d.exists():
                    d.rename(d.with_name(f"{d.name}.failed-{t['jobs'][-1]}"))
                continue
            t["i"] += 1
            if t["i"] >= len(t["steps"]):
                t["state"] = "done"; print(f"done {t['arm']} {t['x']}"); continue
        step = t["steps"][t["i"]]
        x, tag = t["x"], TAGS[t["arm"]]
        if step == "spin8":
            jid = S.cold(8, *x, arm=t["arm"], go=go, tag=tag, spin_only=True); rid = C.run_id(8, *x, tag=tag)
        elif step == "l7":
            jid = S.cold(7, *x, arm=t["arm"], go=go, tag=tag); rid = C.run_id(7, *x, tag=tag)
        else:
            lev = int(step[1]); src = _source_for(t, step)
            if src is None:
                print(f"no source for {step} at {x}"); continue
            jid = S.warm(lev, src, *x, arm=t["arm"], go=go, tag=tag); rid = C.run_id(lev, *x, tag=tag)
        if jid:
            t["jobs"].append(jid); t["runs"].append(rid); t["state"] = "running"
            print(f"submitted {step} {rid} job {jid}")
    if go:
        save_tasks(tasks)


def status():
    for t in load_tasks():
        print(t["arm"], t["x"], t["steps"], "step", t["i"], t["state"], t["runs"][-1:] or "")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd"); ap.add_argument("arm", nargs="?"); ap.add_argument("--go", action="store_true")
    a = ap.parse_args()
    {"status": lambda: status(), "advance": lambda: advance(a.go), "propose": lambda: propose(a.arm, a.go)}[a.cmd]()
