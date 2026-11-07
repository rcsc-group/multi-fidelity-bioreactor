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
TAGS = {"l9only": "bo1", "chain": "bo3"}
SEED9 = ["c2d_l9_rpm16.35_th2.33", "c2d_l9_rpm19.14_th2.89", "c2d_l9_rpm23.14_th2.13", "c2d_l9_rpm36.35_th5.24"]
QK, QT = "kLa_1T_10", "tau_95"


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
    return pd.DataFrame(out, columns=["level", "rpm", "theta", "run", "y"])


def arm_data(arm, qoi):
    seed = rows(qoi, SEED9, {9})
    own = rows(qoi, [TAGS[arm] + "_"], {9})
    d9 = pd.concat([seed, own], ignore_index=True)
    if arm == "l9only":
        return {9: d9}
    real8 = C.load_db(qoi); real8 = real8[real8.level == 8]
    d8 = pd.concat([real8, rows(qoi, [TAGS[arm] + "_l8_"], {8})], ignore_index=True)
    d7 = rows(qoi, ["bo_l7_", TAGS[arm] + "_l7_"], {7})
    return {7: d7, 8: d8, 9: d9}


def fit_arm(arm):
    K, T = arm_data(arm, QK), arm_data(arm, QT)
    if arm == "l9only":
        return bo.fit_single(K[9], QK), bo.fit_single(T[9], QT), K, T
    return bo.fit_chain(K[7], K[8], K[9], QK)[0], bo.fit_chain(T[7], T[8], T[9], QT)[0], K, T


def job_state(j):
    out = subprocess.run(["sacct", "-X", "-n", "-P", "-j", str(j), "-o", "State"], capture_output=True, text=True).stdout
    return (out.split() or ["?"])[0]


def propose(arm, go):
    cfg = json.loads(CFG.read_text())
    lncap, q = np.log(cfg["cap"]), int(cfg["q"])
    mk, mt, K, T = fit_arm(arm)
    best, _ = bo.best_feasible(K[9], T[9], cfg["cap"])
    G = bo.grid()
    a = bo.cei(*bo.gauss(mk, G), *bo.gauss(mt, G), best, lncap)
    B = bo.pick_batch(a, G, q)
    tasks = load_tasks()
    busy = [t for t in tasks if t["arm"] == arm and t["state"] == "running"]
    if busy:
        print(f"{arm}: {len(busy)} task(s) still running; no new proposal"); return
    rec = bo.recommend(mk, mt, K[9], lncap)
    print(f"{arm}: n9={len(K[9])} best feasible obs={best} recommend={rec}")
    for x in B:
        x = [round(float(x[0]), 3), round(float(x[1]), 3)]
        if arm == "l9only":
            steps = ["spin8", "l9"]
        else:
            lev = bo.next_level(np.array([x]), K) or 9
            if lev > 7:   # snap to the existing lower-level row near x, so the warm source has the same forcing
                d = K[lev - 1]; P = d[["rpm", "theta"]].to_numpy(float)
                j = int(np.argmin(np.linalg.norm(bo.scaled(P) - bo.scaled(np.array([x])), axis=1)))
                x = [round(float(P[j, 0]), 3), round(float(P[j, 1]), 3)]
            steps = {7: ["l7", "l8"], 8: ["l8"], 9: ["l9"]}[lev]     # one cascade step; L7 and L8 together (cheap)
        t = dict(arm=arm, x=x, steps=steps, i=0, jobs=[], runs=[], state="new")
        print(f"  propose {x} steps {t['steps']}")
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
    d = pd.concat([C.load_db(QT), rows(QT, ["bo_", "bo3_"], {7, 8})], ignore_index=True)
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
                t["state"] = "failed"; print(f"FAILED {t['runs'][-1]} ({st})"); continue
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
