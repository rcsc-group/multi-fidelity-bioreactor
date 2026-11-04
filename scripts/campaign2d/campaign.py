"""2D (rpm, theta) multi-fidelity campaign driver. DRAFT, not approved. Nothing here submits by default.

Subcommands
  plan                      print the design, budget and wall-clock estimate (no compute, no sacct)
  submit-stage L8|L9|L10    print the sbatch commands only (default). Real submission needs --go.
  status                    coverage of the campaign domain and state of recorded campaign jobs (sacct)
  fit                       multi-fidelity fit + 2*sigma/mean map per QoI. Needs a compute node (SLURM_JOB_ID).
  next                      acquisition from the last fit map: next L8-site points for L9/L10
  ledger                    core-h per campaign job (sacct), running total, Kim counterfactual

Rules encoded here
  * mbessa-condo is never used. node2356 is excluded in every sbatch.
  * Warm L9/L10 starts need a settled L8 checkpoint at the SAME (rpm, theta): candidates are L8 sites only.
  * dtmix_0.50 and kLa_1T_25 / kLa_1T_50 at L10 are NOT validated for warm starts (diary 2026-10-24 W9):
    the L10 scalar points for those QoIs are flagged 'cold' in the plan and cost ~3x more.
  * Replay (looped period), amended 2026-11-04 (lead; plan section 3.1): the source of kLa_1T_10, kLa_1T_25 and
    dtmix_0.50 at L8 and L9 for rpm < 35. kLa_1T_50, tau_95, edr_95 and rpm >= 35 use the coupled DNS tail.
    Each replay run passes a flow check (tau_95 and posY_max within REPLAY_FLOW_TOL of its recording DNS).

Usage: uv run python scripts/campaign2d/campaign.py <subcommand> [options]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
STATE = Path(__file__).resolve().parent / "state"
JOBS = STATE / "jobs.jsonl"                       # written ONLY by a real (--go) submission

RPM = (15.0, 37.5)                                # campaign domain = Kim's range (lead, 2026-11-04)
THETA = (2.0, 7.0)
EXCLUDE = "node2356"                              # hard rule: every sbatch
FORBIDDEN_PARTITION = "mbessa-condo"              # hard rule: never planned or submitted
KIM_CORE_H_PER_CASE = 120.0                       # Kim 2025 p.6. PENDING AUDIT (another agent). Do not cite as final.
SEED = 20261104
# Persistent binary (scratch copy of lean-f1c11e0 was purged; kmix params still point there). Pre-flight: see plan.
BINARY_L8 = "/oscar/data/dharri15/eaguerov/bin/BioReactor-mpi-lean-2ef5b07"   # amended 2026-11-04: edr_95 build
# Partitions (sacctmgr 2026-11-04, user eaguerov): batch under account default, batch2 under account batch2.
# Value = (account, cpu cap for this user's QOS). mbessa-condo is never listed here (hard rule).
PARTITIONS = {"batch": ("default", 312), "batch2": ("batch2", 512)}

# Central estimates (core-h per point) with [low, high]. Sources in docs/campaign2d_plan.md section 4.
COST = {
    "L8_cold":        (55.0, 45.0, 70.0),     # sacct 7156006: 3h13m x 16 ranks = 52; diary kmix ~57
    "L9_warm":        (220.0, 190.0, 260.0),  # sacct 7847552 / 7847553: 6.0 h, 7.0 h x 32 ranks
    "L9_cold":        (580.0, 500.0, 650.0),  # diary 2026-10-19: fig9_l9_rpm32.5, 18.2 h x 32 ranks
    "L10_full":       (1500.0, 900.0, 3200.0),  # sacct 7512023 1500, 7512031 1162, 7720889 1500; 17.5 rpm ~3200
    "L10_hydro":      (500.0, 250.0, 960.0),  # 12 warm + 3 window cycles x 17-64 core-h/cycle (l10_cost_per_cycle)
    "L10_replay":     (927.0, 900.0, 950.0),  # sacct 7893240: 1d 4h59m x 32 ranks (22 replay cycles; NOT cheaper)
    "edr95_backfill_L10": (105.0, 60.0, 190.0),  # 3 window cycles from the end checkpoint
    # Replay point at L8 (plan 3.1): DNS to cycle 81 (sacct 8184447 = 29, 8184458 = 11) plus a 150-cycle replay
    # tail at 1/2.4 of the cold per-cycle cost (55/230 per cycle; 15 core-h). Estimate, not measured.
    "L8_replay":      (35.0, 25.0, 45.0),
}

# Replay rule, amended 2026-11-04 (lead; plan section 3.1). Sources of the scalar QoIs.
REPLAY_LEVELS = (8, 9)
REPLAY_RPM_MAX = 35.0                             # replay only for rpm < 35
REPLAY_QOIS = frozenset({"kLa_1T_10", "kLa_1T_25", "dtmix_0.50"})
COUPLED_ONLY_QOIS = frozenset({"kLa_1T_50", "tau_95", "edr_95"})   # always from the coupled DNS tail
REPLAY_FLOW_TOL = 0.02                            # amended from the pre-registered 1 % (L9 offset 1.0-1.7 %)
REPLAY_SPEEDUP = 2.4                              # replay tail vs coupled tail (diary 2026-10-31)

STAGES = {
    "L8": dict(ranks=16, walltime="16:00:00", n_new=24, protocol="cold, 80-cycle spin-up, release, 150 cycles"),
    "L9": dict(ranks=32, walltime="40:00:00", n_new=12, protocol="warm from L8 site (12 cycles), release, 40 cycles"),
    "L10": dict(ranks=32, walltime="47:00:00", n_new=8, protocol="warm from L8/L9 site (12 cycles), release, 40 cycles; chained"),
}

QOIS_CRITERION = ["tau_95", "edr_95", "kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "dtmix_0.50"]
EDR95_PRESENT = False                             # checked 2026-11-04: NOT saved per run (see plan section 1)


# ----------------------------------------------------------------- design
def domain_grid(n_rpm: int = 31, n_theta: int = 11) -> np.ndarray:
    r = np.linspace(*RPM, n_rpm)
    t = np.linspace(*THETA, n_theta)
    return np.array([(a, b) for a in r for b in t])


def lhs_maximin(n: int, seed: int = SEED, tries: int = 300) -> np.ndarray:
    """Maximin Latin hypercube on the campaign domain. Deterministic for a given seed."""
    rng = np.random.default_rng(seed)
    best, best_d = None, -1.0
    for _ in range(tries):
        u = np.column_stack([(rng.permutation(n) + rng.random(n)) / n for _ in range(2)])
        d = np.min([np.linalg.norm(u[i] - u[j]) for i in range(n) for j in range(i)])
        if d > best_d:
            best, best_d = u, d
    rpm = RPM[0] + best[:, 0] * (RPM[1] - RPM[0])
    th = THETA[0] + best[:, 1] * (THETA[1] - THETA[0])
    return np.column_stack([rpm, th])


def l8_sites() -> list[tuple[float, float]]:
    """(rpm, theta) of every finished L8 run in runs/ (params.json fidelity == 8, results.json present)."""
    out = []
    for d in sorted((ROOT / "runs").glob("*")):
        pf, rf = d / "params.json", d / "results.json"
        if not (pf.exists() and rf.exists()):
            continue
        try:
            p = json.loads(pf.read_text())
        except json.JSONDecodeError:
            continue
        if p.get("fidelity") != 8:
            continue
        rpm = round(float(p["omega_b"]) * 60.0 / (2 * math.pi), 3)
        out.append((rpm, float(p["theta_max"][0])))
    return sorted(set(out))


def run_id(level: int, rpm: float, theta: float, tag: str) -> str:
    return f"{tag}_l{level}_rpm{rpm:.2f}_th{theta:.2f}"


def sbatch_argv(level: str, rid: str, rpm: float, theta: float, partition: str = "batch",
                account: str = "default") -> tuple[list[str], str]:
    st = STAGES[level]
    argv = ["sbatch", "--no-requeue", f"--partition={partition}", f"--account={account}", "--nodes=1",
            f"--time={st['walltime']}",   # no --mem-per-cpu: srun --mem=0 in the template conflicts with it (smoke 8340540)
            "--cpus-per-task=1", f"--ntasks={st['ranks']}", f"--exclude={EXCLUDE}",
            f"--job-name=c2d-{level.lower()}-{rid}"]
    return argv, f"# run_id={rid} rpm={rpm:g} theta={theta:g} protocol: {st['protocol']}"


def assign_partitions(ranks: list[int], mode: str = "split") -> list[tuple[str, str]]:
    """(partition, account) per job. mode batch|batch2 puts every job on one partition. mode split fills batch
    up to its 312-cpu cap, then batch2 up to its 512-cpu cap (plan 2). Jobs past both caps are queued on the
    partition with the lower fill, and the dry run says so."""
    if mode in PARTITIONS:
        return [(mode, PARTITIONS[mode][0])] * len(ranks)
    used = {p: 0 for p in PARTITIONS}
    out = []
    for r in ranks:
        choice = next((p for p in ("batch", "batch2") if used[p] + r <= PARTITIONS[p][1]), None)
        if choice is None:
            choice = min(PARTITIONS, key=lambda p: used[p] / PARTITIONS[p][1])   # over both caps: waits
        used[choice] += r
        out.append((choice, PARTITIONS[choice][0]))
    return out


def plan_parts() -> dict[str, float]:
    """Central core-h of the planned campaign, by item. Used by `plan` and `ledger`."""
    return {
        "L8 new (24, cold)": STAGES["L8"]["n_new"] * COST["L8_cold"][0],
        "L9 new (12, warm)": STAGES["L9"]["n_new"] * COST["L9_warm"][0],
        "L9 cold surcharge (4)": 4 * (COST["L9_cold"][0] - COST["L9_warm"][0]),
        "L10 full (8)": STAGES["L10"]["n_new"] * COST["L10_full"][0],
        "L10 hydro-only (6)": 6 * COST["L10_hydro"][0],
        "edr_95 backfill": 20 * COST["edr95_backfill_L10"][0] + 120 * 0.75,
    }


def scalar_source(level: int, rpm: float, qoi: str) -> str:
    """'replay' or 'coupled' for one QoI at one point (plan section 3.1, amended 2026-11-04)."""
    if level in REPLAY_LEVELS and rpm < REPLAY_RPM_MAX and qoi in REPLAY_QOIS:
        return "replay"
    return "coupled"


def point_needs_coupled(level: int, rpm: float, qois) -> bool:
    """True if any QoI of the point must come from the coupled DNS tail (so the replay saves nothing)."""
    return any(scalar_source(level, rpm, q) == "coupled" for q in qois)


def replay_saving(pts, level: int = 8, qois=None) -> float:
    """Core-h saved by the replay for the points `pts` (rpm, theta). Zero for any point that still needs the
    coupled tail. With the criterion QoIs (kLa_1T_50 is reported at every point) it is zero for every point."""
    qois = list(QOIS_CRITERION if qois is None else qois)
    saved = 0.0
    for rpm, _th in pts:
        if not point_needs_coupled(level, rpm, qois):
            saved += COST[f"L{level}_cold"][0] - COST[f"L{level}_replay"][0]
    return saved


def flow_check_ok(dev_tau95: float, dev_posY: float) -> bool:
    """Replay flow check: relative deviation of tau_95 and posY_max from the recording DNS, both within REPLAY_FLOW_TOL."""
    return abs(dev_tau95) <= REPLAY_FLOW_TOL and abs(dev_posY) <= REPLAY_FLOW_TOL


# ----------------------------------------------------------------- plan
def cmd_plan(_args) -> None:
    pts = lhs_maximin(STAGES["L8"]["n_new"])
    print("DOMAIN  rpm", RPM, " theta", THETA, "deg  (Kim 2025 range; theta 9 data exist but are outside the domain)")
    print(f"L8 existing sites inside domain: {sum(1 for r, t in l8_sites() if RPM[0] <= r <= RPM[1] and THETA[0] <= t <= THETA[1])}")
    print(f"L8 new LHS points (maximin, seed {SEED}): {len(pts)}")
    for rpm, th in pts:
        print(f"   rpm {rpm:6.2f}  theta {th:4.2f}")
    print("\nCOST MODEL (core-h per point, central [low, high]):")
    for k, (c, lo, hi) in COST.items():
        print(f"   {k:20s} {c:8.0f}  [{lo:.0f}, {hi:.0f}]")
    parts = plan_parts()
    tot = sum(parts.values())
    print("\nCENTRAL ESTIMATE (core-h):")
    for k, v in parts.items():
        print(f"   {k:26s} {v:8.0f}")
    print(f"   {'TOTAL':26s} {tot:8.0f}")
    print(f"WALL at 288 cores (priority QOS 312 cap; sacct queue median 0.2 h, p90 5.6 h): "
          f"~{tot / 288:.0f} h of full-cluster time; critical path 7-10 d with three adaptive rounds (plan section 2)")
    print(f"\nedr_95 saved per run: {EDR95_PRESENT}  (stand-in in fit: edr_mean_t, NOT the 95th percentile)")
    sav = replay_saving(pts, level=8)
    print(f"Replay saving, L8 new points (plan 3.1; kLa_1T_50 keeps the coupled tail): {sav:.0f} core-h")
    print(f"Replay flow check tolerance: {REPLAY_FLOW_TOL:.0%} (amended from 1 %, plan 3.1)")
    print(f"Kim counterfactual: total / {KIM_CORE_H_PER_CASE:g} core-h per case (PENDING AUDIT)")


# ----------------------------------------------------------------- submit
def cmd_submit(args) -> None:
    level = args.level
    if level not in STAGES:
        raise SystemExit(f"unknown level {level}")
    if args.go and args.dry_run:
        raise SystemExit("--go and --dry-run are exclusive")
    if not args.go and not args.dry_run:
        print("No --go: printing only (dry run).")
    if level == "L8":
        pts = lhs_maximin(STAGES["L8"]["n_new"])
        todo = [(float(r), float(t)) for r, t in pts]
    else:
        nxt = STATE / f"next_{level}.json"
        if not nxt.exists():
            raise SystemExit(f"no acquisition batch yet: run `next` after `fit` ({nxt} missing)")
        todo = [tuple(map(float, p)) for p in json.loads(nxt.read_text())["points"]]
    ranks = STAGES[level]["ranks"]
    placement = assign_partitions([ranks] * len(todo), args.partition)
    for (rpm, th), (part, acct) in zip(todo, placement):
        rid = run_id(int(level[1:]), rpm, th, tag="c2d")
        argv, note = sbatch_argv(level, rid, rpm, th, partition=part, account=acct)
        print(" ".join(argv), "#", note, f"binary={BINARY_L8}")
        if args.go:
            if level != "L8":
                raise SystemExit("--go is wired for L8 only; L9/L10 need the warm builder reviewed first")
            _real_submit_l8(rid, rpm, th, part, acct)
    for p in PARTITIONS:
        n = sum(1 for q, _ in placement if q == p)
        print(f"# partition {p}: {n} jobs x {ranks} ranks = {n * ranks} cpu concurrently "
              f"(cap {PARTITIONS[p][1]} cpu, plan 2)")
    if not args.go:
        print(f"(dry run: {len(todo)} commands printed; nothing submitted)")


def _real_submit_l8(rid: str, rpm: float, th: float, partition: str, account: str) -> None:
    """Real submission path. Reuses simulate.submit_slurm (partition and account passed through)."""
    from scripts import simulate
    p = json.loads((ROOT / "runs/kmix_l8_rpm32.5/params.json").read_text())
    p.update(run_id=rid, omega_b=2 * math.pi * rpm / 60.0, theta_max=[th, 0.0, 0.0], _binary=BINARY_L8)
    jid = simulate.submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs",
                                walltime=STAGES["L8"]["walltime"], ntasks=STAGES["L8"]["ranks"],
                                exclude=EXCLUDE, partition=partition, account=account)
    STATE.mkdir(parents=True, exist_ok=True)
    with JOBS.open("a") as f:
        f.write(json.dumps(dict(job=str(jid), level="L8", run=rid, rpm=rpm, theta=th,
                                ranks=STAGES["L8"]["ranks"], partition=partition, account=account,
                                submitted=datetime.now(timezone.utc).isoformat())) + "\n")


# ----------------------------------------------------------------- status / ledger
def _sacct(ids: list[str]) -> dict[str, tuple[str, str, str]]:
    if not ids:
        return {}
    out = subprocess.run(["sacct", "-X", "-n", "-P", "-j", ",".join(ids), "-o", "JobID,Elapsed,AllocCPUS,State"],
                         capture_output=True, text=True).stdout
    res = {}
    for line in out.splitlines():
        jid, el, cpus, st = (line.split("|") + ["", "", "", ""])[:4]
        res[jid] = (el, cpus, st)
    return res


def _hours(s: str) -> float:
    days = 0
    if "-" in s:
        d, s = s.split("-", 1)
        days = int(d)
    h, m, sec = (s.split(":") + ["0", "0"])[:3]
    return days * 24 + int(h) + int(m) / 60 + float(sec) / 3600


# ----------------------------------------------------------------- warm-start settled gate (2026-11-04)
# Every warm-started run (params t_checkpoint > 0) is gated. Plan section 1.4 / 3. Status is PASS, FAILED-GATE
# (excluded from the fit) or PENDING (data not yet sufficient). The gate is evaluated by `status`.
WARM_CYCLES = 12                       # warm-up cycles after the source checkpoint (plan 1.4)
GATE_DRIFT_TOL = 0.02                  # max relative change per cycle, tau_95 and EDR, last 3 complete cycles
GATE_HEIGHT_TOL = 0.01                 # free-surface height range at release vs the source run
GATE_FS_LINE = "fs sanity [post-solid]: OK"   # must be in the run's .err (BioReactor3D CLAUDE.md, restart rule)
GATES_JSON = STATE / "gates.json"


def _rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b else math.inf


def settled_gate(t, tau95, edr, period, warm_end, release_t, height_range, src_height_range, err_text):
    """Settled-state gate for one warm-started run. Returns (status, detail).

    Drift: cycle means of tau_95 and EDR over the last 3 complete cycles after warm_end. The two cycle-to-cycle
    relative changes must each be below GATE_DRIFT_TOL. Height: free-surface range at release vs the source run,
    relative change below GATE_HEIGHT_TOL. fs sanity: GATE_FS_LINE must be in err_text.
    None for a check means its data are not there yet (PENDING). Any False gives FAILED-GATE.
    """
    t = np.asarray(t, float)
    tau95 = np.asarray(tau95, float)
    edr = np.asarray(edr, float)
    det = {"drift_ok": None, "height_ok": None, "fs_ok": None}
    post = t >= warm_end
    if post.any():
        n_done = int(math.floor((t[post].max() - warm_end) / period + 1e-9))   # complete cycles after warm-up
    else:
        n_done = 0
    if n_done >= 3:
        ks = range(n_done - 3, n_done)
        means = {"tau95": [], "edr": []}
        for k in ks:
            m = post & (t >= warm_end + k * period) & (t < warm_end + (k + 1) * period)
            means["tau95"].append(float(tau95[m].mean()))
            means["edr"].append(float(edr[m].mean()))
        drifts = {}
        for name, vals in means.items():
            drifts[name] = [_rel(vals[i + 1], vals[i]) for i in range(len(vals) - 1)]
        worst = max(max(d) for d in drifts.values())
        det.update(drift_ok=bool(worst < GATE_DRIFT_TOL), drift_max=worst, drift=drifts)
    if release_t is not None and height_range is not None and src_height_range is not None:
        dh = _rel(height_range, src_height_range)
        det.update(height_ok=bool(dh < GATE_HEIGHT_TOL), height_rel=dh)
    if err_text is not None:
        det["fs_ok"] = bool(GATE_FS_LINE in err_text)
    vals = [det["drift_ok"], det["height_ok"], det["fs_ok"]]
    if any(v is False for v in vals):
        status = "FAILED-GATE"
    elif any(v is None for v in vals):
        status = "PENDING"
    else:
        status = "PASS"
    return status, det


def warm_run_gate(run_dir: Path) -> tuple[str, dict]:
    """Read one warm-started run and its source run (params _source_run, a runs/ dir name) and gate it."""
    from scripts.postprocess import _load_dat, _t_scales   # reuse the repo's readers and nondim scales
    params = json.loads((run_dir / "params.json").read_text())
    _, period = _t_scales(params)
    t_ck = float(params["t_checkpoint"])
    warm_end = t_ck + WARM_CYCLES * period
    release_t = t_ck + period * float(params.get("n_mix_cycles", 0))
    sh = _load_dat(run_dir / "shear_stress.dat")
    t, tau95 = sh[:, 1], sh[:, 2]
    edr = sh[:, 9] if sh.shape[1] >= 10 else np.full_like(tau95, math.nan)
    height_range = src_height = None
    vol = run_dir / "vol_frac_interf.dat"
    if vol.exists():
        v = _load_dat(vol)
        h = v[:, 4] - v[:, 5]                              # posY_max - posY_min
        m = (v[:, 1] >= release_t) & (v[:, 1] < release_t + period)
        if m.any():
            height_range = float(h[m].mean())
    src_name = params.get("_source_run")
    if src_name:
        src = ROOT / "runs" / src_name
        sv = _load_dat(src / "vol_frac_interf.dat")
        tend = float(sv[:, 1].max())
        ms = sv[:, 1] >= tend - period
        src_height = float((sv[ms, 4] - sv[ms, 5]).mean())
    err_text = None
    if params.get("_err_file") and Path(params["_err_file"]).exists():
        err_text = Path(params["_err_file"]).read_text(errors="replace")
    # `release_t` is passed only when the run has reached it, so the release check can be PENDING.
    reached = bool(t.max() >= release_t)
    return settled_gate(t, tau95, edr, period, warm_end, release_t if reached else None,
                        height_range, src_height, err_text)


def cmd_status(_args) -> None:
    sites = l8_sites()
    inside = [(r, t) for r, t in sites if RPM[0] <= r <= RPM[1] and THETA[0] <= t <= THETA[1]]
    print(f"L8 finished sites in domain: {len(inside)} (all levels: {len(sites)})")
    jobs = [json.loads(x) for x in JOBS.read_text().splitlines()] if JOBS.exists() else []
    print(f"campaign jobs recorded: {len(jobs)}")
    if jobs:
        st = _sacct([j["job"] for j in jobs])
        for j in jobs:
            print(f"   {j['job']:>10s} {j['run']:32s} {st.get(j['job'], ('?', '?', 'not in sacct'))[2]}")
    gates = {}
    for d in sorted((ROOT / "runs").glob("*")):
        pf = d / "params.json"
        if not pf.exists():
            continue
        try:
            if float(json.loads(pf.read_text()).get("t_checkpoint", 0.0)) <= 0.0:
                continue
            status, det = warm_run_gate(d)
        except (KeyError, IndexError, ValueError, FileNotFoundError, json.JSONDecodeError):
            status, det = "PENDING", {"reason": "data not readable yet"}
        gates[d.name] = status
        print(f"GATE {d.name:24s} {status}")
    STATE.mkdir(parents=True, exist_ok=True)
    GATES_JSON.write_text(json.dumps(gates, indent=1))
    print(f"warm-started runs gated: {len(gates)}; FAILED-GATE (excluded from fit): "
          f"{sum(v == 'FAILED-GATE' for v in gates.values())}")


def cmd_ledger(_args) -> None:
    jobs = [json.loads(x) for x in JOBS.read_text().splitlines()] if JOBS.exists() else []
    st = _sacct([j["job"] for j in jobs])
    total = 0.0
    per = {}
    for j in jobs:
        if j["job"] not in st:
            continue
        el, cpus, state = st[j["job"]]
        ch = _hours(el) * int(cpus)
        per.setdefault(j["level"], [0.0, 0])
        per[j["level"]][0] += ch
        per[j["level"]][1] += 1
        total += ch
    for lev, (ch, n) in sorted(per.items()):
        print(f"{lev}: {ch:10.1f} core-h over {n} jobs")
    print(f"TOTAL campaign core-h (sacct Elapsed x AllocCPUS): {total:.1f}")
    print(f"Kim counterfactual: {total / KIM_CORE_H_PER_CASE:.1f} Kim cases at {KIM_CORE_H_PER_CASE:g} core-h "
          "(PENDING AUDIT of the 120 core-h figure)")
    print(f"No budget cap (lead, 2026-11-04). Plan central estimate: {sum(plan_parts().values()):.0f} core-h.")


# ----------------------------------------------------------------- fit / next
def fit(level: str) -> None:
    """Multi-fidelity fit (reuses FullBayesKRRLRGP and CachedLooKRR). Refuses to run outside a compute job."""
    if "SLURM_JOB_ID" not in os.environ:
        raise SystemExit("fit is heavy: run it inside a compute job (sbatch ... --wrap 'uv run python ...'), "
                         "never on the login node")
    from scripts.test_t3_learning_curve import CachedLooKRR  # noqa: F401  (LF KRR, reused)
    from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP
    domain = np.array([list(RPM), list(THETA)])
    raise NotImplementedError("wire: collect ln(QoI) per site (results.json + hydro_dataset_v2 row builder), "
                              "LF = CachedLooKRR(design_space=domain) on L8, HF = FullBayesKRRLRGP(domain, lf_model=lf) "
                              "on the level's sites; 2*sigma/mean = 2*sqrt(exp(v)-1) with v the ln-variance of the "
                              "mixture (w.(sc^2 nu/(nu-2)) + Var_w(mu)). Written only after the user approves the plan.")


def cmd_fit(args) -> None:
    fit(args.level)


def cmd_next(args) -> None:
    raise SystemExit("next needs the fit map (fit output). Not run in this draft. Rule: argmax over L8 sites of "
                     "max_q(2sigma/mean) divided by the cost of the cheapest level that informs it; batch 2 at L10, "
                     "4 at L9.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan").set_defaults(func=cmd_plan)
    s = sub.add_parser("submit-stage")
    s.add_argument("level", choices=["L8", "L9", "L10"])
    s.add_argument("--dry-run", action="store_true")
    s.add_argument("--partition", choices=["batch", "batch2", "split"], default="split",
                   help="batch (account default, 312 cpu), batch2 (account batch2, 512 cpu), or split (fill batch, then batch2)")
    s.add_argument("--go", action="store_true", help="real submission; explicit approval only")
    s.set_defaults(func=cmd_submit)
    sub.add_parser("status").set_defaults(func=cmd_status)
    f = sub.add_parser("fit")
    f.add_argument("level", choices=["L8", "L9", "L10"])
    f.set_defaults(func=cmd_fit)
    sub.add_parser("next").set_defaults(func=cmd_next)
    sub.add_parser("ledger").set_defaults(func=cmd_ledger)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
