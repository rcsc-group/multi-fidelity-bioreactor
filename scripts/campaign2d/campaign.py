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
import re
import math
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from scripts import paths  # noqa: E402  (machine-specific folders)
STATE = Path(__file__).resolve().parent / "state"
JOBS = STATE / "jobs.jsonl"                       # written ONLY by a real (--go) submission

RPM = (15.0, 37.5)                                # campaign domain = Kim's range (lead, 2026-11-04)
THETA = (2.0, 7.0)
EXCLUDE = "node2356,node2410"                    # hard rule: every sbatch; node2410 added 2026-11-11 (load ~194, L8 ran 10x slow)
FORBIDDEN_PARTITION = "mbessa-condo"              # hard rule: never planned or submitted
KIM_CORE_H_PER_CASE = 120.0                       # Kim 2025 p.6. PENDING AUDIT (another agent). Do not cite as final.
SEED = 20261104
# Persistent binary (scratch copy of lean-f1c11e0 was purged; kmix params still point there). Pre-flight: see plan.
BINARY_L8 = f"{paths.DATA}/bin/BioReactor-mpi-lean-2ef5b07"
# L9 warm restarts: the same source (2ef5b07) built with -DCROSS_LEVEL_WARMSTART=1 (lead 2026-11-05).
BINARY_L9 = f"{paths.DATA}/bin/replay-gate/BioReactor-mpi-lean-2ef5b07-xlevel"   # amended 2026-11-04: edr_95 build
# Partitions (sacctmgr 2026-11-04, user eaguerov): batch under account default, batch2 under account batch2.
# Value = (account, cpu cap for this user's QOS). mbessa-condo is never listed here (hard rule).
PARTITIONS = {"batch": ("default", 312), "batch2": ("batch2", 512)}

# Central estimates (core-h per point) with [low, high]. Sources in internal/docs/campaign2d_plan.md section 4.
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
COUPLED_ONLY_QOIS = frozenset({"tau_95", "edr_95"})   # from the coupled DNS tail at every rpm
# kLa_1T_50 (lead, 2026-11-04): taken only from runs that have a coupled tail anyway (rpm >= 35, L10). It is
# reported, not gated, and it never forces a coupled tail elsewhere.
REPORTED_ONLY_QOIS = frozenset({"kLa_1T_50"})
REPLAY_FLOW_TOL = 0.02                            # amended from the pre-registered 1 % (L9 offset 1.0-1.7 %)
REPLAY_SPEEDUP = 2.4                              # replay tail vs coupled tail (diary 2026-10-31)

STAGES = {
    "L8": dict(ranks=16, walltime="16:00:00", n_new=24, protocol="cold, 80-cycle spin-up, release, 150 cycles"),
    "L9": dict(ranks=32, walltime="12:00:00", n_new=12, protocol="warm from L8 site (12 cycles), release, 40 cycles"),
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
    """'replay', 'coupled' or 'reported' for one QoI at one point (plan section 3.1, amended 2026-11-04).
    'reported' (kLa_1T_50) is taken from a coupled tail where one exists; it never forces one."""
    if qoi in REPORTED_ONLY_QOIS:
        return "reported"
    if level in REPLAY_LEVELS and rpm < REPLAY_RPM_MAX and qoi in REPLAY_QOIS:
        return "replay"
    return "coupled"


def point_needs_coupled(level: int, rpm: float, qois) -> bool:
    """True if any QoI of the point must come from the coupled DNS tail (so the replay saves nothing)."""
    return any(scalar_source(level, rpm, q) == "coupled" for q in qois)


def replay_saving(pts, level: int = 8, qois=None) -> float:
    """Core-h saved by the replay for the points `pts` (rpm, theta). Zero for any point that still needs the
    coupled tail. With the rule as decided (tau_95 and edr_95 coupled-only) it is zero for every point."""
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
    alt = replay_saving(pts, level=8, qois=[q for q in QOIS_CRITERION if q not in COUPLED_ONLY_QOIS])
    print(f"Replay saving, L8 new points, decided rule (tau_95, edr_95 coupled-only): {sav:.0f} core-h")
    print(f"   alternative, NOT decided (tau_95, edr_95 also replay-sourced, rpm < 35): {alt:.0f} core-h")
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
        data = json.loads(nxt.read_text())
        if not data.get("sources") or len(data["sources"]) != len(data["points"]):
            raise SystemExit(f"{nxt} has no warm source per point: cold L9/L10 refused (lead rule 2026-11-04)")
        todo = [tuple(map(float, p)) for p in data["points"]]
    ranks = STAGES[level]["ranks"]
    # Lead rule (2026-11-05): every MPI job runs on batch only (no batch2: UCX wireup failure on gpu4001).
    placement = assign_partitions([ranks] * len(todo), "batch")
    srcs = data["sources"] if level != "L8" else [None] * len(todo)
    for i, ((rpm, th), (part, acct)) in enumerate(zip(todo, placement)):
        rid = run_id(int(level[1:]), rpm, th, tag="c2d")
        argv, note = sbatch_argv(level, rid, rpm, th, partition=part, account=acct)
        print(" ".join(argv), "#", note, f"binary={BINARY_L8}")
        if args.go:
            check_mpi_partition(part)
            if level == "L8":
                _real_submit_l8(rid, rpm, th, part, acct)
            elif level == "L9":
                _real_submit_warm(rid, rpm, th, srcs[i], part, acct)
            else:
                raise SystemExit("L10 warm submission is not wired yet")
    for p in PARTITIONS:
        n = sum(1 for q, _ in placement if q == p)
        print(f"# partition {p}: {n} jobs x {ranks} ranks = {n * ranks} cpu concurrently "
              f"(cap {PARTITIONS[p][1]} cpu, plan 2)")
    if not args.go:
        print(f"(dry run: {len(todo)} commands printed; nothing submitted)")


N_WARM_L9, N_POST_L9 = 12, 40          # plan 1.3: warm-up cycles, then release, then post cycles


def check_mpi_partition(partition: str) -> None:
    """MPI jobs run on batch only. batch2 gave a UCX wireup failure on gpu4001 (job 8340753, 2026-11-05)."""
    if partition != "batch":
        raise ValueError(f"MPI jobs run on batch only (got {partition}: UCX wireup failure on gpu4001)")


def l9_warm_params(src_params: dict, source_run: str, t_ck: float, rpm: float, theta: float, binary: str) -> dict:
    """L9 warm start (plan 1.3, W9 type): the completed L8 source at the same (rpm, theta), fidelity 9, the same binary,
    t_checkpoint = the source's last logged time, N_WARM_L9 warm-up cycles, release, N_POST_L9 post cycles."""
    from scripts.postprocess import _t_scales
    p = {k: v for k, v in src_params.items() if not k.endswith("_prev")}   # the source's own restart keys
    for k in ("_canonical_run_dir", "next_run_id", "_err_file", "_source_run"):
        p.pop(k, None)
    _, T = _t_scales(p)
    p.update(run_id=run_id(9, rpm, theta, tag="c2d"), fidelity=9, t_checkpoint=float(t_ck),
             n_mix_cycles=N_WARM_L9, t_end=round((N_WARM_L9 + N_POST_L9) * T, 6),
             _binary=binary, _source_run=source_run)   # forcing kept exactly as the source (no rpm change)
    # The source forcing (pitfall "restart-ramp", diary 2026-11-13): without it the solver ramps the rocking
    # again from ZERO after the restart; the first 8 c2d_l9 runs had none and were still 8-19 % above the L8
    # velocity at the release. A binary that reads the dump's label takes it from checkpoint.phase and refuses
    # *_prev next to it; an older binary needs *_prev = the source's forcing.
    from scripts.simulate import binary_reads_label
    if not binary_reads_label(binary):
        z3 = [0.0, 0.0, 0.0]
        p.update(omega_b_prev=float(src_params["omega_b"]), theta_max_prev=list(src_params["theta_max"]),
                 phi_angular_prev=list(src_params.get("phi_angular", z3)),
                 amplitude_h_prev=list(src_params.get("amplitude_h", z3)),
                 phi_horizontal_prev=list(src_params.get("phi_horizontal", z3)),
                 omega_h_prev=float(src_params.get("omega_h", 0.0)))
    if abs(float(src_params["omega_b"]) * 60.0 / (2 * math.pi) - rpm) > 0.01 or \
            abs(float(src_params["theta_max"][0]) - theta) > 0.01:
        raise ValueError("source forcing does not match (rpm, theta) of the point")
    return p


def _real_submit_warm(rid: str, rpm: float, th: float, source: str, partition: str, account: str) -> None:
    """Real L9 warm submission. Refuses a non-batch partition and a cross-binary source. Stamps _err_file after sbatch
    (the job's log name is only known then), since the settled gate reads it."""
    from scripts import simulate
    check_mpi_partition(partition)
    require_warm_sources([[rpm, th]], {(rpm, th): source}, BINARY_L9, runs_root=ROOT / "runs")
    sd = ROOT / "runs" / source
    p = l9_warm_params(json.loads((sd / "params.json").read_text()), source, source_checkpoint_time(source),
                       rpm, th, BINARY_L9)
    if p["run_id"] != rid:
        raise SystemExit(f"run id mismatch {p['run_id']} != {rid}")
    jid = simulate.submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs",
                                walltime=STAGES["L9"]["walltime"], ntasks=STAGES["L9"]["ranks"],
                                template=ROOT / "config" / "slurm_mpi_template.sh",
                                checkpoint=str((sd / "checkpoint.dump").resolve()),
                                exclude=EXCLUDE, partition=partition, account=account, cpus=1, mem="4G")
    rp = ROOT / "runs" / rid / "params.json"
    stamped = warm_provenance(json.loads(rp.read_text()), source, ROOT / "logs" / f"slurm_{jid}.err")
    rp.write_text(json.dumps(stamped, indent=1))
    STATE.mkdir(parents=True, exist_ok=True)
    with JOBS.open("a") as f:
        f.write(json.dumps(dict(job=str(jid), level="L9", run=rid, rpm=rpm, theta=th, ranks=STAGES["L9"]["ranks"],
                                partition=partition, account=account, source=source,
                                t_checkpoint=p["t_checkpoint"], submitted=datetime.now(timezone.utc).isoformat())) + "\n")
    print(f"    job={jid}")


def _real_submit_l8(rid: str, rpm: float, th: float, partition: str, account: str) -> None:
    """Real submission path. Reuses simulate.submit_slurm (partition and account passed through)."""
    from scripts import simulate
    p = json.loads((ROOT / "runs/kmix_l8_rpm32.5/params.json").read_text())
    p.update(run_id=rid, omega_b=2 * math.pi * rpm / 60.0, theta_max=[th, 0.0, 0.0], _binary=BINARY_L8)
    jid = simulate.submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs",
                                walltime=STAGES["L8"]["walltime"], ntasks=STAGES["L8"]["ranks"],
                                exclude=EXCLUDE, partition=partition, account=account,
                                cpus=1, mem="4G")  # MPI: 1 cpu per rank, 4G per rank (template default); the
                                # submit_slurm defaults (cpus=4, 12G) gave 64-cpu jobs on 2026-11-04 (cancelled)
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
GATE_BLOCK_CYCLES = 3                  # drift/height use 3-cycle block means (2026-11-07: L9 3T limit cycles)
GATE_HEIGHT_TOL = 0.01                 # free-surface height range at release vs the source run
GATE_FS_LINE = "fs sanity [post-solid]: OK"   # must be in the run's .err (BioReactor3D CLAUDE.md, restart rule)
GATES_JSON = STATE / "gates.json"


def _rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b else math.inf


def settled_gate(t, tau95, edr, period, warm_end, release_t, height_range, src_height_range, err_text):
    """Settled-state gate for one warm-started run. Returns (status, detail).

    Drift: means of tau_95 and EDR over the last 3 blocks of GATE_BLOCK_CYCLES complete cycles after warm_end. The two block-to-block
    relative changes must each be below GATE_DRIFT_TOL. Height: free-surface range at release vs the source run,
    relative change below GATE_HEIGHT_TOL (src_height_range = the run's own last cycle since 2026-11-07). fs sanity: GATE_FS_LINE must be in err_text.
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
    B = GATE_BLOCK_CYCLES
    if n_done >= 3 * B:
        ks = range(n_done - 3 * B, n_done, B)
        means = {"tau95": [], "edr": []}
        for k in ks:
            m = post & (t >= warm_end + k * period) & (t < warm_end + (k + B) * period)
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


def last_cycle_height(v, period: float, n: int = 1) -> float:
    """Mean free-surface range (posY_max - posY_min) over the last n complete cycles of a vol_frac_interf array."""
    tend = float(v[:, 1].max())
    m = v[:, 1] >= tend - n * period
    return float((v[m, 4] - v[m, 5]).mean())


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
        m = (v[:, 1] >= release_t) & (v[:, 1] < release_t + GATE_BLOCK_CYCLES * period)
        if m.any():
            height_range = float(h[m].mean())
        # Reference = the run's OWN wave height in its last complete cycle (2026-11-07). The old reference was the
        # source run's wave; a cross-level source (L8 -> L9) gives a different wave near the 35.2 rpm mode by
        # grid effect alone, so steady runs failed. The check now asks: is the wave the same at release and at end?
        if float(v[:, 1].max()) >= release_t + 2 * GATE_BLOCK_CYCLES * period:
            src_height = last_cycle_height(v, period, GATE_BLOCK_CYCLES)
    err_text = None
    if params.get("_err_file") and Path(params["_err_file"]).exists():
        err_text = Path(params["_err_file"]).read_text(errors="replace")
    # `release_t` is passed only when the run has reached it, so the release check can be PENDING.
    reached = bool(t.max() >= release_t)
    return settled_gate(t, tau95, edr, period, warm_end, release_t if reached else None,
                        height_range, src_height, err_text)


def gate_all_runs(verbose: bool = False) -> dict[str, str]:
    """Gate every warm-started run in runs/ (t_checkpoint > 0) and write state/gates.json. Returns name -> status."""
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
        if verbose:
            print(f"GATE {d.name:24s} {status}")
    STATE.mkdir(parents=True, exist_ok=True)
    GATES_JSON.write_text(json.dumps(gates, indent=1))
    return gates


def exclude_failed_gate(d, gates: dict):
    """Drop the rows of FAILED-GATE runs (warm starts that did not settle). Rows without a gate entry (cold runs,
    or PASS / PENDING) are kept. `d` needs a 'run' column (the run folder name)."""
    bad = {k for k, v in gates.items() if v == "FAILED-GATE"}
    return d[~d["run"].isin(bad)].reset_index(drop=True)


def warm_provenance(params: dict, source_run, err_file) -> dict:
    """Params of a warm submission with the keys the settled gate reads: _source_run (the source run folder, for the
    height check) and _err_file (the .err that must hold the fs sanity line). Refuses anything else, so a warm run
    cannot be submitted un-gateable (plan 1.4 / 3)."""
    if float(params.get("t_checkpoint", 0.0)) <= 0.0:
        raise ValueError("not a warm start: t_checkpoint must be > 0 for a warm submission")
    if not source_run:
        raise ValueError("warm submission needs _source_run (the gate reads the source run's last cycle)")
    if not err_file:
        raise ValueError("warm submission needs _err_file (the gate reads the fs sanity line from it)")
    return {**params, "_source_run": str(source_run), "_err_file": str(err_file)}


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
    gates = gate_all_runs(verbose=True)
    print(f"warm-started runs gated: {len(gates)}; FAILED-GATE (excluded from fit): "
          f"{sum(v == 'FAILED-GATE' for v in gates.values())}")


def cmd_check_grid(_args) -> None:
    """Hard guard (lead 2026-11-05): every L9/L10 job in jobs.jsonl whose first logstats line is not the uniform level
    grid (total cells = 4**level, per-rank #Cells x ranks) is cancelled, and the cancellation is printed."""
    jobs = [json.loads(x) for x in JOBS.read_text().splitlines()] if JOBS.exists() else []
    scratch = Path(f"{paths.SCRATCH}/mpi_runs")
    for j in jobs:
        if j["level"] not in ("L9", "L10"):
            continue
        ok, total = check_warm_grid(scratch / j["run"] / "logstats.dat", int(j["level"][1:]), j["ranks"])
        if ok is None:
            print(f"{j['job']} {j['run']}: no logstats yet")
        elif ok:
            print(f"{j['job']} {j['run']}: grid OK, total {total}")
        else:
            print(f"{j['job']} {j['run']}: GRID WRONG, total {total} (want {4 ** int(j['level'][1:])}): cancelling")
            subprocess.run(["scancel", j["job"]], check=False)


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


# ----------------------------------------------------------------- fit (baseline, 2026-11-04)
# Model per QoI, on ln(QoI): LF = LooKRR (CachedLooKRR, reused) trained on L8 in the domain;
# HF = FullBayesKRRLRGP (reused, mfbml_local/krr_lr_gpr_fullbayes.py) trained on L9 and L10 in the domain.
# Criterion: c(x) = 2 sigma/mean on the linear scale = 2 sqrt(exp(v) - 1), v = sum_g w_g sc_g^2 nu/(nu-2) + Var_g(mu_g).
# Map: dense 101 x 101 grid over the domain, then the top 10 grid points refined by L-BFGS-B (bounds = domain).
DOMAIN = np.array([list(RPM), list(THETA)])
DB_KLA = ROOT / "experiments/multifidelity/kla_mix_dataset.csv"
DB_HYDRO = ROOT / "experiments/multifidelity/hydro_dataset_v2.json"
DOCS_FIT = ROOT / "docs" / "campaign2d"
# edr_95 is NOT saved for the existing runs (plan 1): the baseline uses edr_mean_t, labelled as a stand-in.
QOI_SOURCE = {
    "tau_95": ("hydro", "tau_95_t"),
    "edr_95": ("hydro", "edr_mean_t"),
    "kLa_1T_10": ("kla", "kLa_1T_10"),
    "kLa_1T_25": ("kla", "kLa_1T_25"),
    "kLa_1T_50": ("kla", "kLa_1T_50"),
    "dtmix_0.50": ("kla", "dtmix_0.50"),
    "dtmix_0.25": ("mix", "dtmix_0.25"),   # 2026-11-09: BO objective; computed from tr_oxy.dat (mix_value), cached
}
MIX_CACHE = ROOT / "experiments/multifidelity/dtmix_cache.json"


def mix_value(run_dir, col):
    """dtmix at one threshold from the run's tracer file via postprocess._compute_mixing_metrics (bulk columns when
    present, else the old method), cached by run name and tr_oxy.dat mtime. None if no tracer file."""
    run_dir = Path(run_dir)
    tr = run_dir / "tr_oxy.dat"
    if not (tr.exists() and (run_dir / "params.json").exists()):
        return None
    try:   # concurrent readers/writers (12 benchmark jobs + autopilot, 2026-11-09): a torn read means 'no cache'
        cache = json.loads(MIX_CACHE.read_text()) if MIX_CACHE.exists() else {}
    except (json.JSONDecodeError, OSError):
        cache = {}
    key, mt = run_dir.name, tr.stat().st_mtime
    if key not in cache or cache[key].get("mtime") != mt:
        from scripts.postprocess import _compute_mixing_metrics
        params = json.loads((run_dir / "params.json").read_text())
        try:
            r = _compute_mixing_metrics(run_dir, params)
        except Exception:  # noqa: BLE001  (a broken old run gives no value, never a crash of the loop)
            r = {}
        cache[key] = {"mtime": mt, **{k: (None if v is None or not math.isfinite(v) else float(v))
                                      for k, v in r.items() if k.startswith("dtmix")}}
        tmp = MIX_CACHE.with_name(f"{MIX_CACHE.name}.{os.getpid()}.tmp")
        tmp.write_text(json.dumps(cache, indent=0)); os.replace(tmp, MIX_CACHE)   # atomic on one filesystem
    v = cache[key].get(col)
    return None if v is None else float(v)
CRIT_QOIS = ["tau_95", "edr_95", "kLa_1T_10"]        # gated (<= 0.30); the others are reported, not gated
MAP_N = 101


def dedup_runs(d):
    """One row per run folder (first kept). experiments/multifidelity/kla_mix_dataset.csv lists six runs in two sets
    ('design' and 't5_grid'); without this the fit counts them twice (found 2026-11-04)."""
    return d.drop_duplicates(subset="run", keep="first").reset_index(drop=True)


def campaign_row(run_dir: Path, col: str, src: str):
    """One finished campaign run as a fit row (level, rpm, theta, run, y), or None if it has not finished.
    kla: the value from results.json. hydro: hydro_dataset_v2.qois (the database definition, so the stand-in
    edr_mean_t stays the same QoI for every row)."""
    run_dir = Path(run_dir)
    pf, rf = run_dir / "params.json", run_dir / "results.json"
    if not (pf.exists() and rf.exists()):
        return None
    p = json.loads(pf.read_text())
    if src == "kla":
        y = json.loads(rf.read_text()).get(col)
        if y is None:
            return None
    elif src == "mix":
        y = mix_value(run_dir, col)
        if y is None:
            return None
    else:
        if not (run_dir / "shear_stress.dat").exists():
            return None
        from scripts.hydro_dataset_v2 import qois
        y = qois(run_dir.name)[col]
    rpm = round(float(p["omega_b"]) * 60.0 / (2 * math.pi), 3)
    return dict(level=int(p["fidelity"]), rpm=rpm, theta=float(p["theta_max"][0]), run=run_dir.name, y=float(y))


def edr95_late(run_dir):
    """edr_95 for the fit (2026-11-05): median over the last 3 whole cycles of edr_pct.dat (raw column 2), scaled to
    W/m^3 as postprocess.py (rho U^3 / L). Used for c2d runs and for the e95 backfill restarts; both end at the last
    checkpoint (backfill_window_test.py: late vs pre-release EDR agrees within 0.63 % on 4 database runs). None if the
    file is absent (the run predates the 2026-11-04 binary)."""
    import pandas  # noqa: F401  (keeps the import order of load_db)
    from scripts.postprocess import _load_dat, _t_scales
    run_dir = Path(run_dir)
    f, pf = run_dir / "edr_pct.dat", run_dir / "params.json"
    if not (f.exists() and pf.exists()):
        return None
    p = json.loads(pf.read_text())
    T_bio, T = _t_scales(p)
    a = p.get("geometry", {}).get("a", 0.25)
    escale = 1000.0 * (a / T_bio) ** 3 / a
    arr = _load_dat(f)
    if arr is None or arr.ndim != 2 or arr.shape[0] < 3:
        return None
    t, v = arr[:, 1], arr[:, 2]
    m = t >= t.max() - 3.0 * T
    return float(np.median(v[m]) * escale)


def load_edr95_db():
    """edr_95 rows (2026-11-05): the tau_95 row set (database and c2d) with the value edr95_late of the run folder:
    the run itself for c2d runs, the _e95 backfill restart for database runs. No stand-in. Rows without a value drop."""
    import pandas as pd
    base = load_db("tau_95")[["level", "rpm", "theta", "run"]].drop_duplicates("run")
    out = []
    for r in base.itertuples(index=False):
        folder = ROOT / "runs" / (r.run if r.run.startswith("c2d_") else r.run + "_e95")
        y = edr95_late(folder)
        if y is not None and y > 0:
            out.append(dict(level=int(r.level), rpm=float(r.rpm), theta=float(r.theta), run=r.run, y=y))
    return pd.DataFrame(out, columns=["level", "rpm", "theta", "run", "y"])


def load_db(qoi: str):
    """(level, rpm, theta, value) rows of the existing database for one QoI, domain-restricted, levels 8-10."""
    if qoi == "edr_95":
        return load_edr95_db()
    import pandas as pd
    kla = pd.read_csv(DB_KLA)
    src, col = QOI_SOURCE[qoi]
    if src == "kla":
        d = kla[["level", "rpm", "theta", "run", col]].rename(columns={col: "y"})
    elif src == "mix":
        d = kla[["level", "rpm", "theta", "run"]].copy()
        d["y"] = [mix_value(ROOT / "runs" / r, col) for r in d.run]
        d["y"] = d.y.astype(float)
    else:
        h = pd.DataFrame(json.loads(DB_HYDRO.read_text()))
        d = h[["run", "level", "rpm", col]].merge(kla[["run", "theta"]], on="run")
        d = d[["level", "rpm", "theta", "run", col]].rename(columns={col: "y"})
    # Finished campaign runs (runs/c2d_*) are not in the database files: add them here (2026-11-04).
    crow = [r for r in (campaign_row(p, col, src) for p in sorted((ROOT / "runs").glob("c2d_*"))) if r is not None]
    if crow:
        d = pd.concat([d, pd.DataFrame(crow)], ignore_index=True)
    d = dedup_runs(d.dropna())
    d = d[d.level.isin([8, 9, 10]) & (d.rpm >= RPM[0]) & (d.rpm <= RPM[1])
          & (d.theta >= THETA[0]) & (d.theta <= THETA[1]) & (d.y > 0)]
    return d.reset_index(drop=True)


# ----------------------------------------------------------------- LF/HF fit with a noise floor (2026-11-04)
# Lead decision: (a) restart selection by marginal likelihood among restarts with noise_std >= the run-to-run floor;
# (b) the optimiser's noise lower bound is that floor. Old baseline (seed 42, one start, unfloored, LOO) is superseded:
# the LOO objective has several local optima and the fit's seed chose between them (diary 2026-11-04).
NOISE_FLOOR_LN = {"kLa_1T_10": 0.095, "kLa_1T_25": 0.189, "dtmix_0.50": 0.104, "tau_95": 0.017,   # diary 2026-10-15 E0
                  "dtmix_0.25": 0.104}   # 2026-11-09: not measured for chi 0.25; the dtmix_0.50 value is used
LF_SEEDS = tuple(range(1, 11))


def noise_floor(qoi: str) -> float:
    """Run-to-run ln sd floor. Without a measured floor: 0.5 x the smallest measured floor (lead, 2026-11-04)."""
    return NOISE_FLOOR_LN.get(qoi, 0.5 * min(NOISE_FLOOR_LN.values()))


_LF_CLS = None


def _lf_class():
    """FastLooKRR (closed-form LOO, production LF) with the noise lower bound at self.noise_floor (ln units)."""
    global _LF_CLS
    if _LF_CLS is None:
        from scripts.fast_loo_krr import FastLooKRR

        class FlooredLOO(FastLooKRR):
            noise_floor = 0.0

            def _bound_definition_for_optimization(self, noise_bound=(10 ** -3, 10)):
                b = super()._bound_definition_for_optimization(list(noise_bound))
                if self.noise_data and self.params_optimize:
                    b[-1, 0] = max(b[-1, 0], self.noise_floor / self.y_std)
                return b

        _LF_CLS = FlooredLOO
    return _LF_CLS


def make_lf(seed: int, floor: float):
    lf = _lf_class()(design_space=DOMAIN, params_optimize=True, noise_data=True, optimizer_restart=0, seed=seed)
    lf.noise_floor = float(floor)
    return lf


def lf_log_ml(lf) -> float:
    """Log marginal likelihood of the LF kernel model on its training data (scaled ln outputs)."""
    X, y = lf.sample_x_scaled, np.asarray(lf.sample_y_scaled, float).ravel()
    K = lf._training_kernel_matrix(scaled_x=X, scaled_noise_std=lf.noise_std / lf.y_std)
    L = np.linalg.cholesky(K)
    a = np.linalg.solve(L, y)
    return float(-0.5 * a @ a - np.sum(np.log(np.diag(L))) - 0.5 * len(y) * np.log(2 * np.pi))


def select_restart(rows) -> int:
    """Index of the restart with the best marginal likelihood among rows with ok=True (noise_std >= floor)."""
    ok = [i for i, r in enumerate(rows) if r["ok"]]
    if not ok:
        raise ValueError("no LF restart with noise_std at or above the floor")
    return max(ok, key=lambda i: rows[i]["lml"])


def hf_grid(floor: float, ystd: float):
    """FullBayes (length, noise ratio, log prior) grid with the noise ratio floor/ystd as its lower bound.
    If that bound exceeds the prior's upper limit 0.3, the upper limit is raised to 4 x the bound (reported)."""
    from mfbml_local.krr_lr_gpr_fullbayes import default_grid
    eta_lo = max(1e-3, float(floor) / float(ystd))
    return default_grid(eta_lo=eta_lo, eta_hi=max(0.3, 4.0 * eta_lo))


def fit_lf_selected(X, Y, floor: float):
    """Run the LF optimisation from LF_SEEDS (one start each), keep the best-LML restart with noise >= floor."""
    models, rows = [], []
    for s in LF_SEEDS:
        np.random.seed(s)
        lf = make_lf(seed=s, floor=floor)
        lf.train(X, Y, portion_test=0.2)
        noise = float(np.asarray(lf.noise_std).ravel()[0])
        rows.append(dict(seed=s, noise_std=noise, lml=lf_log_ml(lf), ok=noise >= floor - 1e-12,
                         lengthscales=np.round(np.asarray(lf.kernel.param, float), 4).tolist()))
        models.append(lf)
    idx = select_restart(rows)
    return models[idx], rows, models


def fit_qoi_info(d, qoi: str):
    """LF (floored, seed selected by LML) on L8, then FullBayesKRRLRGP on L9+L10 with the floored grid.
    Returns (HF model, info) with the per-seed LF table and the HF evidence per seed (same rule, for the check)."""
    from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP
    lo, hi = d[d.level == 8], d[d.level >= 9]
    if len(lo) < 5 or len(hi) < 2:
        raise SystemExit(f"too few points: L8 {len(lo)}, L9+L10 {len(hi)}")
    floor = noise_floor(qoi)
    Xl, Yl = lo[["rpm", "theta"]].to_numpy(float), np.log(lo.y.to_numpy(float)).reshape(-1, 1)
    Xh, Yh = hi[["rpm", "theta"]].to_numpy(float), np.log(hi.y.to_numpy(float))
    lf, rows, models = fit_lf_selected(Xl, Yl, floor)
    grid = hf_grid(floor, float(np.std(Yh)))
    m = FullBayesKRRLRGP(DOMAIN, lf_model=lf, grid=grid)
    m.train(Xh, Yh)
    per_seed = []
    for r, lf_s in zip(rows, models):
        ms = FullBayesKRRLRGP(DOMAIN, lf_model=lf_s, grid=grid)
        ms.train(Xh, Yh)
        lw = ms.logprior + ms.logml
        per_seed.append(dict(seed=r["seed"], lf_lml=r["lml"], lf_noise=r["noise_std"], hf_log_evidence=float(
            np.log(np.sum(np.exp(lw - lw.max()))) + lw.max())))
    sel = [r["seed"] for r, _ in zip(rows, models)][select_restart(rows)]
    return m, dict(floor=floor, selected_seed=sel, lf_rows=rows, hf_per_seed=per_seed,
                   hf_eta_floor=floor / float(np.std(Yh)), hf_eta_hi=float(grid[1].max()))


def fit_qoi(d, qoi: str):
    return fit_qoi_info(d, qoi)[0]


def crit(m, X: np.ndarray, chunk: int = 400) -> np.ndarray:
    """2 sigma/mean (linear scale) of the model at rows X (k, 2) = (rpm, theta). Chunked: predict_nodes is (G, n, m)."""
    out = np.empty(len(X))
    for s in range(0, len(X), chunk):
        mu, sc, nu, w = m.predict_dist(X[s:s + chunk])
        wm = w @ mu
        v = (w[:, None] * sc ** 2 * (nu / (nu - 2.0))).sum(0) + (w[:, None] * (mu - wm) ** 2).sum(0)
        out[s:s + chunk] = 2.0 * np.sqrt(np.expm1(v))
    return out


def max_crit(m) -> dict:
    """Max of crit over the domain: dense MAP_N x MAP_N grid, then L-BFGS-B from the top 10 grid points."""
    from scipy.optimize import minimize
    r = np.linspace(*RPM, MAP_N)
    t = np.linspace(*THETA, MAP_N)
    RR, TT = np.meshgrid(r, t, indexing="ij")
    G = np.column_stack([RR.ravel(), TT.ravel()])
    cg = crit(m, G)
    top = np.argsort(cg)[::-1][:10]
    best_x, best_v = G[top[0]].copy(), float(cg[top[0]])
    for i in top:
        res = minimize(lambda x: -float(crit(m, x[None])[0]), G[i], method="L-BFGS-B",
                       bounds=[tuple(RPM), tuple(THETA)])
        if -res.fun > best_v:
            best_v, best_x = float(-res.fun), np.asarray(res.x, float)
    return dict(max=best_v, rpm=float(best_x[0]), theta=float(best_x[1]), grid_max=float(cg.max()),
                grid_argmax=[float(G[np.argmax(cg)][0]), float(G[np.argmax(cg)][1])],
                grid=cg.reshape(MAP_N, MAP_N), r=r, t=t)


def plot_map(qoi: str, res: dict, d, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    im = ax.pcolormesh(res["r"], res["t"], res["grid"].T, shading="auto", cmap="viridis", vmin=0, vmax=max(0.5, res["max"]))
    for lev, mk, sz in ((8, ".", 10), (9, "o", 18), (10, "s", 18)):
        s = d[d.level == lev]
        ax.scatter(s.rpm, s.theta, marker=mk, s=sz, facecolors="none" if lev > 8 else "k", edgecolors="w" if lev == 8 else "w", linewidths=0.6)
    ax.set_xlabel("rpm")
    ax.set_ylabel("theta (deg)")
    ax.set_title(qoi, fontsize=9)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("2 sigma / mean", fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fit(_args=None) -> None:
    """Baseline fit for every QoI on the existing database (domain subset). Refuses to run outside a compute job."""
    if "SLURM_JOB_ID" not in os.environ:
        raise SystemExit("fit is heavy: run it inside a compute job (srun --overlap or sbatch), never on the login node")
    out = {}
    gates = gate_all_runs()
    for qoi in QOI_SOURCE:
        d_all = load_db(qoi)
        d = exclude_failed_gate(d_all, gates)
        n_excl = len(d_all) - len(d)
        m, finfo = fit_qoi_info(d, qoi)
        res = max_crit(m)
        png = DOCS_FIT / f"fit0_{qoi.replace('.', '_')}.png"
        plot_map(qoi, res, d, png)
        out[qoi] = dict(max=res["max"], rpm=res["rpm"], theta=res["theta"], grid_max=res["grid_max"],
                        grid_argmax=res["grid_argmax"], n_L8=int((d.level == 8).sum()),
                        n_HF=int((d.level >= 9).sum()), n_failed_gate_excluded=n_excl, floor=finfo["floor"], selected_seed=finfo["selected_seed"], lf_seeds=finfo["lf_rows"], hf_per_seed=finfo["hf_per_seed"],
                        gated=qoi in CRIT_QOIS, png=str(png))
        print(f"{qoi:12s} max 2s/mean = {res['max']:.3f} at rpm {res['rpm']:.2f} theta {res['theta']:.2f} "
              f"(grid max {res['grid_max']:.3f}) L8 n={out[qoi]['n_L8']} HF n={out[qoi]['n_HF']} "
              f"{'GATED' if qoi in CRIT_QOIS else 'reported'}")
    STATE.mkdir(parents=True, exist_ok=True)
    (STATE / "fit0.json").write_text(json.dumps(out, indent=1))


def cmd_fit(args) -> None:
    fit(args)


# ----------------------------------------------------------------- next (acquisition, 2026-11-04)
NEXT_SEP = 0.05          # min distance between picks, as a fraction of the domain size in each axis (Euclidean)
# Measured warm medians (sacct cost audit, cost_table.tsv, warm rows, core-h per run, corehc_all): L9 n=54 97.7,
# L10 n=82 764.5 (2026-11-04). Lead rule: every L9 and L10 point warm-starts from a COMPLETED L8 run at the same
# (rpm, theta). Cold L9/L10 points are refused.
WARM_MEDIAN_CORE_H = {9: 97.7, 10: 764.5}
NEXT_COST_L9 = WARM_MEDIAN_CORE_H[9]


def l8_source_map(runs_root=None) -> dict:
    """(rpm, theta) -> run folder of a COMPLETED L8 run: params fidelity 8, results.json and checkpoint.dump present."""
    root = Path(runs_root) if runs_root else ROOT / "runs"
    out = {}
    for d in sorted(root.glob("*")):
        pf = d / "params.json"
        if not (pf.exists() and (d / "results.json").exists() and (d / "checkpoint.dump").exists()):
            continue
        try:
            p = json.loads(pf.read_text())
        except json.JSONDecodeError:
            continue
        if p.get("fidelity") != 8:
            continue
        key = (round(float(p["omega_b"]) * 60.0 / (2 * math.pi), 3), float(p["theta_max"][0]))
        out.setdefault(key, d.name)
    return out


def binary_sha256(path) -> str | None:
    """sha256 of a binary file, or None if the file is not on disk (the hash is then unknown)."""
    import hashlib
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return None
    return hashlib.sha256(p.read_bytes()).hexdigest()


ALLOWED_EXTRA_FLAGS = frozenset({"-DCROSS_LEVEL_WARMSTART=1"})   # the only flag a same-source binary may differ by


def build_info(binary):
    """Build sidecar <binary>.build.json: {commit, flags}. None if absent (the source is then unknown)."""
    if not binary:
        return None
    p = Path(str(binary) + ".build.json")
    return json.loads(p.read_text()) if p.exists() else None


def same_source(a, b) -> bool:
    """True if both binaries come from the same source commit and differ only by ALLOWED_EXTRA_FLAGS (lead, 2026-11-05).
    A binary without a sidecar has an unknown source and is not the same source."""
    ia, ib = build_info(a), build_info(b)
    if ia is None or ib is None or ia.get("commit") != ib.get("commit"):
        return False
    fa = set(ia.get("flags", [])) - ALLOWED_EXTRA_FLAGS
    fb = set(ib.get("flags", [])) - ALLOWED_EXTRA_FLAGS
    return fa == fb and not ((set(ia.get("flags", [])) | set(ib.get("flags", []))) - ALLOWED_EXTRA_FLAGS - fa)


def same_source_map(target_binary, runs_root=None) -> dict:
    """l8_source_map restricted to sources whose binary is the same source as the target (see same_source)."""
    root = Path(runs_root) if runs_root else ROOT / "runs"
    out = {}
    for key, name in l8_source_map(root).items():
        sb = json.loads((root / name / "params.json").read_text()).get("_binary")
        if same_source(sb, target_binary):
            out[key] = name
    return out


same_binary_map = same_source_map   # old name, kept for callers


def require_warm_sources(points, smap: dict, target_binary, runs_root=None) -> list:
    """Source run folder for each (rpm, theta) point. Raises if a point has no completed L8 source (cold refused), or if
    a source is not the same source as the target binary (cross-source refused, lead 2026-11-05)."""
    root = Path(runs_root) if runs_root else ROOT / "runs"
    if not Path(target_binary).exists():
        raise ValueError(f"target binary {target_binary} is not on disk")
    missing = [p for p in points if (float(p[0]), float(p[1])) not in smap]
    if missing:
        raise ValueError(f"no completed warm source L8 run at {missing}: cold L9/L10 points are refused")
    srcs = [smap[(float(p[0]), float(p[1]))] for p in points]
    for name in srcs:
        sb = json.loads((root / name / "params.json").read_text()).get("_binary")
        if not same_source(sb, target_binary):
            raise ValueError(f"source {name} binary {sb} is not the same source as the target (source refused)")
    return srcs


def grid_ok(total_cells: int, level: int) -> bool:
    """A uniform level-L grid in 2D has 4**L leaves. Anything else is not the level-L grid (lead 2026-11-05)."""
    return int(total_cells) == 4 ** int(level)


def first_cells_per_rank(text: str):
    """Per-rank #Cells of the first logstats line, or None."""
    m = re.search(r"#Cells:\s*(\d+)", text)
    return int(m.group(1)) if m else None


def check_warm_grid(logstats_path, level: int, ranks: int) -> tuple:
    """(ok, total) from the first logstats line: total = per-rank #Cells x ranks (per-rank count, lead 2026-11-05)."""
    p = Path(logstats_path)
    if not p.exists():
        return None, None
    first = next((ln for ln in p.read_text(errors="replace").splitlines() if "#Cells" in ln), None)
    per = first_cells_per_rank(first or "")
    if per is None:
        return None, None
    total = per * int(ranks)
    return grid_ok(total, level), total


def source_checkpoint_time(run: str) -> float:
    """Last logged time of the source run (the warm start time). Same definition as submit_continue.py."""
    import numpy as np
    return float(np.loadtxt(ROOT / "runs" / run / "shear_stress.dat", skiprows=1)[-1, 1])


def select_next(sites, scores: dict, k: int, sep: float = NEXT_SEP) -> list:
    """Greedy top-k by score. A site within `sep` (domain-normalised) of an already picked site is skipped, so one
    ridge does not take the whole batch. Returns (rpm, theta, score) tuples, best first."""
    picked = []
    for s in sorted(sites, key=lambda s: -scores[s]):
        if len(picked) >= k:
            break
        if all(math.hypot((s[0] - p[0]) / (RPM[1] - RPM[0]), (s[1] - p[1]) / (THETA[1] - THETA[0])) >= sep
               for p in picked):
            picked.append((float(s[0]), float(s[1]), float(scores[s])))
    return picked


def cmd_next(args) -> None:
    """Acquisition for L9. Candidates: finished L8 sites in the domain (a warm L9 start needs the L8 checkpoint at
    the same rpm and theta). Score(x) = max over the gated QoIs of 2 sigma/mean at x, divided by the L9 warm cost.
    Gated QoIs and FAILED-GATE exclusion as in fit. Writes state/next_L9.json unless --dry-run. Compute node only."""
    if "SLURM_JOB_ID" not in os.environ:
        raise SystemExit("next is heavy: run it inside a compute job (srun --overlap or sbatch), never on the login node")
    gates = gate_all_runs()
    same = same_source_map(BINARY_L9)   # warm sources must be the same binary (lead 2026-11-04)
    sites = [(r, t) for r, t in sorted(same) if RPM[0] <= r <= RPM[1] and THETA[0] <= t <= THETA[1]]
    if not sites:
        raise SystemExit("no completed L8 source with the target binary in the domain: no L9 candidates")
    X = np.array(sites, float)
    worst = np.zeros(len(X))
    per_qoi = {}
    for qoi in CRIT_QOIS:
        d_all = load_db(qoi)
        d = exclude_failed_gate(d_all, gates)
        m, _finfo = fit_qoi_info(d, qoi)
        cq = crit(m, X)
        per_qoi[qoi] = cq
        worst = np.maximum(worst, cq)
        print(f"{qoi:12s} candidate max 2s/mean {cq.max():.3f} (n HF {int((d.level >= 9).sum())}, "
              f"FAILED-GATE rows dropped {len(d_all) - len(d)})")
    scores = {s: float(worst[i]) / NEXT_COST_L9 for i, s in enumerate(sites)}
    picks = select_next(sites, scores, k=args.k)
    print(f"candidates: {len(sites)} L8 sites; continuous-domain max is in fit0.json (not re-evaluated here)")
    for rpm, th, sc in picks:
        i = sites.index((rpm, th))
        print(f"  pick rpm {rpm:6.2f} theta {th:4.2f}  score {sc:.3e}  "
              + "  ".join(f"{q} {per_qoi[q][i]:.3f}" for q in CRIT_QOIS))
    if args.dry_run:
        print("(dry run: state/next_L9.json not written)")
        return
    smap = same_source_map(BINARY_L9)
    pts = [[rpm, th] for rpm, th, _ in picks]
    srcs = require_warm_sources(pts, smap, BINARY_L8)
    out = dict(level="L9", cost_per_point=NEXT_COST_L9, cost_basis="measured warm median, corehc_all, n=54",
               expected_core_h=round(NEXT_COST_L9 * len(pts), 1),
               source="warm from the completed L8 run at the same (rpm, theta)",
               points=pts, sources=srcs, checkpoint_t=[source_checkpoint_time(s) for s in srcs],
               scores=[sc for *_, sc in picks],
               qois=CRIT_QOIS, fit_gates=dict(n_gate_entries=len(gates),
                                              n_failed=sum(v == "FAILED-GATE" for v in gates.values())))
    (STATE / "next_L9.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {STATE / 'next_L9.json'}")


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
    sub.add_parser("fit").set_defaults(func=cmd_fit)
    nx = sub.add_parser("next")
    nx.add_argument("--k", type=int, default=4, help="L9 batch size (plan 1.3: 4)")
    nx.add_argument("--dry-run", action="store_true", help="print the picks, do not write state/next_L9.json")
    nx.set_defaults(func=cmd_next)
    sub.add_parser("ledger").set_defaults(func=cmd_ledger)
    sub.add_parser("check-grid").set_defaults(func=cmd_check_grid)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
