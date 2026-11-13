# Feature guardrails and the feature ledger — design

Status: draft for approval (2026-11-13). Owner: the code owner (Claude); the user takes high-level decisions only.

## 1. Why

Since August the solver and its workflow got many mechanisms. Some save compute (warm starts, record/replay,
the lean build, chaining). Others repair a defect that cost days (face fields after restore, the cs-masked
averages, the tracer carried across a restart). Each one was checked once, when it was new, and the checks
are spread over diary entries, one-off scripts and a few `medium` tests that each run their own simulation.
Two questions have no single answer today:

1. **Does each mechanism still work?** A refactor (like the 2026-11 src split) can break one without a
   failing test. The golden gate proves "nothing changed"; it does not say "this mechanism is correct".
2. **What is each one worth?** For a compute saver: cost saved against error added, at the level where we
   use it. Some answers exist (diary: exact replay is only 1.8x faster than coupled; L8 -> L10 warm start is
   within 5 % after ~12 cycles), but in prose, at different levels, never side by side.

Success: one command on a compute node produces a fresh bundle of runs; `pytest` then says, per mechanism,
"works" or "broken", in a few seconds and without computing; and one table + one figure show, per compute
saver, the cost it saves and the error it adds, at L6 (cheap) and at L8-L10 (where we use it).

## 2. Goals, non-goals, anti-goals

Goals
- A guard per mechanism (list in section 3), each with a stated pass bar and a demonstrated power
  (a planted defect makes it fail).
- Compute is shared: one bundle of short runs at one cheap level, many tests read it (the user's rule:
  "run one L_W and reuse it for many tests").
- A ledger of effects per compute saver: core-hours saved and error added per QoI (tau_95, EDR, kLa,
  mixing time), from existing production runs first, new runs only for gaps.

Non-goals
- New physics validation against Kim et al. (exists separately).
- Running the bundle in GitHub CI (no Basilisk/SLURM there). CI keeps its L3 `medium` tests.
- Merging before the DFD talk; changing solver physics.

Anti-goals (the work fails if any of these happen, even if tests pass)
- A test that runs a simulation when someone types `pytest` (2026-11-06: ~90 jobs from one pytest run).
- A guard that **skips silently** when its data is missing (seen today in test_figstyle). Missing or stale
  bundle data must be one visible failure, not N quiet skips.
- A tolerance chosen after seeing the result. Bars are written before the first bundle run, from the
  noise floor (two identical runs) where one exists.
- Bundle runs on a login node, or a bundle that needs more than one batch node for ~30 min.

## 3. Inventory: the mechanisms

Compute savers (guard + ledger row):

| # | Mechanism | Saves | Existing check |
|---|---|---|---|
| S1 | Exact flux record/replay (`RECORD_FLUX`/`REPLAY_FLUX`) | the flow solve for scalar runs | golden (bitwise vs itself), test_replay_gate (hpc), F3b |
| S2 | Looped-period replay (low fidelity, `replay_lf.py`) | long scalar runs | diary F2/periodic-flow results; test_replay_lf (params only) |
| S3 | Cross-level warm start L(k) -> L(k+1) and L8 -> L10 (`CROSS_LEVEL_WARMSTART`) | spin-up at the fine level | test_cross_level_warmstart (medium), x8/W9 diary tests |
| S4 | Cross-case warm start ("checkpointing", chain.py: condition k -> k+1, velocity rescale, phase fix) | spin-up per new condition | settling studies (diary), test_chain (unit) |
| S5 | Same-condition chaining: continuation restart, early dump before the wall limit, phase sidecar, autoextend | lost work at timeouts | test_restart_continue, test_graceful_timeout (medium), test_autoextend |
| S6 | Lean build (`EXTRA_TRACERS=0`) | ~40 % per step | golden (lean vs itself only) |
| S7 | MPI decomposition (rank count) | wall time | test_mpi_checkpoint_parity (medium) |
| S8 | Early tracer release after a warm start (12 instead of 80 cycles) | cycles before release | R12, W9 diary tests |

Correctness repairs (guard only, no ledger row): S9 cs/fs recomputed after restore (`assert_fs_sane`);
S10 tracer and oxygen carried across a restart; S11 cs-weighted, bag-masked averages; S12 stop on the
period boundary (`periods.py`); S13 campaign source rule and grid guard; S14 tracer seed in full liquid
cells (`SEED_KIM`); S15 cross-level pre-refine order (`CROSS_LEVEL_PREREFINE_SOLID`).

## 4. How (decided)

### 4.1 Produce once, assert many

`tools/bundle/` (grown from `tools/golden/`, which already builds 5 variants and runs 6 L6 cases in ~5 min):

- `bundle.sbatch SRC WORK` builds every variant once and runs a small DAG of cases. Each case reuses its
  parent's outputs (checkpoint, record file) instead of recomputing them.
- Output: `WORK/bundle.json` with, per case, the binary hash, the source commit, the outputs' paths, the
  QoIs, and the measured wall time (cost).
- Tests (new marker `bundle`, in the default run) read `bundle.json` through one fixture. The fixture
  fails once, loudly, when the bundle is missing or was built from another source commit than the
  working tree; it never computes.

Level: **L6** for every case (L5 only as the coarse source of S3). Reason: the guards test mechanisms
(does restart carry the tracer, does replay reproduce the scalars), which do not depend on the level; L6 is
the cheapest level that has a stable two-phase flow in the bag and is already proven in the golden gate.
L4 is too coarse for a meaningful mixing time; that is the measured limit, not a guess (to be confirmed in
step 1 of the plan).

Cases (draft; 1 node, serial except S7):

| Case | From | Guards |
|---|---|---|
| `cold` | — | reference: 3 periods, release after 1, record the flux for the last period, checkpoint at the end |
| `cold2` | — | same as `cold`, a second time: the noise floor (expected bitwise equal, OMP=1) |
| `cont` | `cold` at 2 periods | S5, S9, S10: restart + 1 period vs the last period of `cold` |
| `kill` | `cold`, wall limit forced | S5: early dump + phase sidecar, then resume, vs `cold` |
| `replay` | `cold`'s record | S1: scalars vs `cold` |
| `loop` | `cold`'s record, looped | S2: error vs `cold` |
| `lean` | — | S6: shared outputs vs `cold` |
| `l5` -> `xl6` | L5 cold -> L6 warm | S3, S15: vs `cold` after k cycles |
| `caseB_cold`, `caseB_warm` | `cold` at rpm A -> rpm B | S4, S8: warm vs cold at rpm B |
| `mpi4` | — (MPI node) | S7: vs `cold` |

S11-S14 keep their static/unit tests (they need no run); S11 also gets a check on `cold`'s outputs.

### 4.2 Bars and power

- Where two identical runs are bitwise equal (OMP=1), the bar is "bitwise" or a tolerance stated per
  output before the first run (restart and replay are not expected to be bitwise; bars from the noise
  floor of `cold` vs `cold2` and from the diary's production measurements).
- Every guard gets one planted defect (as for the golden gate): e.g. skip the fs recompute, drop the
  tracer from the dump list, shift the replay by one step, swap the restart phase. The guard must fail on
  it. Power results go in the diary.

### 4.3 The ledger

`experiments/feature_ledger.json` + `scripts/figures/fig_feature_ledger.py` (one figure: per mechanism,
core-hours saved vs error added, one point per QoI, L6 bundle value next to the production value).

- Production rows come from existing runs first (zero compute): warm vs cold pairs (xlevel vs l10c, W9 vs
  cold L9, x8 continuation), replay (F3b, periodic-flow results), lean (logged wall times), timeouts
  recovered (autoextend logs). A script reads the run folders; it does not copy numbers from the diary.
- Gaps (no production pair exists) are listed with a core-hour estimate and go to the user for a go/no-go
  before any run.
- The L6 value next to the production value also tests the idea behind the bundle: does the cheap level
  predict the effect at the fine level? If not for some mechanism, the ledger says so.

### 4.4 Retire duplicate compute

`medium` tests whose check moves into the bundle (e.g. test_mpi_checkpoint_parity, test_restart_continue)
stay in CI at L3 until the bundle tests have passed their power checks, then are reduced to what CI alone
can check.

## 5. Constraints

- Compute: bundle on one batch node, target <= 30 min (~10-15 core-h per run of the bundle); MPI case on
  a node allowed for MPI (not batch2). Exclude node2356, node2410. Never mbessa-condo.
- No running chain may be disturbed: work on a branch; merge after the campaigns end (as lean/phase2).
- Basilisk rules: no bare `return;` in events; recompute cs/fs after restore; qcc from
  /oscar/data/dharri15/eaguerov/basilisk.

## 6. Left to the builder

Exact case lengths (periods), output-by-output tolerances (written before the first run), the planted
defects, the ledger's JSON schema, and which existing production runs pair up for each ledger row.

## 7. My calls (check these)

- L6 bundle, L5 as the coarse source; L4 rejected (to confirm in step 1).
- The bundle replaces, not adds to, the per-test simulations of the `medium` tests (after power is shown).
- Ledger from existing production runs first; new production runs only with your go.
- S2 (looped replay) is in the ledger even though the diary already judged it "not a cost saver": the
  ledger should show dead ends too, so nobody re-tries them.
