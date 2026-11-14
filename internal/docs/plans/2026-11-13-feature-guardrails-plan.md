# Feature guardrails and feature ledger — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** One batch job builds every solver variant once and runs a small chain of L6 cases. Tests that only
read its results then say per mechanism "works" or "broken". A ledger says per compute saver how much it
saves and how much error it adds, at L6 and in production, and whether each production number is still
current.

**Architecture:** `tools/bundle/` (grown from `tools/golden/golden.py`: `build`, `run`, `digest` are reused)
defines the cases as a small dependency graph and writes `bundle.json`. `tests/bundle/` reads it through one
fixture. `tools/ledger/` pairs production runs by source (campaign `build_info`/`same_source`), measures
cost and error, and uses bundles on old and new commits to mark rows current or stale.

**Tech stack:** Python 3.11 (uv), pytest, Basilisk qcc at /oscar/data/dharri15/eaguerov/basilisk/src, SLURM.

**Spec:** `internal/docs/specs/2026-11-13-feature-guardrails-design.md` (read it first; section numbers below refer to it).

## Global constraints

- Never compute on a login node. Bundle and mutant runs only through `sbatch` (batch partition; exclude
  node2356,node2410; never mbessa-condo; MPI case not on batch2). `OMP_NUM_THREADS=1` for every run.
- `pytest` never starts a simulation or a job. Bundle tests run as
  `MFB_BUNDLE=<work>/bundle.json uv run python -m pytest tests/bundle` (no `-m`).
- Pass bars are written in `tools/bundle/bars.yaml` and committed **before** the first full bundle run
  (anti-goal: bars after results). The only exception is the `cold`/`cold2` noise floor, which the bars cite.
- Basilisk: no bare `return;` in an event; cs/fs recomputed after restore; mutants change a source *copy*,
  never `src/`.
- Diary entry for every run result and every power result; commit + push after each task.
- Work on branch `lean/phase2` in `/oscar/data/dharri15/eaguerov/Github/mfb-lean`; nothing merges while
  chains run from the campaign-2d checkout.

## Review focus

1. Bundle file missing, or built from another source than the working tree: one visible **failure**, never
   N skips (Task 3, `test_stale_bundle_fails`, `test_missing_bundle_fails`).
2. Bundle tests not requested (no `MFB_BUNDLE`): one module-level skip with the reason, so the default unit
   run on a laptop or CI stays green (Task 3, `test_not_requested_skips_once`).
3. A case crashes in the middle of the graph: `bundle.json` records the error, its children are recorded
   as "not run (parent X failed)", and their tests fail with that text (Task 2, `test_failed_parent_marks_children`).
4. A production run whose binary has no sidecar and an unknown SHA-256: excluded from the ledger and
   counted in its "provenance unknown" line (Task 6, `test_unknown_binary_is_counted_not_used`).
5. A pair whose two runs come from different sources: listed as "confounded", never used for a number
   (Task 6, `test_cross_source_pair_is_confounded`).

---

### Task 1: Case graph and parameters (pure Python)

**Files:** Create `tools/bundle/__init__.py`, `tools/bundle/cases.py`; test `tests/test_bundle_cases.py`.

**Interfaces — produces:**
- `@dataclass(frozen=True) class Case: name: str; build: str; parent: str | None; restart: str | None; level: int; ranks: int`
  (`restart`: which file of the parent to restart from: `"checkpoint.dump"`, `"replay_start.dump"`, or None)
- `BUILDS: dict[str, list[str]]` = golden's five + `"crosslevel": ["-DCROSS_LEVEL_WARMSTART=1"]` + `"mpi": ["-D_MPI=1", "-DEXTRA_TRACERS=0"]`
- `CASES: dict[str, Case]` with exactly: `cold, cold2, half, cont, kill, kill_resume, replay, loop, lean, l5, xl6, caseB_cold, caseB_warm, mpi4` (spec 4.1 table plus `half`, `kill_resume`)
- `order(cases) -> list[str]` (parents before children; raises `ValueError` on a cycle or an unknown parent)
- `params(name: str, T: float) -> dict` (golden's `base_params` + per-case keys)

Case parameters (decided here; T = the solver's period from `scripts.postprocess._t_scales`):
- `cold`: build `record`, L6, `t_end = 3T`, `n_mix_cycles = 1`, `n_rec_periods = 1`, `rec_dir = "."`. (One run is both
  the reference and the flux record. Task 2 checks that `record` and `default` builds give bitwise-equal
  flow outputs, from the golden gate's existing evidence `base` vs `record`; if they do not, `cold` uses
  `default` and a separate `record` case is added.)
- `cold2`: identical to `cold` (noise floor; expected bitwise except `logstats.dat`).
- `half`: as `cold` but `t_end = 2T`; the solver dumps `checkpoint.dump` at t_end. Its outputs over [0, 2T]
  must equal `cold`'s (a free guard: stopping earlier does not change the history).
- `cont`: parent `half`, restart `half/checkpoint.dump` with `restart_continue = 1, t_checkpoint = 2T,
  n_mix_cycles = 0, t_end = T` (golden's `restart` keys); compared with `cold` over (2T, 3T].
- `kill`: build `default`, as `cold` but `BIOREACTOR_MAXRUNTIME_S` set so the early dump fires at about
  1.5T (value set in Task 2 from `cold`'s measured wall time); then `kill_resume` restarts its
  dump with `restart_continue = 1` to 3T.
- `replay`: build `replay`, restart `cold/replay_start.dump`, golden `replay` keys with `rec_dir = <cold>`.
- `loop`: as `replay` with the record looped for 2 periods (`scripts.replay_lf.replay_params` keys).
- `lean`: build `lean`, as `cold` without record keys, `t_end = 3T`.
- `l5`: build `default`, `fidelity = 5`, `t_end = 2T`. `xl6`: build `crosslevel`, restart `l5/checkpoint.dump`
  at `fidelity = 6` with `t_checkpoint = 2T, omega_b_prev, theta_max_prev` (as tests/verification/test_cross_level_warmstart.py:102-113), `t_end = T`.
- `caseB_cold`: as `cold` at 25 rpm. `caseB_warm`: restart `half/checkpoint.dump` at 25 rpm with
  `omega_b_prev`/`theta_max_prev` of 32.5 rpm (chain.py:207-213 keys), `n_transition_cycles = 1`, release 1 period later.
- `mpi4`: build `mpi`, `ranks = 4`, as `lean` (compare with `lean`, same flags but MPI).

- [ ] Step 1: tests `test_order_puts_parents_first`, `test_cycle_raises`, `test_unknown_parent_raises`,
  `test_every_case_has_a_known_build`, `test_restart_cases_set_t_checkpoint` (every case with `restart`
  has `t_checkpoint` and `restart_continue` or the `*_prev` keys), `test_caseB_warm_rescales_from_32p5`.
- [ ] Step 2: run `uv run python -m pytest tests/test_bundle_cases.py` in a batch job (unit.sbatch pattern): FAIL (no module).
- [ ] Step 3: implement `tools/bundle/cases.py`.
- [ ] Step 4: same command: PASS.
- [ ] Step 5: commit "bundle: case graph and parameters".

### Task 2: Bundle runner and `bundle.json`

**Files:** Create `tools/bundle/bundle.py`, `tools/bundle/bundle.sbatch`, `tools/bundle/README.md`;
test `tests/test_bundle_runner.py`. Reuse `tools/golden/golden.py:build, run, digest` by import.

**Interfaces:**
- Consumes: Task 1 `CASES, BUILDS, order, params`.
- Produces: `bundle.py run SRC WORK [--only CASE ...]` and `bundle.json` schema:
  `{"src_hash": str, "git_commit": str|None, "builds": {name: {"flags": [..], "sha256": str}},
    "cases": {name: {"status": "ok"|"failed"|"not run", "error": str, "wall_s": float, "dir": str,
                     "files": {relpath: sha256}, "qoi": {name: float}}}}`
- `src_hash(src: Path) -> str`: SHA-256 over sorted `src/*.c, src/*.h` contents (identity without git).
- `qoi(run_dir: Path) -> dict`: `scripts.postprocess` results for the run (tau_95, edr mean, kLa_1T_10,
  dtmix_0.50 where defined; missing ones absent, not NaN).
- Runner contract: a failed case is recorded, its descendants become `"not run"` with `error = "parent <X> failed"`,
  the runner continues with independent cases and exits 1 at the end.

- [ ] Step 1: tests with a fake runner function injected (`run_case=` parameter): `test_failed_parent_marks_children`,
  `test_src_hash_ignores_file_order`, `test_bundle_json_schema` (keys above).
- [ ] Step 2: FAIL. Step 3: implement. Step 4: PASS (batch job).
- [ ] Step 5: smoke on a compute node: `sbatch tools/bundle/bundle.sbatch $PWD/src /oscar/scratch/eaguerov/bundle/head --only cold cold2`;
  expected: both `ok`, `cold` and `cold2` file hashes equal except `logstats.dat`. Record wall time per case in the diary.
- [ ] Step 6: level check (spec 4.1, L4 rejection): add `--fidelity` override, run `cold` at L4 and L6 once;
  diary: dtmix_0.50 and kLa_1T_10 defined (not NaN) at L6, and whether at L4. Keep L6 unless L6 fails this.
- [ ] Step 7: commit "bundle: runner, bundle.json, smoke + level check".

### Task 3: Bars file and the bundle fixture (before any full run)

**Files:** Create `tools/bundle/bars.yaml`, `tests/bundle/__init__.py`, `tests/bundle/conftest.py`;
test `tests/test_bundle_fixture.py`.

**Interfaces — produces:** fixture `bundle` (session) -> dict (parsed bundle.json); fixture `bars` -> dict;
helper `case(bundle, name) -> dict` that fails with the recorded error if status != "ok".

Bars (copy values, cite the source in a comment next to each):
- `cold` vs `cold2`: bitwise, all files but `logstats.dat` (golden gate 2026-11-12).
- `lean` vs `cold`: bitwise for `normf.dat`, `shear_stress.dat`, `vol_frac_interf.dat`, `edr_pct.dat`, `tr_oxy.dat`
  oxygen columns. Pre-registered as bitwise; a failure is a finding to explain, not a reason to loosen.
- `half` vs `cold` over [0, 2T]: bitwise on the common rows. `cont` vs `cold` over (2T, 3T], `kill_resume` vs `cold`, `mpi4` vs `lean`: the tolerances of
  tests/verification/test_restart_continue.py and test_mpi_checkpoint_parity.py (copy the numbers, cite file:line).
- `replay` vs `cold`: the pre-registered F1 criteria of tests/test_replay_gate.py (copy, cite).
- `xl6` vs `cold` after T: tests/verification/test_cross_level_warmstart.py bars (copy, cite).
- `caseB_warm` vs `caseB_cold`: 5 % on tau_95 and EDR mean over the last period (diary x8 bar, 2026-10).
- `loop`: no pass bar (ledger only); guard = finite C* and it ran.

- [ ] Step 1: tests: `test_not_requested_skips_once` (no MFB_BUNDLE → collected bundle tests report 1 skip
  reason), `test_missing_bundle_fails`, `test_stale_bundle_fails` (bundle src_hash != `src_hash(src/)` → failure
  naming both hashes), `test_failed_case_fails_with_its_error`.
- [ ] Step 2: FAIL. Step 3: implement conftest + write bars.yaml. Step 4: PASS.
- [ ] Step 5: commit "bundle: pre-registered bars and fixture" (this commit must precede Task 4's run).

### Task 4: Guards and the first full bundle

**Files:** Create `tests/bundle/test_restart.py` (S5, S9, S10, S12), `tests/bundle/test_replay.py` (S1, S2),
`tests/bundle/test_warm.py` (S3, S4, S8, S15), `tests/bundle/test_builds.py` (S6, S7), `tests/bundle/test_averaging.py` (S11).
Each test compares the cases of spec 4.1 with the bar from `bars.yaml`; S9 also reads `restart_diagnostic_*.txt`
(fs sanity line), S12 checks the last time in `normf.dat` equals `periods.py`'s stop time, S11 checks the
liquid area in `shear_stress.dat` is below the bag area (not half the domain).

- [ ] Step 1: write the tests (they fail without a bundle: run with `MFB_BUNDLE` pointing at a missing file → 1 failure).
- [ ] Step 2: full run on a node: `sbatch tools/bundle/bundle.sbatch $PWD/src /oscar/scratch/eaguerov/bundle/head` (+ MPI node job for `mpi4`).
- [ ] Step 3: `MFB_BUNDLE=/oscar/scratch/eaguerov/bundle/head/bundle.json uv run python -m pytest tests/bundle` in a batch job.
  Every failure is investigated as a possible real defect first (systematic-debugging), never by changing a bar.
- [ ] Step 4: diary (results per guard, wall time per case, total core-h); commit "bundle: guards S1-S15, first bundle".

### Task 5: Power — planted defects

**Files:** Create `tools/bundle/mutants.py`, `tools/bundle/mutants.sbatch`.

**Interfaces:** `MUTANTS: dict[str, Mutant]`, `Mutant(file: str, old: str, new: str, cases: list[str], must_fail: list[str])`
(exact-string patch of a `src/` copy; `old` must occur once). `mutants.py run NAME WORK`: copy src, patch, run
`bundle.py run --only <cases + parents>`, run pytest on `must_fail` tests, exit 0 iff each named test fails.

Mutants (the builder finds the exact `old` string in the named file):
- `no_fs_recompute` (init.h, restore_geometry: skip the cs/fs recompute) → S9 test fails.
- `drop_tracer_dump` (the dump list in checkpoint.h: drop `c2`) → S10 test fails.
- `replay_off_by_one` (replay.h: read the record one step late) → S1 test fails.
- `no_phase_sidecar` (checkpoint.h: do not write `checkpoint.phase`) → S5 kill test fails.
- `prerefine_old` (flag `-DCROSS_LEVEL_PREREFINE_SOLID=0` on `crosslevel`) → S15 test fails.
- `no_rescale` (restart.h/init.h: skip the velocity rescale when omega changes) → S4 test fails.
- `lean_drops_oxygen_mask` is NOT used (no realistic defect found yet); S6's power = bitwise bar (any change fails).

- [ ] Step 1: one mutant end-to-end in a batch job; expected: exit 0 and the named test FAILED in the log.
- [ ] Step 2: all mutants; a mutant whose guard passes = a guard with no power: fix the guard (not the mutant), rerun.
- [ ] Step 3: diary table mutant → guard → failed yes/no; commit "bundle: power of every guard".

### Task 6: Ledger provenance and pairs

**Files:** Create `tools/ledger/__init__.py`, `tools/ledger/provenance.py`, `tools/ledger/pairs.yaml`;
test `tests/test_ledger_provenance.py`.

**Interfaces:**
- Consumes: `scripts.campaign2d.campaign.build_info`, `same_source`.
- `provenance(run_dir: Path, known: dict[str, dict]) -> dict` → `{"commit", "flags", "status": "sidecar"|"sha-match"|"unknown"}`
  (`known`: sha256 → sidecar dict, built from every `*.build.json` under /oscar/data/dharri15/eaguerov/bin).
- `pair_status(a: dict, b: dict, mechanism_flag: str | None) -> "same-source"|"confounded"`.
- `pairs.yaml`: per mechanism, list of `{with: run, without: run, qoi: [...]}` from the diary's pairs
  (xlevel vs l10c; w9_l9 vs cold L9 at 17.5/25/32.5; F3b record vs replay; x8e continuation; lean vs full
  wall times; timeouts recovered by autoextend). Run names copied from the diary entries cited in the file.

- [ ] Step 1: tests `test_sidecar_wins`, `test_sha_match_without_sidecar`, `test_unknown_binary_is_counted_not_used`,
  `test_cross_source_pair_is_confounded`, `test_cross_level_flag_only_is_same_source`.
- [ ] Step 2: FAIL. Step 3: implement. Step 4: PASS. Step 5: commit.

### Task 7: Ledger numbers, carry-over, figure

**Files:** Create `tools/ledger/ledger.py`, `tools/ledger/carry.py`, `scripts/figures/fig_feature_ledger.py`;
output `experiments/feature_ledger.json`; test `tests/test_ledger_numbers.py`.

**Interfaces:**
- `cost_core_h(run_dir) -> float` (wall time from logstats.dat or sacct × ranks); `rel_error(a, b, qoi) -> float`.
- `ledger.py build` → `feature_ledger.json`: per mechanism row `{saves_core_h, error: {qoi: rel}, L6: {...},
  status: "current"|"stale"|"confounded", source_commit, runs: [...]}`.
- `carry.py COMMIT_BUNDLE HEAD_BUNDLE MECHANISM -> "current"|"stale: <first differing output>"`: compares the
  mechanism's guard outputs between the two bundles with the Task 3 bars.

- [ ] Step 1: unit tests on synthetic run folders for `cost_core_h`, `rel_error`, and carry (equal → current,
  one changed file → stale naming it).
- [ ] Step 2: FAIL → implement → PASS.
- [ ] Step 3: bundle on each production source commit in the ledger (start with 2ef5b07): `git worktree add` the
  commit, `sbatch tools/bundle/bundle.sbatch <worktree>/src /oscar/scratch/eaguerov/bundle/<commit>`; cases the
  old source cannot run are recorded "not run" and their rows stay "stale (no old bundle)".
- [ ] Step 4: `ledger.py build`, figure (figstyle; colour = level, marker = QoI; no text in the figure).
  Rows whose numbers need a new production pair: list with core-h estimate → AskUserQuestion before any run.
- [ ] Step 5: diary + commit "ledger: first build".

### Task 8: Retire duplicate compute, docs

**Files:** Modify the `medium` tests that the bundle now covers (test_restart_continue, test_mpi_checkpoint_parity,
test_cross_level_warmstart, test_graceful_timeout's run tests); create `docs_site/developer/guardrails.md` (+ mkdocs nav).

- [ ] Step 1: only for guards with shown power (Task 5): keep the CI-level (L3) check of each medium test,
  remove its duplicate assertions that the bundle covers at L6; note in each test's docstring which bundle test covers it.
- [ ] Step 2: docs page: how to run the bundle, the tests, the mutants; how to read the ledger. mkdocs --strict in batch.
- [ ] Step 3: unit suite in batch (expect the same pass count as before minus removed assertions, 0 new skips); diary; commit; push.
