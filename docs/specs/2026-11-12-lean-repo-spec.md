# Spec: a lean, cold-reader-friendly repository

Status: DRAFT for review (2026-11-12). Scope: `multi-fidelity-bioreactor` only (BioReactor3D is not touched).
Internal document: lives in `docs/` (never published; `mkdocs.yml` builds `docs_site/` only).

## 1. Goal

A person who has never seen this repo (a new student, a reviewer, Elvis in two years) can, in one afternoon:

1. understand what the code simulates and why (README, one page);
2. build and run one simulation, and rebuild one figure, with one command each;
3. find the code for any feature in at most two hops (repo map -> folder README -> file);
4. learn every hard-won lesson (a bug that cost days) without reading the git history.

The repo will outlast the DFD talk. Nothing is removed because "DFD does not need it"; things are removed only
because nothing uses them AND their lesson is kept elsewhere.

## 2. Principles (apply to every change)

- P1. No behaviour change. Every code move is checked by an experiment (section 7) before it is merged.
- P2. Lessons are kept, not deleted. A long comment that tells history moves to the lessons page; the code keeps a
  one-line pointer to it (`// LESSON fs-after-restore: see docs_site/.../pitfalls.md#fs-after-restore`).
- P3. Strict boundary between public and internal text (section 6). Internal text never enters the mkdocs nav.
- P4. Nothing is lost. Before any deletion: grep all consumers (also OUTSIDE the repo: the DFD deck `talks/dfd2026/tools`
  imports `scripts.postprocess`, `scripts.bo2d`, `scripts.campaign2d`, `scripts.hydro_dataset_v2`), then tag
  `archive/pre-lean-2026-11` so the old tree is one checkout away.
- P5. Follow Basilisk examples for all C changes (project CLAUDE.md). Splitting a Basilisk file into local headers
  that define events is the canonical pattern (e.g. `two-phase.h`, `tension.h`).
- P6. Small steps, one pull request each, each with its own check. No "big bang" reorganisation.

## 3. Measured state (2026-11-12)

| Area | Now | Problem |
|---|---|---|
| tracked files | 1,298; `.git` 150 MB | big media in history (largest blobs: GIF 3.3 MB, MP4s 1.6-2.0 MB, one JSON 2.8 MB) |
| `src/BioReactor.c` | 3,759 lines (Kim upstream: 830); 1,554 comment lines (41 %); 120 `#if`s over 30 switches | `main()` ~450 lines, `event init` ~660, `event normcal` ~480, record/replay ~460 |
| `src/draw3.h`, `view3.h` | 1,836 + 749 lines | copies of canonical `draw.h`/`view.h`; `draw3.h` differs by ~8 lines only |
| `scripts/` | 607 files, 59k lines; 422 `.py` files at top level | 71 `submit_*`, 15 `resubmit_*`, 52 `test_*` (outside `tests/`), 45 `plot_*`, 32 `diag_*`; 111 top-level scripts named by no other file |
| `experiments/` | 397 files, 225k lines (137 JSON, 134 PNG, 64 PDF, 7 MP4) | results, figures and videos mixed; most figures can be rebuilt |
| `config/` | 34 YAML; 19 named by no other file | smoke/diag one-offs next to the canonical configs |
| docs | `docs/` 73 files (internal); `docs_site/` 42 files (public mkdocs) | site unchanged for 30+ days; media duplicated between `docs/` and `experiments/` |
| root | `diary.md` 2,053 lines + `diary_archive_to_2026-10-07.md` 9,368 lines; `runs` symlink | archive at root |

Audit script and full table: `scratchpad repo_audit.py` (to be committed as `tools/repo_audit.py`, see W8).

## 4. Target layout

```
README.md            what it is, 5-minute quick start, repo map (links to every folder README)
src/                 Basilisk solver: a short driver + one local header per concern (W1)
  BioReactor.c       driver: includes, parameters, main(); target < 500 lines
  geometry.h  forcing.h  restart.h  replay.h  init.h  diagnostics.h  outputs.h  ...
  README.md          event order, compile switches table, how to add a feature
mfb/                 importable Python library (postprocess, campaign, bo, mfbml_local, figstyle)
scripts/             entry points only, grouped: run/  analyse/  figures/  admin/ ; each folder has a README
tests/               ALL tests (unit, medium, hpc); no test_*.py elsewhere
config/              canonical configs + examples/ ; one-offs archived
experiments/         small results only (JSON/CSV) + a manifest; figures rebuilt by `make figures`
docs_site/           PUBLIC: tutorials, how-to, reference, explanation, developer pitfalls (lessons)
docs/                INTERNAL: specs, plans, literature reviews, meetings, diary archive  (README says "internal")
```

## 5. Workstreams

### W1. `src/BioReactor.c` -> driver + local headers (largest value, highest risk)

Split by concern, in this order (lowest risk first). Each header keeps its events; include order in the driver is
kept identical, because Basilisk runs same-named events in registration order (the file says so at the replay
event, line ~3446) and the order of `event vof`/`event init` actions matters.

| Header | Content now | Lines (approx.) |
|---|---|---|
| `outputs.h` | snapshots, uv field, phase sidecar, movies (`movies_output`, `movies_output_tau`) | 3542-3760, 2661-2800 |
| `diagnostics.h` | `normcal`, `logstats`, `pressure_diagnostics`, `streaming_vorticity` | 2880-3440, 2326-2436 |
| `replay.h` | flux record / replay subsystem + its event | 416-877, 3443-3540 |
| `restart.h` | checkpoint, phase-offset logic, `assert_fs_sane`, restart diagnostics | from `main()` ~1000-1100, 1341-1455, 2800-2870 |
| `forcing.h` | `acceleration` event, rocking motion, kick | 2310-2660 |
| `scalars.h` | tracer and oxygen setup | 2148-2308 |
| `init.h` | `event init` (split into named static functions inside: geometry, fields, restore) | 1459-2123 |

Then, inside `main()`: move parameter derivation into `params_derive()` (in `params_read.h` or `geometry.h`).

Compile switches (30). For each switch: (a) search builds, configs, tests, run logs and `runs/*/params.json`
compile records; (b) classify LIVE / RESEARCH (used once, kept for reproducibility of a diary result) / DEAD.
DEAD branches are removed after their lesson is written down. Candidates with no reference outside the C file
(rough search): `REC_CHECKSUM`, `STREAMING_VORTICITY`, `KICK_TEST`, `KICK_EPS`. Record/replay stays (decided).
`src/README.md` gets a table: switch, default, what it does, who uses it.

Stale headers: `draw3.h`/`view3.h` are copies of canonical `draw.h`/`view.h` (project rule: check `git log`, diff
function by function, use canonical + restore project additions only). Expected: both files go, ~2,600 lines.

Comments: history and "this cost four days" text -> lessons page (W6). The code keeps the WHY in one or two lines
plus the pointer. Target comment share ~20 %, not lower: comments are wanted.

### W2. `scripts/` -> library + grouped entry points

1. Separate the importable library (modules imported by other code: `postprocess`, `hydro_dataset_v2`, `bo2d`,
   `campaign2d`, `mfbml_local`, `figstyle`, ...) from scripts you run. Library moves to a package; the old import
   paths keep thin re-export shims until the deck and every caller have moved (P4). No shim removal before the talk.
2. Entry points go to `scripts/run/`, `scripts/analyse/`, `scripts/figures/`, `scripts/admin/`, each with a README
   (one line per script: what it makes, which input).
3. One-off scripts of finished campaigns (`submit_*`, `resubmit_*`, `diag_*` for closed bugs): archived (tag +
   delete) unless a diary entry or a figure still needs them; the diary entry gets the tag name.
4. `scripts/test_*.py` (52): ~~pytest tests -> `tests/`~~ **Amended 2026-11-13:** these are numerical
   experiments (one per diary hypothesis; many submit SLURM jobs), not pytest tests. They stay out of `tests/`
   (pytest would collect them). They keep their names because the diary cites them.

### W3. `tests/`

One test tree. Mark levels with the existing markers (unit / medium / hpc / slow). Add the golden-run test of
section 7 as `medium` (CI only, per project rule: no `pytest -m` on OSCAR by hand).

### W4. `config/`

Keep canonical configs (`bo_config.yaml`, `chain_config.yaml`, ...) and move smoke examples to `config/examples/`.
The 19 configs named by nothing: archive unless a diary result cites them.

### W5. `experiments/`, figures and media

- Keep in git: small results (JSON/CSV) that figures read, plus `experiments/MANIFEST.md` (result -> script -> figure).
- Out of git (new commits only; history NOT rewritten): MP4, GIF and rebuildable PNG/PDF. Store under
  `/oscar/data/dharri15/eaguerov/mfb_artifacts/` with an `rclone copy` to Drive. `.gitignore` blocks new media.
- `make figures` rebuilds every repo figure from the kept results, with `scripts/figstyle.py` (one visual grammar).
- Remove duplicates between `docs/` and `experiments/` (e.g. the same `08_l8_matrix...mp4` in both).

### W6. Documentation (strict public/internal boundary)

Public (`docs_site/`, mkdocs, built with `--strict` in CI):
tutorial (build, run, postprocess one case), how-to (submit a sweep, chain a run, rebuild figures), reference
(params.json fields, compile switches, outputs and units), explanation (physics, numerics, multi-fidelity method),
**developer pitfalls** (technical lessons written neutrally: face fields not dumped, bare `return` in events,
volume vs wall metrics, ...). Refresh the stale pages against the current code.

Internal (`docs/`, never published): specs, plans, literature reviews, meetings, diary archive, project history
(dates, costs, who found what). `docs/README.md` says "internal". A CI check fails when a `docs_site/` page links
into `docs/` or when `mkdocs.yml` nav names a `docs/` file.

Lesson routing rule: the technical content of a lesson -> public pitfalls page; the story (dates, cost, dead ends)
-> internal history page that links to the public one.

### W7. Root and repo map

Root keeps: README, LICENSE, Makefile, pyproject/uv files, mkdocs.yml, diary.md, BACKLOG.md. The diary archive moves
to `docs/diary/` (fix its one referrer). README gains a repo map and a "where do I start" table by role
(run simulations / change the solver / rebuild a figure / extend the surrogate).

### W8. Keep it lean (so it does not grow back)

`tools/repo_audit.py` in CI as a report (not a gate): unreferenced files, top-level script count, largest blobs.
Pre-commit: block new media files over 1 MB outside an allow-list.

## 6. Public vs internal: what goes where

| Content | Public `docs_site/` | Internal `docs/` |
|---|---|---|
| how to run, params, outputs | yes | no |
| physics and method explanation | yes | no |
| technical pitfalls (neutral) | yes | no |
| bug stories, dates, compute cost, dead ends | no | yes |
| plans, specs, meetings, PI notes, literature notes | no | yes |
| talk material | no | yes (or outside the repo) |

## 7. Verification gate (applies to W1, and to W2 for figure scripts)

- Golden run: L6, short window (e.g. 2 rocking cycles), on a compute node (sbatch), for every LIVE switch
  combination (default; TRACER+OXYGEN; RECORD_FLUX; REPLAY_FLUX; restart from a checkpoint).
- Pass = every output file byte-identical to the pre-change binary, plus the dump. If not byte-identical (e.g. a
  reordering that changes floating-point sums), the step is reverted, not "accepted within tolerance".
- Test power: the golden run must FAIL when an event is removed or reordered (checked once by a deliberate
  mutation before the gate is trusted).
- Figure scripts: pixel-identical output (as done for the deck fit caches).

## 8. Phases

| Phase | Work | Risk | When |
|---|---|---|---|
| 0 | golden-run harness + its mutation test; `tools/repo_audit.py`; archive tag | none | now |
| 1 | W1 headers in table order; stale `draw3.h`/`view3.h`; switch audit | medium (gated) | now; no binary deploy to live campaigns without the gate |
| 2 | W3, W4, W7, W8 | low | now |
| 3 | W6 docs split and refresh; lessons pages | low | now |
| 4 | W2 library + script regrouping (shims kept) | medium (imports) | after the DFD talk for anything the deck imports |
| 5 | W5 media out of git | low | after phase 4 |
| 6 | remove shims | low | after the deck is archived |

Done when: `BioReactor.c` < 500 lines and no header > 700; top-level `scripts/` < 15 files; no `test_*.py` outside
`tests/`; every folder has a README; mkdocs strict build passes; repo audit shows 0 unreferenced scripts outside
`archive`; golden runs byte-identical.

## 9. Decisions (user, 2026-11-12)

1. Library: `mfb/` package at the repo root (`from mfb.postprocess import ...`); old paths keep re-export shims.
2. Lessons: split. Technical content -> public "developer pitfalls" page (no dates, no costs); the story -> internal
   history page; the two link to each other.
3. `talks/dfd2026/` in the repo: no change before the talk; after it, the final deck goes to an internal archive
   outside the repo and the repo keeps a link.
4. Archive: git tag `archive/pre-lean-2026-11` only, then delete; each diary entry that cites a deleted file gets the
   tag name. No `archive/` folder.
5. Record/replay stays in the solver (the repo outlives the DFD talk).
6. Results: keep small results in git, move big media out, no history rewrite.
