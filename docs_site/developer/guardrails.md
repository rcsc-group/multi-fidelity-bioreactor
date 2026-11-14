# Feature guardrails and the feature ledger

The solver has several mechanisms that save compute: restarting from a dump, warm-starting a finer level
from a coarser one, warm-starting a new condition from an old one, recording the flow once and replaying
it for the scalars, a lean build without the extra tracers, MPI. Each one can break silently: the run still
finishes, and only the numbers are wrong. This page explains the three tools that guard them.

| Tool | Question it answers | Where |
|---|---|---|
| Feature bundle | Does each mechanism still work on this source? | `tools/bundle/`, `tests/bundle/` |
| Mutants | Would the guards notice if a mechanism broke? | `tools/bundle/mutants.py` |
| Feature ledger | What does each mechanism save, and what error does it add in production? | `tools/ledger/` |

## The feature bundle

One batch job builds every solver variant once (default, lean, record, replay, cross-level), then runs a
small chain of L6 cases in which each case restarts from another one's dump. For example, `half` stops at
4 periods and `cont` continues it, and `cold` is the unbroken run they must reproduce. The cases and their
parents are in `tools/bundle/cases.py`. The job writes `bundle.json`, with the status, error, file hashes and
quantities of every case.

```bash
# from the repository root, on the login node (the jobs run on compute nodes)
sbatch tools/bundle/bundle.sbatch     $PWD/src /oscar/scratch/$USER/bundle/mywork
sbatch tools/bundle/bundle_mpi.sbatch $PWD/src /oscar/scratch/$USER/bundle/mywork   # the 4-rank case
```

The serial job takes about one hour. The MPI case runs in its own job, because MPI is not available on
every partition, and it writes `mpi4.json` next to `bundle.json`.

If a case fails, the bundle records its error, and its children are recorded as "not run (parent X
failed)". The other cases still run.

## The guards

The guards are ordinary pytest tests that **only read** `bundle.json`. They never start a run:

```bash
MFB_BUNDLE=/oscar/scratch/$USER/bundle/mywork/bundle.json uv run python -m pytest tests/bundle
```

- Without `MFB_BUNDLE`, every guard is skipped with one reason, so the default unit run stays green.
- A bundle built from other sources than `src/` (its `src_hash` differs) is a failure, not a skip: a guard
  never passes on old evidence. After any change in `src/`, run the bundle again.

| Guard | Mechanism | Pass condition (bar) |
|---|---|---|
| S1, S2 | flux record/replay | exact replay equals the coupled run to 1e-6; a looped replay runs its whole window |
| S3, S15 | cross-level warm start (L5 → L6) | matches the cold L6 run |
| S4 | warm start to a new condition (32.5 → 25 rpm) | matches a cold run at the new condition within max(5 %, 2 sd), mean of the last 3 periods |
| S5 | same-condition restart, early dump at the wall-clock limit | reproduces the unbroken run |
| S6, S7 | lean build, MPI | lean equals default byte for byte; MPI within 5 % of serial |
| S8 | tracer release after a warm start | released after the transition cycle, not before |
| S9, S10 | restore | cs/fs recomputed; tracer and oxygen come back unchanged |
| S11 | averaging region | the liquid area, not half the domain |
| S12 | stop time | the run stops on the period boundary that `scripts/periods.py` predicts |
| S13 | restart label | the rocking angle and its rate do not jump at a restart; ambiguous restarts are refused |

Every bar is in `tools/bundle/bars.yaml`, with its `source`: the test or diary entry the number comes from.
The bars were written before the first bundle ran. A failing guard is a finding to explain, not a reason to
change its bar. The one bar that changed after a result (S4) records why in `window_source`. That change
came from an extension run whose pass rule was written before the run.

## Mutants: do the guards have power?

A guard that passes on broken code is worse than no guard. `tools/bundle/mutants.py` plants one realistic
defect in a copy of `src/` (never in `src/` itself), runs only the cases that the defect touches, and
checks that the named guards fail:

```bash
sbatch -J ramp_from_rest tools/bundle/mutants.sbatch ramp_from_rest /oscar/scratch/$USER/bundle/mut
```

| Mutant | Defect | Caught by |
|---|---|---|
| `no_fs_recompute` | no cs/fs rebuild after restore | S9 |
| `drop_tracer_dump` | tracer missing from the dump | S10 |
| `replay_off_by_one` | replayed flux one step late | S1 |
| `no_phase_sidecar` | dump written without its label | S5, S12 |
| `prerefine_old` | old pre-refine order at a level change | S3 |
| `ramp_from_rest` | restart ramps the rocking from zero | S13, S5 |
| `no_rescale` | no velocity rescale at an rpm change | S4 (see the diary for the latest result) |

The patches are exact strings, so a patch that no longer matches the source fails loudly
(`tests/test_bundle_mutants.py`). A mutant therefore cannot silently test nothing.

## The feature ledger

The bundle says whether a mechanism works at L6. The ledger says what it does in production, from pairs of
real runs at the same condition, one with the mechanism and one without (`tools/ledger/pairs.yaml`).

```bash
# in a batch job (it calls sacct)
uv run python -m tools.ledger.ledger --old OLD/bundle.json --new HEAD/bundle.json experiments/feature_ledger.json
```

Each row (one mechanism) has:

- `error`: the median relative change of each quantity, with / without − 1, over the clean pairs.
- `saves_core_h`: the median compute saved per pair, from sacct (all segments, time-outs and failures).
  This is only meaningful for mechanisms that replace compute (warm starts, replay).
- `status`:
    - `confounded`: no clean pair. The two runs of each pair came from different solver sources
      (`tools/ledger/provenance.py`), so their difference is not the mechanism's.
    - `current`: the production pairs came from an old source, and the mechanism's bundle outputs are byte
      for byte the same on that source and on HEAD (`tools/ledger/carry.py`). The production number still
      holds.
    - `stale: <case/file>`: the bundle outputs differ, starting at that file. The production number was
      measured on code that behaves differently now; measure it again before you rely on it.

To build the old-source bundle, check out the production commit in a worktree and run the bundle on its
`src/` (for example `git worktree add /oscar/scratch/$USER/tmp/wt_2ef5b07 2ef5b07`). Cases that the old
source cannot run are recorded as "not run", and rows that need them stay stale.

The figure `scripts/figures/fig_feature_ledger.py` plots the rows.
