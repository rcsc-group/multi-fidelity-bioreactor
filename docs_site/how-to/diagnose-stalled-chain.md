# How to diagnose a stalled or failed chain

A chain that isn't progressing usually falls into one of a small number of
known failure modes. Check these in order before assuming something novel
is wrong.

## 1. Is the job actually still running?

```bash
squeue -u $USER
```

If nothing's there and the chain isn't done, either the job finished and
failed to auto-submit its successor, or it was never queued. Check the
SLURM log:

```bash
tail -50 logs/slurm_<jobid>.out
```

Look for the self-chaining line the template prints on success:
`Submitted next segment: <run_id> (job <jobid>)`. If it's missing, the job
either crashed before reaching that point, or the next segment's
`params.json` wasn't staged where the template expected it — check
`_experiment_dir`/`_canonical_run_dir` in the failed segment's own
`params.json` actually point at this project's `runs/` directory, not a
stale path from a different clone or a moved repo.

## 2. `checkpoint.dump` existing does NOT mean the segment finished

This is the single most common false positive. A job killed seconds into a
restart still has a `checkpoint.dump` on disk — it's just an untouched copy
of *that segment's own input*, not its output. The only reliable proof a
segment actually completed is a well-formed `results.json`, written
exclusively by `postprocess.py` at the very end of a real run.

If you're staging the next segment's checkpoint by hand instead of letting
the chain self-submit, use `stage_segment.py` rather than a manual `cp` —
it refuses to stage from a predecessor that doesn't have a genuine
`results.json`:

```bash
uv run python scripts/stage_segment.py <prev_run_id> <next_run_id>
```

## 3. Restart segment starting at `t=0` instead of `t>0`

```bash
head -1 runs/<seg1_id>/logstats.dat
```

If this shows `t: 0`, the restore didn't happen — the binary was given a
`t_checkpoint` of 0 (fresh-run path), or the wrong `checkpoint.dump` was
staged. Check `t_checkpoint` in that segment's `params.json` matches what
the predecessor actually reached (its own `t_dump_checkpoint`, visible as
the argument to the `checkpoint: writing checkpoint.dump at t=...` line in
its SLURM log).

## 4. Job runs forever, never reaches `t_end`

Basilisk's event system won't exit `run()` on its own if every registered
event is unconditional (e.g. `acceleration(i++)`) — this bit the project
once already. The fix is that `dump_checkpoint`'s event handler explicitly
returns `1` (Basilisk's `event_stop` idiom) once it's written the final
checkpoint. If you've modified `BioReactor.c`'s event handlers and a job
that should finish in minutes is still running hours later with no
progress in `logstats.dat`, check that whichever event is meant to end the
run still returns `1`.

## 5. `results.json` has NaN kLa values

Means `tr_oxy.dat` has no data, which happens when `t_end < t_mix` — the run
ended before oxygen injection ever started. Check `n_mix_cycles` and
`t_end` are consistent for the fidelity/frequency you're using; see
[params.json reference](../reference/params.md) and the note on `t_buffer`
sizing in the [Glossary](../glossary.md).

## Still stuck?

`tau_chain_watchdog.py` and `theta_chain_watchdog.py` implement exactly this
diagnosis loop programmatically (self-healing: they detect a completed
segment with no queued successor and resubmit it). Reading their source is
often faster than re-deriving the check by hand.
