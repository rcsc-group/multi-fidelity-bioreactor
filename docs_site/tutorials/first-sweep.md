# Your first sweep

This chains two fidelity-3 runs together with checkpoint restart — the
mechanism every real sweep in this project is built on. See
[Checkpoint restart and warm-start chains](../explanation/checkpoint-restart.md)
for why it works the way it does; this page is just running it.

## Run the smoke-test chain

```bash
uv run python scripts/chain.py config/chain_config_smoke.yaml
```

Expected output:

```
  [seg 0] run=46e365e1  omega_b=3.14  n_mix=5  t_end≈12.7  → job 4153341
  [seg 1] run=bedb814c  omega_b=6.28  n_mix=3  t_end≈11.6  → job 4153342

Chain submitted: ['4153341', '4153342']
```

Two SLURM jobs, `--dependency=afterok` between them — segment 1 won't start
until segment 0 finishes successfully. Both are fidelity 3, so this whole
chain finishes in under a minute of actual compute once it starts running.

## Verify the restart actually happened

This is the check worth internalizing, because it's the one thing that tells
you checkpoint restart is doing its job rather than silently starting cold:

```bash
head -1 runs/46e365e1/logstats.dat   # segment 0 -- should start at t=0
head -1 runs/bedb814c/logstats.dat   # segment 1 -- should start at t>0
```

Ours:

```
$ head -1 runs/46e365e1/logstats.dat
i: 0 t: 0 dt: 0.0507282 #Cells: 64 ...

$ head -1 runs/bedb814c/logstats.dat
i: 2278 t: 13.2 dt: 0.00666667 #Cells: 64 ...
```

Segment 1 picks up at `t=13.2` — wherever segment 0's `checkpoint.dump` left
off — not at `t=0`. If you ever see a restart segment starting at `t=0`,
something is wrong with the chain (see
[Diagnose a stalled or failed chain](../how-to/diagnose-stalled-chain.md)).

## Check both segments produced results

```bash
uv run python scripts/postprocess.py runs/46e365e1/
uv run python scripts/postprocess.py runs/bedb814c/
```

Each writes its own `results.json` — a chain doesn't aggregate results for
you; see [Output files reference](../reference/output-files.md).

## Next

- [Sweep one parameter](../how-to/sweep-one-parameter.md) — do this for real, at production fidelity
- [Sweep any parameter combination](../how-to/sweep-json-multi-param.md) — when one swept parameter isn't enough
