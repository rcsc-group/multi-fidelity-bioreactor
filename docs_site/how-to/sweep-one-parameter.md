# How to sweep one parameter with checkpoint restart

Use this when you want to sweep **one** parameter (e.g. rocking frequency
`omega_b`) across several values, chained via checkpoint restart so each
segment after the first starts from an already-developed flow instead of
cold. For why this is worth doing and how the restart ramp works, see
[Checkpoint restart and warm-start chains](../explanation/checkpoint-restart.md).

!!! tip
    To sweep multiple parameters at once, or write configs in JSON instead
    of YAML, use [Sweep any parameter combination](sweep-json-multi-param.md) instead.

## 1. Write a chain config

Copy `config/chain_config.yaml` and edit:

```yaml
# ── Fixed params (same for all segments) ────────────────────────────────────
fidelity: 7
geometry: {a: 0.20, b: 0.09, n: 2.0}   # superellipse; n=2 → ellipse
fill_level: 0.4

# ── Timing ───────────────────────────────────────────────────────────────────
n_mix_cycles: 80         # rocking cycles before O₂ injection (segment 0)
n_transition_cycles: 10  # rocking cycles before O₂ re-injection (segments 1+)
t_buffer: 150            # non-dim kLa measurement window

# ── Sweep: one segment per value ─────────────────────────────────────────────
sweep:
  parameter: omega_b
  values: [3.14159, 6.28318]   # 0.5 Hz, 1.0 Hz

# ── Base motion (non-swept params, fixed across all segments) ─────────────────
motion:
  n_harmonics: 1
  theta_max:      [5.0, 0.0, 0.0]
  phi_angular:    [0.0, 0.0, 0.0]   # index 0 must stay 0.0
  omega_h:        6.28
  amplitude_h:    [0.02, 0.0, 0.0]
  phi_horizontal: [0.3, 0.0, 0.0]

videos: false          # set true to use BioReactor-video and render MP4s
submit: true
walltime: "04:00:00"
```

**Supported sweep parameters:**

| Name | Maps to |
|------|---------|
| `omega_b` | rocking angular frequency (rad/s) |
| `omega_h` | horizontal translation frequency (rad/s) |
| `theta_max_0` / `_1` / `_2` | `theta_max[0..2]` |
| `amplitude_h_0` / `_1` / `_2` | `amplitude_h[0..2]` |
| `phi_angular_1` / `_2` | `phi_angular[1..2]` (index 0 is always 0) |
| `phi_horizontal_0` / `_1` / `_2` | `phi_horizontal[0..2]` |

## 2. Smoke-test before submitting long jobs

```bash
uv run python scripts/chain.py config/chain_config_smoke.yaml
```

Then check (see [Your first sweep](../tutorials/first-sweep.md) for what this looks like in full):

- `runs/<seg1_id>/logstats.dat` starts at `t > 0`, not `t=0`
- `runs/<seg1_id>/results.json` has finite kLa values

## 3. Submit for real

```bash
python scripts/chain.py config/chain_config.yaml
```

Prints one line per segment with run IDs and SLURM job IDs. SLURM
`--dependency=afterok` ensures each segment starts only after the previous
one succeeds.

## 4. Collect results

Each segment writes independently to `runs/<run_id>/results.json` — there's
no aggregation script for Workflow B specifically; collect manually or use
`collect_results.py` (see [Scripts reference](../reference/scripts.md)).
