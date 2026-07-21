# Workflow B — Chained parameter sweep (YAML)

Use this when you want to sweep **one** parameter (e.g. rocking frequency `omega_b`)
across several values using a YAML config file.
Instead of starting each run cold (80-cycle transient), each segment (SLURM job)
restores the checkpoint from the previous one — the flow is already developed,
so only ~10 transition cycles are needed.

!!! tip
    To sweep multiple parameters simultaneously, or to write configs in JSON,
    use [Workflow E](json-sweep.md) instead.

```
seg 0  ──── fresh start, 80 mix cycles ────────────────► checkpoint.dump
seg 1  ── restore checkpoint, 10 transition cycles ────► checkpoint.dump
seg 2  ── restore checkpoint, 10 transition cycles ────► checkpoint.dump
...
```

## Step 1 — Write a chain config

Copy and edit `config/chain_config.yaml`:

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
  parameter: omega_b           # any scalar motion param (see supported list below)
  values: [3.14159, 6.28318]   # 0.5 Hz, 1.0 Hz

# ── Base motion (non-swept params, fixed across all segments) ─────────────────
motion:
  n_harmonics: 1
  theta_max:      [5.0, 0.0, 0.0]
  phi_angular:    [0.3, 0.0, 0.0]
  omega_h:        6.28318
  amplitude_h:    [0.02, 0.0, 0.0]
  phi_horizontal: [0.3, 0.0, 0.0]

# ── Output ───────────────────────────────────────────────────────────────────
videos: false          # set true to use BioReactor-video and render MP4s

# ── SLURM ─────────────────────────────────────────────────────────────────────
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

## Step 2 — Submit

```bash
python scripts/chain.py config/chain_config.yaml
```

Prints one line per segment with run IDs and SLURM job IDs.
SLURM `--dependency=afterok` ensures each segment starts only after the
previous one succeeds.

## Step 3 — Results

Each segment writes independently to `runs/<run_id>/results.json`.
No aggregation script — collect manually or extend `postprocess.py`.

## Before submitting long jobs — always smoke-test first

```bash
# Low-fidelity two-segment dry run (takes ~1 min)
python scripts/chain.py config/chain_config_smoke.yaml
```

Check that:

- `runs/<seg1_id>/logstats.dat` starts at `t > 0` (not at `t=0`)
- `runs/<seg1_id>/results.json` has finite kLa values
