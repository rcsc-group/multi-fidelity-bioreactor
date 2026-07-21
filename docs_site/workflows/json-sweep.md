# Workflow E — JSON multi-param sweep

Use this when you want to sweep over **any combination of parameters** — including
geometry, fidelity, timing, and all motion parameters — using a single JSON file.

Like Workflow B, runs are chained via checkpoint restart within groups of
compatible simulations (same fidelity and geometry). Across groups with different
fidelities or geometries, jobs run independently in parallel.

## How lists are interpreted

Write a normal `params.json` file. Any parameter can be swept by giving it a
**list of values** instead of a single value:

| How you write it | What it means |
|-----------------|---------------|
| `"omega_b": 3.14` | Fixed value — same for every simulation |
| `"omega_b": [3.14, 6.28, 9.42]` | Sweep: try all three values |
| `"fidelity": [5, 7]` | Sweep: run once at fidelity 5, once at fidelity 7 |
| `"theta_max": [7.0, 0.0, 0.0]` | **Fixed** 3-element vector — NOT a sweep |
| `"theta_max": [[5.0,0,0], [7.0,0,0]]` | Sweep: two different rocking amplitudes |
| `"geometry": [{"a":0.2,"b":0.09,"n":2}, {"a":0.25,"b":0.071,"n":8}]` | Sweep: two bag shapes |

**Vector parameters** (`theta_max`, `phi_angular`, `amplitude_h`, `phi_horizontal`)
are only treated as a sweep when their elements are themselves lists (nested lists).
A plain list of three numbers is always interpreted as a single vector value.

## Zip vs. cartesian expansion

- **All swept lists the same length N** → **zip**: N simulations, element k of
  each list goes into simulation k together.
  Example: `"omega_b":[1,2,3], "fill_level":[0.4,0.5,0.6]` → 3 sims:
  `(omega_b=1, fill=0.4)`, `(omega_b=2, fill=0.5)`, `(omega_b=3, fill=0.6)`.

- **Lists of different lengths** → **cartesian product**: every combination.
  Example: `"omega_b":[1,2]`, `"fill_level":[0.3,0.5,0.7]` → 6 sims:
  `(1,0.3)`, `(1,0.5)`, `(1,0.7)`, `(2,0.3)`, `(2,0.5)`, `(2,0.7)`.

## Checkpoint grouping

Checkpoint restart is only valid between simulations that share the same
**fidelity and geometry** (the computational grid must be identical).
The sweep runner automatically groups simulations by `(fidelity, geometry)` and
submits each group as a separate chain. Groups with different fidelities or
geometries run in parallel with no checkpoint between them.

## Step 1 — Write a sweep config

Copy and edit `config/sweep_example.json`:

```json
{
  "fidelity":    3,
  "geometry":    {"a": 0.25, "b": 0.071, "n": 2.0},
  "fill_level":  0.5,
  "n_harmonics": 1,
  "theta_max":   [7.0, 0.0, 0.0],
  "phi_angular": [0.0, 0.0, 0.0],
  "omega_h":     0.0,
  "amplitude_h": [0.0, 0.0, 0.0],
  "phi_horizontal": [0.0, 0.0, 0.0],

  "omega_b": [3.14159, 6.28318],

  "_sweep": {
    "n_mix_cycles":        3,
    "n_transition_cycles": 3,
    "t_buffer":            5.0,
    "walltime":            "00:10:00",
    "submit":              false
  }
}
```

The `"_sweep"` key holds sweep-control options that are **not** simulation parameters:

| Option | Meaning |
|--------|---------|
| `n_mix_cycles` | Rocking cycles before O₂ injection for the **first** segment of each group (fresh start) |
| `n_transition_cycles` | Rocking cycles before O₂ re-injection for **restart** segments (flow is already developed) |
| `t_buffer` | Length of the kLa measurement window in non-dim time (see the [Glossary](../glossary.md) for sizing guidance) |
| `walltime` | SLURM time limit per segment (`HH:MM:SS`) |
| `cpus` | CPUs per job (OpenMP threads). Default: 4. Use 16 for fidelity ≥ 7 |
| `mem` | Memory per job (e.g. `"16G"`). Default: `"12G"` |
| `submit` | `true` → submit via `sbatch`; `false` → write `params.json` files only (dry run) |

!!! note
    If you include `n_mix_cycles` in the **body** of the JSON (not inside `_sweep`)
    as a list, each simulation uses its own value from the sweep — the `_sweep.n_mix_cycles`
    and `_sweep.n_transition_cycles` chain defaults are then ignored for that sweep.

## Step 2 — Dry-run first

Set `"submit": false` and check that params files are written correctly:

```bash
python scripts/sweep.py config/sweep_example.json
```

The script prints one line per segment. Check that `t_end > n_mix_cycles × T_period`
(otherwise kLa will be NaN) and that restart segments have `t_checkpoint > 0`.

## Step 3 — Submit

Set `"submit": true` and run:

```bash
uv run python scripts/sweep.py config/my_sweep.json
```

Output:
```
Group 0 (fidelity=3, a=0.25, b=0.071, n=2.0) — 2 segment(s)
  [seg 0] run=abc12345  omega_b=3.142  n_mix=3  t_end≈6.8  → next:def67890
  [seg 1] run=def67890  omega_b=6.283  n_mix=3  t_end≈6.8  → last
  → submitted seg-0 as job 1234567 (chain self-submits from here)
```

Only seg-0 of each chain is submitted to SLURM upfront. Each segment automatically
submits its successor at the end of its SLURM script. This means at most one job per
chain is ever in the queue at a time — chains progress at their own pace without
waiting on each other's SLURM state.

## Step 4 — Verify results

After all jobs complete:

- Restart segment: `runs/<seg1_id>/logstats.dat` must start at `t > 0` (not at `t=0`)
- All segments: `runs/<id>/results.json` must contain finite (non-NaN) kLa values
