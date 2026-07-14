# Rocking Bioreactor 2D — Simulation Suite

Two-phase CFD solver for a rocking bioreactor, implemented in [Basilisk](http://basilisk.fr/).
Developed at the Harris Lab (Brown University) in collaboration with the Cimpeanu group (Warwick).
Publication: [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538) | preprint: [arXiv: 2504.05421](https://arxiv.org/abs/2504.05421)

---

## Table of contents

0. [Glossary](#0-glossary)
1. [Setup](#1-setup)
2. [Scripts at a glance](#2-scripts-at-a-glance)
3. [Workflow A — Single run](#3-workflow-a--single-run)
4. [Workflow B — Chained parameter sweep (YAML)](#4-workflow-b--chained-parameter-sweep-yaml)
5. [Workflow C — Batch sampling](#5-workflow-c--batch-sampling)
6. [Workflow D — Multi-fidelity Bayesian optimisation](#6-workflow-d--multi-fidelity-bayesian-optimisation)
7. [Workflow E — JSON multi-param sweep](#7-workflow-e--json-multi-param-sweep)
8. [Video generation](#8-video-generation)
9. [params.json reference](#9-paramsjson-reference)
10. [Output files reference](#10-output-files-reference)
11. [Fidelity guide](#11-fidelity-guide)
12. [Project structure](#12-project-structure)
13. [Test suite](#13-test-suite)

---

## 0. Glossary

Plain-language definitions of every technical term used in this document.
If you are new to CFD or HPC, read this section first.

**kLa (mass-transfer coefficient)**
The rate at which oxygen dissolves from the air into the liquid in the bag,
measured per unit time. Higher kLa = better mixing and faster oxygenation.
`kLa_25` is the value when the liquid has reached 25 % of full oxygen saturation
and is the standard industrial metric. Units: h⁻¹ (converted from the solver's
internal non-dimensional rate via `3600 / T_bio`, so it is directly comparable
to literature values such as Kim et al. 2024).

**Rocking period (T)**
The time it takes the bag to complete one full back-and-forth rock.
Related to rocking frequency: T = 2π / ω_b. Shorter period = faster rocking.

**Angular frequency (ω_b, omega\_b)**
How fast the bag rocks, measured in radians per second (rad/s).
1 Hz = 2π ≈ 6.28 rad/s. A typical bioreactor runs at 0.3–2 Hz.

**Non-dimensional time**
Simulation time is scaled by the characteristic sloshing timescale T_bio = L / U_bio,
where L is the bag length and U_bio is the typical sloshing velocity.
This makes results independent of the exact bag size or fluid speed, so
one non-dim time unit means roughly "one sloshing timescale has passed."
Physical time in seconds = t × T_bio.

**Fidelity / grid level**
How fine the computational mesh (grid) is. Grid size = 2^fidelity × 2^fidelity cells.
Higher fidelity = more cells = more accurate results, but slower and more memory-intensive.
Fidelity 3 (8×8) takes seconds; fidelity 7 (128×128) takes hours.
See [§11 Fidelity guide](#11-fidelity-guide) for a full table.

**Checkpoint / checkpoint restart**
A snapshot of the full simulation state (fluid velocity, volume fraction, dissolved oxygen,
etc.) saved to a binary file called `checkpoint.dump`.
A later simulation can *restore* this snapshot and continue from that exact point,
skipping the warm-up phase. This is used in chained sweeps to save 70–90 % of compute time.

**Segment**
One SLURM job in a chained sweep. Each segment is a complete, self-contained simulation
that starts either from scratch (segment 0) or from the checkpoint left by the previous segment.

**Chain**
A sequence of segments linked by checkpoint restart. Each segment automatically submits
its successor when it completes successfully (self-submitting chain). No external scheduler
bookkeeping is required.

**n_mix_cycles**
The number of rocking cycles to run *before* injecting oxygen. Used to let the flow
reach a steady state before taking measurements. Typical value: 80 cycles for a fresh start,
10 for a checkpoint restart (flow is already developed).

**t_buffer**
The duration (in non-dim time) of the kLa measurement window after oxygen injection.
Larger t_buffer gives a longer average and more stable kLa estimate, at the cost of more
compute time. Rule of thumb: t_buffer > ln(2) / kLa_nd_expected, where kLa_nd is the
solver's internal *non-dimensional* rate (not the h⁻¹ value reported in `results.json` —
t_buffer is itself a non-dim duration, so it must be sized against the non-dim rate).
At fidelity 7, kLa_nd ≈ 0.03–0.1, so t_buffer = 30–50 is sufficient. At fidelity 5,
kLa_nd ≈ 0.1–0.5, so t_buffer = 10–30 suffices.

**Superellipse (geometry.n)**
The mathematical shape of the bag cross-section. n=2 is a standard ellipse;
n=8 or more looks like a rounded rectangle. The shape is controlled by
`geometry.a` (half-width), `geometry.b` (half-height), and `geometry.n` (roundness).

**Surrogate model**
A fast mathematical approximation of kLa built by fitting to existing simulation data.
Once trained, it can predict kLa for any parameter combination in milliseconds
(vs. hours for a real simulation). Used in Bayesian optimisation to guide where to sample next.

**Bayesian optimisation (BO)**
An algorithm that iteratively chooses which parameter combination to simulate next,
based on a balance between exploring unknown regions (*exploration*) and refining
around the best results found so far (*exploitation*). Uses the surrogate to evaluate
candidates cheaply. The expected improvement (*EI*) acquisition function quantifies
how promising each candidate is.

**DoE (Design of Experiments)**
An initial set of simulation runs that covers the parameter space broadly before
optimisation begins. Used to build the first surrogate model.
Common strategy: Latin Hypercube Sampling (LHS), which spreads points evenly.

**SLURM**
The job scheduler on OSCAR (Brown University HPC). You submit a job with `sbatch`;
SLURM queues it and runs it on a compute node when resources are available.
Jobs communicate via environment variables and output files — never interact
with SLURM interactively from a login node.

**HPC / OSCAR**
High-Performance Computing cluster at Brown University. Always verify you are on a
*compute node* (not a login node) before running expensive simulations.
Login nodes are shared; compute nodes are allocated exclusively per job.

---

## 1. Setup

### Build the simulation binary

```bash
cd dev/rocking-bioreactor-2d
make build           # standard kLa-only binary  → build/BioReactor
make build-video     # + frame dumps for videos  → build/BioReactor-video
make build-health    # + Poisson diagnostics      → build/BioReactor-health
```

The Makefile calls `/oscar/data/dharri15/eaguerov/basilisk/src/qcc` (built under
the persistent data allocation, not scratch — `/oscar/scratch` purges files
unmodified for ~30 days, and Basilisk is only ever read from, never modified,
so it silently looked stale and got swept there in July 2026. The OSCAR spack
module has a broken header path — never use `module load basilisk`).

### Install Python dependencies

```bash
pip install uv      # one-time, if not already installed
uv sync
```

Run all scripts via `uv run python scripts/foo.py ...` or `uv run pytest`.
`uv` manages the virtual environment automatically — no manual activation needed.

---

## 2. Scripts at a glance

| Script | How to run | What it does |
|--------|-----------|--------------|
| `simulate.py` | `uv run python simulate.py params.json --slurm --wait` | Submit or run one simulation; wait for `results.json` |
| `launch.py` | `uv run python launch.py params.json` | Set up a run directory and write a SLURM script, but do **not** submit |
| `chain.py` | `uv run python chain.py config.yaml` | Sweep **one** parameter across values via a YAML config, using checkpoint restart between segments |
| `sweep.py` | `uv run python sweep.py config.json` | Sweep **any combination** of parameters from a JSON file (zip or cartesian); groups runs by checkpoint compatibility; chains self-submit |
| `sample.py` | `uv run python sample.py config.yaml` | Space-filling batch of independent runs (Latin Hypercube or random) |
| `loop.py` | `uv run python loop.py config.yaml` | Full Bayesian optimisation loop (DoE → train surrogate → suggest → repeat) |
| `postprocess.py` | `uv run python postprocess.py runs/my_run/` | Extract kLa, mixing time, and vorticity from run output files → `results.json` |
| `collect_results.py` | `uv run python collect_results.py --sweep config/my_sweep.json` | Aggregate all `results.json` files into a single CSV |
| `plot_heatmaps.py` | `uv run python plot_heatmaps.py` | Generate KPI heatmap figures from all completed sweep results → `experiments/figures/` |
| `render_videos.py` | `uv run python render_videos.py runs/my_run/` | Convert raw frame dumps to MP4 videos (called automatically by SLURM jobs) |
| `suggest.py` | `uv run python suggest.py exp_dir param_space.yaml` | Print the next highest-EI parameter point to stdout |
| `train_surrogate.py` | `uv run python train_surrogate.py exp_dir model.pkl` | Train the multi-fidelity surrogate from existing run data |

---

## 3. Workflow A — Single run

Use this for one-off simulations, debugging, or manual parameter exploration.

### Step 1 — Write a params.json

Create `runs/my_run/params.json` (see [§9](#9-paramsjson-reference) for all fields):

```json
{
  "run_id":         "my_run",
  "fidelity":       7,
  "omega_b":        3.93,
  "n_harmonics":    1,
  "theta_max":      [7.0, 0.0, 0.0],
  "phi_angular":    [0.0, 0.0, 0.0],
  "omega_h":        0.0,
  "amplitude_h":    [0.0, 0.0, 0.0],
  "phi_horizontal": [0.0, 0.0, 0.0],
  "geometry":       {"a": 0.25, "b": 0.071, "n": 8.0},
  "fill_level":     0.5,
  "n_mix_cycles":   80,
  "t_end":          250.0
}
```

### Step 2 — Run

**Local (blocks terminal; use fidelity 3–4 only):**
```bash
make run PARAMS=runs/my_run/params.json
```

**SLURM (non-blocking; production):**
```bash
make submit PARAMS=runs/my_run/params.json
# or, with Python — also waits for results:
python scripts/simulate.py runs/my_run/params.json --slurm --wait --walltime 04:00:00
```

### Step 3 — Postprocess

```bash
uv run python scripts/postprocess.py runs/my_run/
```

Writes `runs/my_run/results.json` with eighteen KPIs (kLa/tau are dimensional;
NaN when the underlying data file is absent or too short, e.g. tau_* are NaN
whenever oxygen transport was disabled for the run, and vice versa):

```json
{
  "kLa_10": 47.87,  "kLa_25": 45.14,  "kLa_50": 29.01,
  "kLa_inst_10": 50.79, "kLa_inst_25": 44.73, "kLa_inst_50": 27.62,
  "dtmix_0.50": 13.66,  "dtmix_0.75": 27.63,  "dtmix_0.95": 86.53,
  "vor_mean": 1.028,    "vel_rms_qss": 0.250, "kla_fit_rmse_25": 0.00046,
  "tau_95_qss": 0.0023, "tau_98_qss": 0.0108, "tau_100_qss": 0.0986,
  "tau_95_max": 0.0028, "tau_98_max": 0.0138, "tau_100_max": 0.1458
}
```

| Key | Description | Unit |
|-----|-------------|------|
| `kLa_10/25/50` | O₂ transfer rate at 10/25/50 % saturation (5-point log-linear fit) | h⁻¹ |
| `kLa_inst_10/25/50` | Same, estimated instantaneously via dC*/dt | h⁻¹ |
| `dtmix_0.50/0.75/0.95` | Time for tracer to reach 50/75/95 % homogeneity | seconds |
| `vor_mean` | Period-averaged mean absolute vorticity (steady streaming strength) | 1/s |
| `vel_rms_qss` | RMS velocity over the quasi-steady-state window | non-dim |
| `kla_fit_rmse_25` | Fit quality (RMSE) of the kLa_25 log-linear regression | non-dim |
| `tau_95/98/100_qss` | Median shear-stress percentile over the QSS window | Pa |
| `tau_95/98/100_max` | Global-max shear-stress percentile over the whole run | Pa |

kLa is converted from the solver's internal non-dimensional rate via `× 3600 / T_bio`
(h⁻¹). tau is converted via `× ρ_w U_bio²` (Pa), where `U_bio = geometry.a / T_bio`
(BioReactor.c sets `rho1=1`, `mu1=1/Re_w`, so `τ_nd = τ_dim / (ρ_w U_bio²)`, not
`τ_dim × T_bio / μ_w`). Both conversions make results directly comparable to
Kim et al. 2024 and other dimensional literature values.

`kLa_25` is the standard industrial metric. `vor_mean` is the hydrodynamic
root cause: stronger steady streaming → faster mixing and higher kLa.

### Alternatively — set up the directory without submitting

`launch.py` writes the run directory and a SLURM script but does **not** submit:

```bash
python scripts/launch.py path/to/params.json [runs_root]
# returns {"run_id": "...", "run_dir": "...", "slurm_script": "..."}
```

---

## 4. Workflow B — Chained parameter sweep (YAML)

Use this when you want to sweep **one** parameter (e.g. rocking frequency `omega_b`)
across several values using a YAML config file.
Instead of starting each run cold (80-cycle transient), each segment (SLURM job)
restores the checkpoint from the previous one — the flow is already developed,
so only ~10 transition cycles are needed.

> **Tip:** To sweep multiple parameters simultaneously, or to write configs in JSON,
> use [Workflow E](#7-workflow-e--json-multi-param-sweep) instead.

```
seg 0  ──── fresh start, 80 mix cycles ────────────────► checkpoint.dump
seg 1  ── restore checkpoint, 10 transition cycles ────► checkpoint.dump
seg 2  ── restore checkpoint, 10 transition cycles ────► checkpoint.dump
...
```

### Step 1 — Write a chain config

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

### Step 2 — Submit

```bash
python scripts/chain.py config/chain_config.yaml
```

Prints one line per segment with run IDs and SLURM job IDs.
SLURM `--dependency=afterok` ensures each segment starts only after the
previous one succeeds.

### Step 3 — Results

Each segment writes independently to `runs/<run_id>/results.json`.
No aggregation script — collect manually or extend `postprocess.py`.

### Before submitting long jobs — always smoke-test first

```bash
# Low-fidelity two-segment dry run (takes ~1 min)
python scripts/chain.py config/chain_config_smoke.yaml
```

Check that:
- `runs/<seg1_id>/logstats.dat` starts at `t > 0` (not at `t=0`)
- `runs/<seg1_id>/results.json` has finite kLa values

---

## 5. Workflow C — Batch sampling

Use this to generate a space-filling design of experiments (DoE) across the parameter
space without running the BO loop.

### Step 1 — Write a sample config

Copy and edit `config/sample_config.yaml`:

```yaml
experiment_dir: experiments/lhs_run_001

sampling:
  strategy: latin    # latin | random | grid | sobol
  n_samples: 20
  seed: 0

fidelity: 5          # 32×32 grid, ~5–15 min each

# t_end is computed per-run as t_mix(params) + t_buffer (adapts to each frequency)
t_buffer: 150

submit: true
walltime: "00:30:00"
```

### Step 2 — Submit

```bash
python scripts/sample.py config/sample_config.yaml
```

All runs are submitted in parallel. Results are stored in an
[f3dasm](https://github.com/bessagroup/f3dasm) `ExperimentData` store at
`experiment_dir/experiment_data/`.

### Step 3 — Inspect results

```python
from f3dasm import ExperimentData
data = ExperimentData.from_file("experiments/lhs_run_001/experiment_data")
print(data)
```

---

## 6. Workflow D — Multi-fidelity Bayesian optimisation

Maximises kLa over the parameter space using a multi-fidelity surrogate
(kernel-ridge regression + linear regression + Gaussian process).
Low-fidelity (fidelity 5) runs build the correlation structure cheaply;
high-fidelity (fidelity 7) runs target the optimum.

```
Phase 1 — Initial DoE (parallel)
  ├── n_lf_init LHS samples at lf_fidelity
  └── n_hf_init LHS samples at hf_fidelity

Phase 2 — BO loop (n_iter sequential iterations)
  ├── train MF surrogate on all existing data
  ├── maximise Expected Improvement over n_candidates random candidates
  └── submit one HF run at the highest-EI point
```

### Step 1 — Write a BO config

Copy and edit `config/bo_config.yaml`:

```yaml
experiment_dir: experiments/bo_run_001

lf_fidelity: 5   # 32×32, ~5–15 min
hf_fidelity: 7   # 128×128, ~1–2 h

n_lf_init: 8     # LF samples for initial DoE (rule of thumb: ≈10× free params)
n_hf_init: 3     # HF samples for initial DoE (seeds the transfer-learning correlation)
n_iter: 20       # BO iterations after DoE

kla_key: kLa_25  # objective: kLa_10 | kLa_25 | kLa_50
n_candidates: 2000

walltime: "02:00:00"
job_timeout: 7200    # seconds to wait for results.json per job
t_buffer: 150
```

The parameter bounds are defined in `config/param_space.yaml`.

### Step 2 — Run the loop

```bash
python scripts/loop.py config/bo_config.yaml
```

The loop is **resumable**: it checks how many runs already exist in
`experiment_dir` and skips the DoE phase if complete. Kill and restart freely.

### Surrogate tools (standalone)

**Train the surrogate manually:**
```bash
python scripts/train_surrogate.py experiments/bo_run_001 surrogate/model.pkl \
    --lf-fidelity 5 --hf-fidelity 7 --kla-key kLa_25
```
Writes a pickled `model.pkl` with a `.predict(X) -> (mean, var)` interface.

**Query the next suggested point:**
```bash
python scripts/suggest.py experiments/bo_run_001 config/param_space.yaml \
    --model-path surrogate/model.pkl --kla-key kLa_25 --n-candidates 2000
```
Prints a JSON params dict to stdout (the highest-EI candidate).

---

## 7. Workflow E — JSON multi-param sweep

Use this when you want to sweep over **any combination of parameters** — including
geometry, fidelity, timing, and all motion parameters — using a single JSON file.

Like Workflow B, runs are chained via checkpoint restart within groups of
compatible simulations (same fidelity and geometry). Across groups with different
fidelities or geometries, jobs run independently in parallel.

### How lists are interpreted

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

### Zip vs. cartesian expansion

- **All swept lists the same length N** → **zip**: N simulations, element k of
  each list goes into simulation k together.
  Example: `"omega_b":[1,2,3], "fill_level":[0.4,0.5,0.6]` → 3 sims:
  `(omega_b=1, fill=0.4)`, `(omega_b=2, fill=0.5)`, `(omega_b=3, fill=0.6)`.

- **Lists of different lengths** → **cartesian product**: every combination.
  Example: `"omega_b":[1,2]`, `"fill_level":[0.3,0.5,0.7]` → 6 sims:
  `(1,0.3)`, `(1,0.5)`, `(1,0.7)`, `(2,0.3)`, `(2,0.5)`, `(2,0.7)`.

### Checkpoint grouping

Checkpoint restart is only valid between simulations that share the same
**fidelity and geometry** (the computational grid must be identical).
The sweep runner automatically groups simulations by `(fidelity, geometry)` and
submits each group as a separate chain. Groups with different fidelities or
geometries run in parallel with no checkpoint between them.

### Step 1 — Write a sweep config

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
| `t_buffer` | Length of the kLa measurement window in non-dim time (see Glossary for sizing guidance) |
| `walltime` | SLURM time limit per segment (`HH:MM:SS`) |
| `cpus` | CPUs per job (OpenMP threads). Default: 4. Use 16 for fidelity ≥ 7 |
| `mem` | Memory per job (e.g. `"16G"`). Default: `"12G"` |
| `submit` | `true` → submit via `sbatch`; `false` → write `params.json` files only (dry run) |

> **Note:** If you include `n_mix_cycles` in the **body** of the JSON (not inside `_sweep`)
> as a list, each simulation uses its own value from the sweep — the `_sweep.n_mix_cycles`
> and `_sweep.n_transition_cycles` chain defaults are then ignored for that sweep.

### Step 2 — Dry-run first

Set `"submit": false` and check that params files are written correctly:

```bash
python scripts/sweep.py config/sweep_example.json
```

The script prints one line per segment. Check that `t_end > n_mix_cycles × T_period`
(otherwise kLa will be NaN) and that restart segments have `t_checkpoint > 0`.

### Step 3 — Submit

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

### Step 4 — Verify results

After all jobs complete:
- Restart segment: `runs/<seg1_id>/logstats.dat` must start at `t > 0` (not at `t=0`)
- All segments: `runs/<id>/results.json` must contain finite (non-NaN) kLa values

---

## 8. Video generation

Videos are generated automatically during every simulation run — no separate step
is required. All SLURM jobs use `BioReactor-video`, which renders frames inline
and encodes them to MP4 at the end of each segment.

### Output files produced per run

| File | Content |
|------|---------|
| `vorticity3.mp4` | Vorticity field (body frame) |
| `volume_fraction3.mp4` | VOF interface (body frame) |
| `oxygen3.mp4` | Dissolved oxygen concentration |
| `tracer.mp4` | Tracer mixing (top-half injection) |

All MP4s land in `runs/<run_id>/` alongside the data files.

### Build requirements

```bash
make build-video    # compiles BioReactor-video (only needed if binary is stale)
```

`ffmpeg` and Basilisk's `ppm2mp4` helper are loaded automatically by the SLURM
template (`module load ffmpeg`). No additional setup is required.

---

## 9. params.json reference

All fields for `runs/<run_id>/params.json`.

### Motion

| Field | Type | Units | Default | Description |
|-------|------|-------|---------|-------------|
| `omega_b` | float | rad/s | required | Fundamental rocking angular frequency (how fast the bag rocks back and forth). 1 Hz = 2π ≈ 6.28 rad/s |
| `n_harmonics` | int | — | 1 | Number of active harmonics (1–3). A harmonic is a frequency component; 1 means pure sinusoidal rocking at `omega_b`. Vectors always padded to length 3 |
| `theta_max` | float[3] | degrees | [7,0,0] | Maximum rocking angle per harmonic. Index 0 is the fundamental (dominant) harmonic. Typical range: 2–15 degrees |
| `phi_angular` | float[3] | rad | [0,0,0] | Phase offset of rocking per harmonic (delays or advances the timing). Index 0 is **always forced to 0** — it is the global time reference |
| `omega_h` | float | rad/s | 0.0 | Horizontal translation frequency. Set to 0.0 to disable horizontal motion |
| `amplitude_h` | float[3] | m | [0,0,0] | Horizontal translation amplitude per harmonic (how far the bag slides sideways) |
| `phi_horizontal` | float[3] | rad | [0,0,0] | Phase offset of horizontal translation per harmonic |

### Geometry

| Field | Type | Units | Default | Description |
|-------|------|-------|---------|-------------|
| `geometry.a` | float | m | 0.25 | Bag half-width (horizontal semi-axis; half the total bag width) |
| `geometry.b` | float | m | 0.071 | Bag half-height (vertical semi-axis; half the total bag height) |
| `geometry.n` | float | — | 8.0 | Superellipse exponent controlling bag shape: n=2 gives an ellipse, n≥8 gives a rounded rectangle |
| `fill_level` | float | fraction | 0.5 | Fraction of bag volume filled with liquid (0 = empty, 1 = full). Typical range: 0.3–0.7 |

### Simulation control

| Field | Type | Units | Default | Description |
|-------|------|-------|---------|-------------|
| `run_id` | string | — | required | Unique label for this run; output files go to `runs/{run_id}/` |
| `fidelity` | int | — | required | Basilisk grid level; the computational mesh is 2^fidelity × 2^fidelity cells. Higher = more accurate and slower. See [§11](#11-fidelity-guide) |
| `n_mix_cycles` | int | — | 80 | Number of complete rocking cycles to run before injecting oxygen. Used to let the flow field reach a steady state before kLa measurement begins |
| `t_end` | float | non-dim | 250.0 | When to stop the simulation (in non-dimensional time). Computed automatically by `simulate.py`, `chain.py`, and `sweep.py`; only set manually for custom runs |

### Checkpoint restart (set automatically by chain.py and sweep.py)

Do not set these manually — they are populated by the sweep scripts.

| Field | Type | Description |
|-------|------|-------------|
| `t_checkpoint` | float | Absolute non-dim time at which the restored checkpoint was saved (0 for fresh runs) |
| `omega_b_prev` | float | Rocking frequency of the segment that wrote the checkpoint; used to smoothly ramp to the new frequency |
| `theta_max_prev` | float[3] | Rocking amplitude of the previous segment |
| `phi_angular_prev` | float[3] | Rocking phase of the previous segment |
| `amplitude_h_prev` | float[3] | Horizontal translation amplitude of the previous segment |
| `phi_horizontal_prev` | float[3] | Horizontal translation phase of the previous segment |
| `omega_h_prev` | float | Horizontal translation frequency of the previous segment |

### Notes

**`t_end` is non-dimensional.** One non-dim time unit ≈ T_bio = L_bio / U_bio seconds,
where U_bio is the characteristic sloshing velocity (function of geometry, fill, and omega_b).
`simulate.py` computes `t_end = t_mix + t_buffer` automatically from `n_mix_cycles` and the
config's `t_buffer`; you only need to set `t_end` manually for custom runs.

**`phi_angular[0]` is always 0.** It is the global time-origin reference for the
rocking phase; it is physically redundant and is overridden at parse time.

**`omega_b` and `omega_h` are independent in the model.** On a physical platform
they are driven by the same motor (ω_h = ω_b). For pure rocking with no horizontal
translation, set `omega_h: 0.0` and all `amplitude_h` to zero.

---

## 10. Output files reference

All files land in `runs/{run_id}/`.

| File | Written by | Columns / content |
|------|-----------|-------------------|
| `params.json` | Python scripts | Input parameters (copied in) |
| `logstats.dat` | BioReactor | `i t dt #Cells wall_time cpu_time` — one row per 0.1 non-dim time |
| `normf.dat` | BioReactor | `i t Omega_liq_avg Omega_liq_rms ... ux ... uy ...` — vorticity and velocity norms in liquid phase |
| `vol_frac_interf.dat` | BioReactor | `i t f_liq_sum f_liq_interf posY_max posY_min` — liquid volume, interface length, interface y-extent |
| `tr_oxy.dat` | BioReactor | `i t oxy_liq_sum oxy_liq_sum2 [c_liq...] c2_liq_sum c2_liq_sum2 ...` — dissolved O₂ and tracer integrals; written from `t_mix` onward |
| `results.json` | postprocess.py | Ten KPIs: `kLa_10/25/50`, `kLa_inst_10/25/50`, `dtmix_0.50/0.75/0.95`, `vor_mean` |
| `checkpoint.dump` | BioReactor | Binary Basilisk dump at the end of each run (for chain restart) |
| `vorticity3.mp4` | BioReactor-video | Vorticity field animation |
| `volume_fraction3.mp4` | BioReactor-video | VOF interface animation (body frame) |
| `oxygen3.mp4` | BioReactor-video | Dissolved oxygen concentration |
| `tracer.mp4` | BioReactor-video | Tracer mixing animation |

All `.dat` files have a one-line header beginning with `i`.
Time `t` is non-dimensional; `t_physical = t × T_bio`.

`kLa` is `NaN` if `tr_oxy.dat` contains no data — this happens when `t_end < t_mix`
(run finished before oxygen injection started).

---

## 11. Fidelity guide

Runtimes are per segment (one SLURM job) with the sweep t_end of ~35–98 nondim units.

| Fidelity | Grid | 4 CPUs | 16 CPUs | Use for |
|----------|------|--------|---------|---------|
| 3 | 8×8 | seconds | — | Smoke tests, CI |
| 4 | 16×16 | ~2 min | — | Quick debugging |
| 5 | 32×32 | ~2 min | ~1 min | BO DoE, low-fidelity surrogate |
| 6 | 64×64 | ~15 min | ~5 min | Mid-fidelity check |
| 7 | 128×128 | ~2 h | ~30 min | Standard production sweeps |
| 8 | 256×256 | ~16 h (approx.) | ~4 h (approx.) | Grid-convergence check vs. fidelity 7 |
| 9 | 512×512 | ~days | ~10–11 h/condition (measured) | High-fidelity reference |
| 10 | 1024×1024 | — | ~2.7–6.6 days/condition (measured) | Publication (Kim et al. 2024/2025) |

For sweeps: use **fidelity 5** (fast, qualitative) or **fidelity 7 at 16 CPUs** (production).
Specify `"cpus": 16` in `_sweep` when using fidelity ≥ 7.

For the optimization suite, `lf_fidelity: 5` and `hf_fidelity: 7` are the
recommended pair. Always smoke-test new workflows at fidelity 3 before
submitting fidelity-7 jobs.

---

## 12. Project structure

```
rocking-bioreactor-2d/
├── src/
│   ├── BioReactor.c          # Basilisk solver (compile-time flags: VIDEOS, DIAGNOSTICS, AMR)
│   ├── params_read.h         # JSON → BioreactorParams struct (jsmn-based)
│   ├── henry_oxy2.h          # Henry's law oxygen transport
│   ├── utils2.h, view3.h     # Basilisk visualization helpers
│   └── jsmn.h                # Minimal JSON parser (MIT, Serge Zaitsev)
│
├── scripts/
│   ├── simulate.py           # Core API: run_local() / submit_slurm() / wait_for_result()
│   ├── postprocess.py        # 18 KPIs from run output files → results.json
│   ├── collect_results.py    # Aggregate all results.json across runs → CSV
│   ├── plot_heatmaps.py      # KPI heatmaps from sweep results → experiments/figures/
│   ├── render_videos.py      # Called automatically by SLURM: frame dumps → mp4
│   ├── launch.py             # set up run directory + SLURM script (no submit)
│   ├── chain.py              # Workflow B: chained sweep (one param, YAML config)
│   ├── sweep.py              # Workflow E: JSON multi-param sweep (zip / cartesian; self-submitting chains)
│   ├── sample.py             # Workflow C: LHS / random / grid / Sobol batch sampling
│   ├── loop.py               # Workflow D: multi-fidelity BO loop
│   ├── suggest.py            # EI acquisition: suggest next HF point
│   └── train_surrogate.py    # Train KRR-LR-GPR multi-fidelity surrogate
│
├── config/
│   ├── slurm_template.sh         # SBATCH script for all simulation runs (videos included)
│   ├── slurm_video_template.sh   # SBATCH script for standalone video reruns
│   ├── param_space.yaml          # Parameter bounds for optimization suite
│   ├── bo_config.yaml            # Workflow D config (multi-fidelity BO)
│   ├── sample_config.yaml        # Workflow C config (batch sampling)
│   ├── chain_config.yaml         # Workflow B config (generic sweep template)
│   ├── chain_config_ellipse.yaml # Workflow B config (ellipse bag, 4 frequencies)
│   ├── chain_config_smoke.yaml   # Workflow B config (fidelity-3 smoke test)
│   ├── sweep_example.json        # Workflow E example / smoke test (2-segment, fidelity 3)
│   ├── sweep_fb_theta.json       # theta_max × omega_b sweep (fidelity 5)
│   ├── sweep_fb_theta_l7.json    # theta_max × omega_b sweep (fidelity 7, 16 CPUs)
│   ├── sweep_fb_fill.json        # fill_level × omega_b sweep (fidelity 5)
│   └── sweep_fb_fill_l7.json     # fill_level × omega_b sweep (fidelity 7, 16 CPUs)
│
├── tests/
│   ├── conftest.py               # Shared fixtures (run_bioreactor, CANONICAL_PARAMS)
│   ├── test_chain.py             # chain.py unit tests (build_chain)
│   ├── test_sweep.py             # sweep.py unit tests (detect, expand, group, build, self-submit)
│   ├── test_sample.py            # sample.py unit tests (row_to_params, run_sampling)
│   ├── test_param_schema.py      # params_read schema tests
│   ├── test_launch.py            # launch.py unit tests
│   ├── test_suggest.py           # suggest.py unit tests
│   ├── test_train_surrogate.py   # surrogate training unit tests
│   ├── test_postprocess.py       # postprocessing unit tests
│   ├── test_simulate.py          # simulate.py unit tests (mocked sbatch)
│   ├── integration/
│   │   └── test_bioreactor_runs.py  # end-to-end output file checks
│   └── verification/
│       ├── test_grid_convergence.py  # velocity RMS convergence L5 vs L6
│       ├── test_mass_conservation.py
│       ├── test_oxygen_monotonicity.py
│       └── ...                       # other physics verification tests
│
├── experiments/
│   └── figures/              # Heatmap PDFs generated by plot_heatmaps.py (committed)
├── surrogate/                # Pickled surrogate models
├── runs/                     # Per-run I/O directories (gitignored)
├── build/                    # Compiled binaries (gitignored)
├── logs/                     # SLURM stdout/stderr logs (gitignored)
├── Makefile                  # make build / build-video / run / submit / clean
└── pyproject.toml            # Python dependencies (managed by uv)
```

---

## 13. Test suite

`pyproject.toml` sets `addopts = -m 'not medium and not hpc'`, so the default
invocation runs only the fast, hermetic unit tests. `medium` and `hpc` tests are
opt-in — pass an explicit `-m` to include them (this overrides addopts, it does
not combine with it).

```bash
# Fast unit tests only (no binary, no SLURM — ~20 s). This is the default.
uv run python -m pytest tests/

# Include numerical verification tests (runs the real binary at fidelity 3, ~2 min each)
uv run python -m pytest tests/ -m medium

# Full suite including SLURM integration tests (submits real jobs, needs cluster allocation)
uv run python -m pytest tests/ -m "medium or hpc"
```

---

## References

- Kim M., Harris D.M., Cimpeanu R. (2025). *Modelling of oxygen transfer in rocking bioreactors.* Int. J. Multiphase Flow. [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538)
- Basilisk CFD: http://basilisk.fr
