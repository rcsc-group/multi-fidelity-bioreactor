# Rocking Bioreactor 2D — Simulation Suite

Two-phase CFD solver for a rocking bioreactor, implemented in [Basilisk](http://basilisk.fr/).
Developed at the Harris Lab (Brown University) in collaboration with the Cimpeanu group (Warwick).
Publication: [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538) | preprint: [arXiv: 2504.05421](https://arxiv.org/abs/2504.05421)

---

## Table of contents

1. [Setup](#1-setup)
2. [Workflow A — Single run](#2-workflow-a--single-run)
3. [Workflow B — Chained parameter sweep](#3-workflow-b--chained-parameter-sweep)
4. [Workflow C — Batch sampling](#4-workflow-c--batch-sampling)
5. [Workflow D — Multi-fidelity Bayesian optimisation](#5-workflow-d--multi-fidelity-bayesian-optimisation)
6. [Video generation](#6-video-generation)
7. [params.json reference](#7-paramsjson-reference)
8. [Output files reference](#8-output-files-reference)
9. [Fidelity guide](#9-fidelity-guide)
10. [Project structure](#10-project-structure)
11. [Test suite](#11-test-suite)

---

## 1. Setup

### Build the simulation binary

```bash
cd dev/rocking-bioreactor-2d
make build           # standard kLa-only binary  → build/BioReactor
make build-video     # + frame dumps for videos  → build/BioReactor-video
make build-health    # + Poisson diagnostics      → build/BioReactor-health
```

The Makefile calls `~/scratch/basilisk/src/qcc` (scratch-compiled; the OSCAR
spack module has a broken header path — never use `module load basilisk`).

### Install Python dependencies

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

All scripts must be run from within `.venv`. Throughout this document
`.venv/bin/python` is abbreviated as `python` for clarity.

---

## 2. Workflow A — Single run

Use this for one-off simulations, debugging, or manual parameter exploration.

### Step 1 — Write a params.json

Create `runs/my_run/params.json` (see [§7](#7-paramsjson-reference) for all fields):

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

### Step 3 — Extract kLa

```bash
python scripts/postprocess.py runs/my_run/
```

Writes `runs/my_run/results.json`:
```json
{"kLa_10": 0.0031, "kLa_25": 0.0028, "kLa_50": 0.0024}
```

`kLa_10/25/50` are the mass-transfer coefficients at 10 %, 25 %, and 50 % oxygen saturation.
`kLa_25` is the standard industry metric.

### Alternatively — set up the directory without submitting

`launch.py` writes the run directory and a SLURM script but does **not** submit:

```bash
python scripts/launch.py path/to/params.json [runs_root]
# returns {"run_id": "...", "run_dir": "...", "slurm_script": "..."}
```

---

## 3. Workflow B — Chained parameter sweep

Use this when you want to sweep one parameter (e.g. rocking frequency) across
several values.  Instead of starting each run cold (80 cycle transient), each
segment restores the checkpoint from the previous one — the flow field is
already developed, so only ~10 transition cycles are needed.

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

## 4. Workflow C — Batch sampling

Use this to generate a space-filling design of experiments across the parameter
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

## 5. Workflow D — Multi-fidelity Bayesian optimisation

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

## 6. Video generation

Videos require the `BioReactor-video` binary and `ffmpeg`.

### Single video run

```bash
make build-video
make submit PARAMS=runs/my_run/params.json TEMPLATE=config/slurm_video_template.sh
```

Or use `chain.py` with `videos: true` — it automatically selects the video template
and calls `render_videos.py` after each segment.

### Render videos from an existing run directory

```bash
python scripts/render_videos.py runs/my_run/
```

Reads `frames/frame_XXXXXX.bin` dumps from the video binary and produces:

| File | View |
|------|------|
| `volume_fraction.mp4` | Body frame (bag at rest, fluid moves) |
| `volume_fraction_lab.mp4` | Lab frame (bag rocks, fluid follows) |

The `frames/` directory is deleted after encoding.
`ffmpeg` must be available (`module load ffmpeg` on OSCAR).

---

## 7. params.json reference

All fields for `runs/<run_id>/params.json`.

### Motion

| Field | Type | Units | Default | Description |
|-------|------|-------|---------|-------------|
| `omega_b` | float | rad/s | required | Fundamental rocking angular frequency. 1 Hz = 2π ≈ 6.28 rad/s |
| `n_harmonics` | int | — | 1 | Number of active harmonics (1–3). Vectors always padded to length 3 |
| `theta_max` | float[3] | degrees | [7,0,0] | Rocking amplitude per harmonic. Index 0 is the fundamental |
| `phi_angular` | float[3] | rad | [0,0,0] | Rocking phase per harmonic. Index 0 is **always forced to 0** at read time |
| `omega_h` | float | rad/s | 0.0 | Horizontal translation frequency. Set to 0 to disable |
| `amplitude_h` | float[3] | m | [0,0,0] | Horizontal translation amplitude per harmonic |
| `phi_horizontal` | float[3] | rad | [0,0,0] | Horizontal translation phase per harmonic |

### Geometry

| Field | Type | Units | Default | Description |
|-------|------|-------|---------|-------------|
| `geometry.a` | float | m | 0.25 | Bag half-width (horizontal semi-axis) |
| `geometry.b` | float | m | 0.071 | Bag half-height (vertical semi-axis) |
| `geometry.n` | float | — | 8.0 | Superellipse exponent: 2 = ellipse, n ≥ 8 ≈ rectangle |
| `fill_level` | float | fraction | 0.5 | Fraction of bag volume filled with liquid (0.3–0.7) |

### Simulation control

| Field | Type | Units | Default | Description |
|-------|------|-------|---------|-------------|
| `run_id` | string | — | required | Unique tag; output goes to `runs/{run_id}/` |
| `fidelity` | int | — | required | Basilisk grid level: grid is 2^fidelity × 2^fidelity cells |
| `n_mix_cycles` | int | — | 80 | Rocking cycles before O₂ is injected. Ensures flow is fully developed |
| `t_end` | float | non-dim | 250.0 | Simulation end time. Computed automatically by `simulate.py` and `chain.py`; set manually only for custom runs |

### Checkpoint restart (set by chain.py; do not set manually)

| Field | Type | Description |
|-------|------|-------------|
| `t_checkpoint` | float | Absolute non-dim time of the restored checkpoint (0 for fresh runs) |
| `omega_b_prev` | float | `omega_b` of the segment that wrote the checkpoint (used to rescale fields) |

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

## 8. Output files reference

All files land in `runs/{run_id}/`.

| File | Written by | Columns / content |
|------|-----------|-------------------|
| `params.json` | Python scripts | Input parameters (copied in) |
| `logstats.dat` | BioReactor | `i t dt #Cells wall_time cpu_time` — one row per 0.1 non-dim time |
| `normf.dat` | BioReactor | `i t Omega_liq_avg Omega_liq_rms ... ux ... uy ...` — vorticity and velocity norms in liquid phase |
| `vol_frac_interf.dat` | BioReactor | `i t f_liq_sum f_liq_interf posY_max posY_min` — liquid volume, interface length, interface y-extent |
| `tr_oxy.dat` | BioReactor | `i t oxy_liq_sum oxy_liq_sum2 c_liq_sum ...` — dissolved O₂ and tracer integrals; written from `t_mix` onward |
| `results.json` | postprocess.py | `{"kLa_10": ..., "kLa_25": ..., "kLa_50": ...}` — kLa at 10/25/50 % saturation |
| `checkpoint.dump` | BioReactor | Binary Basilisk dump at the end of each run (for chain restart) |
| `volume_fraction.mp4` | render_videos.py | Body-frame VOF animation (requires BioReactor-video) |
| `volume_fraction_lab.mp4` | render_videos.py | Lab-frame VOF animation (requires BioReactor-video) |

All `.dat` files have a one-line header beginning with `i`.
Time `t` is non-dimensional; `t_physical = t × T_bio`.

`kLa` is `NaN` if `tr_oxy.dat` contains no data — this happens when `t_end < t_mix`
(run finished before oxygen injection started).

---

## 9. Fidelity guide

| fidelity | grid | typical runtime | use for |
|----------|------|-----------------|---------|
| 3 | 8×8 | seconds | import / smoke test only |
| 4 | 16×16 | ~2 min | quick single-parameter scan |
| 5 | 32×32 | ~5–15 min | low-fidelity surrogate data, BO DoE |
| 6 | 64×64 | ~30 min | mid-fidelity check |
| 7 | 128×128 | ~2–4 h | standard production / HF BO runs |
| 9 | 512×512 | ~days | high-fidelity reference / publication |

For the optimization suite, `lf_fidelity: 5` and `hf_fidelity: 7` are the
recommended pair. Always smoke-test new workflows at fidelity 3 before
submitting fidelity-7 jobs.

---

## 10. Project structure

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
│   ├── postprocess.py        # kLa extraction from tr_oxy.dat → results.json
│   ├── render_videos.py      # frame dumps → volume_fraction*.mp4 via ffmpeg
│   ├── launch.py             # set up run directory + SLURM script (no submit)
│   ├── chain.py              # Workflow B: chained sweep via checkpoint restart
│   ├── sample.py             # Workflow C: LHS / random / grid / Sobol batch sampling
│   ├── loop.py               # Workflow D: multi-fidelity BO loop
│   ├── suggest.py            # EI acquisition: suggest next HF point
│   └── train_surrogate.py    # Train KRR-LR-GPR multi-fidelity surrogate
│
├── config/
│   ├── slurm_template.sh         # SBATCH script for standard kLa runs
│   ├── slurm_video_template.sh   # SBATCH script for video runs (loads ffmpeg)
│   ├── param_space.yaml          # Parameter bounds for optimization suite
│   ├── bo_config.yaml            # Workflow D config (multi-fidelity BO)
│   ├── sample_config.yaml        # Workflow C config (batch sampling)
│   ├── chain_config.yaml         # Workflow B config (generic sweep template)
│   ├── chain_config_ellipse.yaml # Workflow B config (ellipse bag, 4 frequencies)
│   └── chain_config_smoke.yaml   # Workflow B config (fidelity-3 smoke test)
│
├── tests/
│   ├── conftest.py               # Shared fixtures (run_bioreactor, CANONICAL_PARAMS)
│   ├── test_chain.py             # chain.py unit tests (build_chain)
│   ├── test_param_schema.py      # params_read schema tests
│   ├── test_launch.py            # launch.py unit tests
│   ├── test_suggest.py           # suggest.py unit tests
│   ├── test_train_surrogate.py   # surrogate training unit tests
│   ├── test_postprocess.py       # kLa extraction unit tests
│   └── test_simulate.py          # simulate.py unit tests (mocked sbatch)
│
├── experiments/              # f3dasm ExperimentData stores (Workflow C & D output)
├── surrogate/                # Pickled surrogate models
├── runs/                     # Per-run I/O directories (gitignored)
├── build/                    # Compiled binaries (gitignored)
├── logs/                     # SLURM stdout/stderr logs (gitignored)
├── Makefile                  # make build / build-video / run / submit / clean
└── pyproject.toml            # Python dependencies (activate .venv before use)
```

---

## 11. Test suite

```bash
source .venv/bin/activate

# Fast unit tests only (no binary, no SLURM — seconds)
pytest tests/ -m "not medium and not hpc"

# Include numerical health tests (runs binary at fidelity 3, ~5–10 min)
pytest tests/ -m "not hpc"

# Full suite including SLURM smoke tests (needs cluster allocation)
pytest tests/
```

---

## References

- Kim M., Harris D.M., Cimpeanu R. (2025). *Modelling of oxygen transfer in rocking bioreactors.* Int. J. Multiphase Flow. [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538)
- Basilisk CFD: http://basilisk.fr
