# Rocking Bioreactor 2D — Simulation Suite

Two-phase CFD solver for a rocking bioreactor, implemented in [Basilisk](http://basilisk.fr/).
Developed at the Harris Lab (Brown University) in collaboration with the Cimpeanu group (Warwick).
Publication: [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538) | preprint: [arXiv: 2504.05421](https://arxiv.org/abs/2504.05421)

---

## Quick start (OSCAR / HPC)

### 1. Build the binary

```bash
cd dev/rocking-bioreactor-2d
make build
```

The Makefile automatically uses `~/scratch/basilisk/src/qcc` (the scratch-compiled
qcc, which works on OSCAR — the spack module has a broken header path).

### 2. Create a params.json

Copy the minimal template below into `runs/my_run/params.json` and edit as needed:

```json
{
  "run_id":         "my_run",
  "fidelity":       7,
  "t_end":          250.0,
  "omega_b":        3.93,
  "n_harmonics":    1,
  "theta_max":      [7.0, 0.0, 0.0],
  "phi_angular":    [0.0, 0.0, 0.0],
  "omega_h":        0.0,
  "amplitude_h":    [0.0, 0.0, 0.0],
  "phi_horizontal": [0.0, 0.0, 0.0],
  "geometry":       {"a": 0.25, "b": 0.071, "n": 8.0},
  "fill_level":     0.5
}
```

### 3. Run locally (quick test)

```bash
make run PARAMS=runs/my_run/params.json
```

Outputs land in `runs/my_run/`. Use `fidelity: 3` or `4` for fast (~minutes) test runs.

### 4. Submit to SLURM (production run)

```bash
make submit PARAMS=runs/my_run/params.json
```

Or via Python (also waits for results):

```bash
.venv/bin/python scripts/simulate.py runs/my_run/params.json --slurm --wait
```

### 5. Extract kLa

```bash
.venv/bin/python scripts/postprocess.py runs/my_run/
```

Reads `tr_oxy.dat` + `vol_frac_interf.dat`, writes `runs/my_run/results.json`:

```json
{"kLa_10": 0.0031, "kLa_25": 0.0028, "kLa_50": 0.0024}
```

---

## params.json reference

| Field | Units | Default | Description |
|-------|-------|---------|-------------|
| `run_id` | — | required | Unique name; output goes to `runs/{run_id}/` |
| `fidelity` | — | required | Basilisk LEVEL (grid = 2^fidelity per side). 3–4 = fast test, 7 = production, 9 = high-res |
| `t_end` | non-dim | 250.0 | Simulation end time. 1 unit ≈ T_bio = L_bio/U_bio seconds |
| `omega_b` | rad/s | required | Fundamental **angular** rocking frequency (e.g. 3.93 ≈ 37.5 RPM) |
| `n_harmonics` | int | 1 | Number of active harmonics (1–3). Vectors always padded to length 3 |
| `theta_max` | degrees | [7,0,0] | Angular amplitude of each harmonic. `theta_max[0]` is fundamental |
| `phi_angular` | rad | [0,0,0] | Angular phase of each harmonic. `phi_angular[0]` is always forced to 0 |
| `omega_h` | rad/s | 0.0 | Fundamental **horizontal** translation frequency. Set to 0 to disable |
| `amplitude_h` | m | [0,0,0] | Horizontal amplitude of each harmonic |
| `phi_horizontal` | rad | [0,0,0] | Horizontal phase of each harmonic |
| `geometry.a` | m | 0.25 | Bag half-width (horizontal semi-axis) |
| `geometry.b` | m | 0.071 | Bag half-height (vertical semi-axis) |
| `geometry.n` | — | 8.0 | Superellipse exponent. 2 = ellipse, n >= 8 = rectangle |
| `fill_level` | fraction | 0.5 | Fraction of bag volume filled with liquid (0.3–0.7) |

### Note on omega_b vs omega_h

`omega_b` and `omega_h` are **independent** in the current model.
On a physical rocking platform they are driven by the same motor and are
therefore synchronized (ω_h = ω_b). For pure angular rocking (no translation)
set `omega_h: 0.0` and all `amplitude_h` to zero.

### Note on phi_angular[0]

The phase of the fundamental angular harmonic is a global time-origin shift —
it is physically redundant and always forced to 0 at read time regardless of
what is written in params.json.

---

## Output files (in `runs/{run_id}/`)

| File | Columns | Description |
|------|---------|-------------|
| `vol_frac_interf.dat` | `i t f_liq_sum f_liq_interf posY_max posY_min` | Liquid volume, interface length, interface y-extent |
| `tr_oxy.dat` | `i t oxy_liq_sum oxy_liq_sum2 c_liq_sum ...` | Total dissolved O2 and tracer in liquid |
| `normf.dat` | `i t Omega_avg Omega_rms ... ux ... uy ...` | Vorticity and velocity norms in liquid |
| `logstats.dat` | text log | Timestep, cell count, wall-clock time |
| `results.json` | JSON | kLa_10, kLa_25, kLa_50 (written by postprocess.py) |

All `.dat` files have a one-line header (starts with `i`). Time `t` is
non-dimensional (t_physical = t × T_bio).

---

## Fidelity guide

| fidelity | grid | typical runtime | use for |
|----------|------|-----------------|---------|
| 3 | 8×8 | seconds | import test only |
| 4 | 16×16 | ~2 min | quick parameter sweep |
| 6 | 64×64 | ~30 min | low-fidelity surrogate data |
| 7 | 128×128 | ~2–4 h | standard production run |
| 9 | 512×512 | ~days | high-fidelity reference |

---

## Project structure

```
rocking-bioreactor-2d/
├── src/
│   ├── BioReactor.c          # Basilisk solver (edit compile-time flags here)
│   ├── params_read.h         # JSON → BioreactorParams struct (jsmn-based)
│   ├── henry_oxy2.h          # Henry's law oxygen transport
│   ├── utils2.h, view3.h     # Basilisk visualization helpers
│   └── jsmn.h                # Minimal JSON parser (MIT)
├── scripts/
│   ├── simulate.py           # run_local() / submit_slurm() / wait_for_result()
│   └── postprocess.py        # kLa extraction from tr_oxy.dat → results.json
├── config/
│   ├── slurm_template.sh     # SBATCH script (uses $PARAMS env var)
│   └── param_space.yaml      # Bounds for optimization suite
├── tests/
│   ├── conftest.py           # Shared fixtures (run_bioreactor, loaders, CANONICAL_PARAMS)
│   ├── integration/          # Smoke test: binary produces output files
│   ├── verification/         # Numerical health: mass conservation, O2 trend, interface
│   ├── test_postprocess.py   # kLa extraction unit tests
│   └── test_simulate.py      # simulate.py unit tests (mocked sbatch)
├── runs/                     # All simulation I/O (gitignored)
├── build/                    # Compiled binary (gitignored)
├── Makefile                  # make build / make run / make submit
└── pyproject.toml            # Python deps; activate .venv before using scripts
```

---

## Running the test suite

```bash
# activate the project venv first
source .venv/bin/activate

# fast tests only (no CFD)
pytest tests/ -m "not medium and not hpc"

# include numerical health tests (~7 min, runs the binary at fidelity=3)
pytest tests/ -m "not hpc"
```

---

## References

- Kim M., Harris D.M., Cimpeanu R. (2025). *Modelling of oxygen transfer in rocking bioreactors.* Int. J. Multiphase Flow. [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538)
- Basilisk CFD: http://basilisk.fr
