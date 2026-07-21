# Workflow A — Single run

Use this for one-off simulations, debugging, or manual parameter exploration.

## Step 1 — Write a params.json

Create `runs/my_run/params.json` (see the [params.json reference](../reference/params.md) for all fields):

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

## Step 2 — Run

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

## Step 3 — Postprocess

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

## Alternatively — set up the directory without submitting

`launch.py` writes the run directory and a SLURM script but does **not** submit:

```bash
python scripts/launch.py path/to/params.json [runs_root]
# returns {"run_id": "...", "run_dir": "...", "slurm_script": "..."}
```
