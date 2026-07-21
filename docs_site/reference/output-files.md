# Output files reference

All files land in `runs/{run_id}/`.

| File | Written by | Columns / content |
|------|-----------|-------------------|
| `params.json` | Python scripts | Input parameters (copied in) |
| `logstats.dat` | BioReactor | `i t dt #Cells wall_time cpu_time` — one row per 0.1 non-dim time |
| `normf.dat` | BioReactor | `i t Omega_liq_avg Omega_liq_rms ... ux ... uy ...` — vorticity and velocity norms in liquid phase |
| `vol_frac_interf.dat` | BioReactor | `i t f_liq_sum f_liq_interf posY_max posY_min` — liquid volume, interface length, interface y-extent |
| `tr_oxy.dat` | BioReactor | `i t oxy_liq_sum oxy_liq_sum2 [c_liq...] c2_liq_sum c2_liq_sum2 ...` — dissolved O₂ and tracer integrals; written from `t_mix` onward |
| `results.json` | postprocess.py | 19 KPIs — see the table below |
| `checkpoint.dump` | BioReactor | Binary Basilisk dump at the end of each run (for chain restart) |
| `frames/` | BioReactor-video | Raw binary frame dumps (grid + VOF field per timestep) — temporary; `render_videos.py` consumes and deletes this directory |
| `volume_fraction.mp4` | render_videos.py | VOF interface animation, body frame (rotates with the bag) |
| `volume_fraction_lab.mp4` | render_videos.py | Same field, lab frame (fixed camera) |

All `.dat` files have a one-line header beginning with `i`.
Time `t` is non-dimensional; `t_physical = t × T_bio`.

`kLa` is `NaN` if `tr_oxy.dat` contains no data — this happens when `t_end < t_mix`
(run finished before oxygen injection started).

## results.json KPIs

All 19 keys `postprocess.py` writes, verified against a real run
(see [Your first simulation](../tutorials/first-simulation.md)):

| Key | Description | Unit |
|-----|-------------|------|
| `kLa_10/25/50` | O₂ transfer rate at 10/25/50 % saturation (5-point log-linear fit) | h⁻¹ |
| `kLa_inst_10/25/50` | Same, estimated instantaneously via dC*/dt | h⁻¹ |
| `dtmix_0.50/0.75/0.95` | Time for tracer to reach 50/75/95 % homogeneity | seconds |
| `vor_mean` | Period-averaged mean absolute vorticity (steady streaming strength) | 1/s |
| `vel_rms_qss` | RMS velocity over the quasi-steady-state window | non-dim |
| `kla_fit_rmse_25` | Fit quality (RMSE) of the kLa_25 log-linear regression | non-dim |
| `tau_95/98/100_qss` | Median shear-stress percentile over the QSS window | Pa |
| `tau_95/98/100_max` | Max shear-stress percentile over the QSS window (NOT the whole run — see [Checkpoint restart](../explanation/checkpoint-restart.md) for why that boundary matters) | Pa |
| `tau_mean_max` | Max over time of the spatially-averaged shear stress over the QSS window | Pa |

Any key can be `NaN` if its underlying data file is absent or too short —
e.g. every `tau_*` key is `NaN` if oxygen transport was disabled for the run
(and vice versa: `kLa_*` needs `tr_oxy.dat`, unrelated to whether shear
stress was tracked). See [Non-dimensionalization](../explanation/non-dimensionalization.md)
for exactly how the kLa and tau physical-unit conversions work.
