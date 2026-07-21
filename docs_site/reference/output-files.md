# Output files reference

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
