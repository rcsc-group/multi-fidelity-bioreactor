# Your first simulation

This walks through running one simulation from a cold start to a finished
`results.json`, using the real binary at fidelity 3 (fast enough to run on a
login node in under a minute). Every command and output below was actually
run to write this page — you should see the same shape of output, though
your exact CPU-time numbers will differ.

## 1. Build

```bash
uv sync
make build
```

`make build` compiles `build/BioReactor` from `src/BioReactor.c` via
Basilisk's `qcc`. If this is your first build, see [Setup](../setup.md) —
you need Basilisk's `qcc` on your `PATH` first.

## 2. Write a params.json

Fidelity 3 (8×8 cells) is too coarse to trust physically, but that's the
point of this tutorial — it's about the *pipeline*, not the physics. We'll
also cut `n_mix_cycles` down to 5 (from a normal 80) so oxygen injection
starts almost immediately, and `t_end` down to 40, so the whole run finishes
in seconds instead of minutes.

```bash
mkdir -p runs/tutorial_demo
cat > runs/tutorial_demo/params.json <<'EOF'
{
  "run_id": "tutorial_demo",
  "fidelity": 3,
  "omega_b": 3.93,
  "n_harmonics": 1,
  "theta_max": [7.0, 0.0, 0.0],
  "phi_angular": [0.0, 0.0, 0.0],
  "omega_h": 0.0,
  "amplitude_h": [0.0, 0.0, 0.0],
  "phi_horizontal": [0.0, 0.0, 0.0],
  "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
  "fill_level": 0.5,
  "n_mix_cycles": 5,
  "t_end": 40.0
}
EOF
```

## 3. Run it

```bash
cd runs/tutorial_demo
../../build/BioReactor params.json
```

Expected tail of output:

```
checkpoint: writing checkpoint.dump at t=40.13

# Quadtree, 7217 steps, 12.7809 CPU, 14.1 real, 3.28e+04 points.step/s, 85 var
```

`14.1 real` is wall-clock seconds. You should now have six files sitting
next to `params.json`:

```
checkpoint.dump  logstats.dat  normf.dat  shear_stress.dat  tr_oxy.dat  vol_frac_interf.dat
```

See the [output files reference](../reference/output-files.md) for what each one contains.

## 4. Postprocess

```bash
cd ../..
uv run python scripts/postprocess.py runs/tutorial_demo/
```

This writes `runs/tutorial_demo/results.json`. Ours came out to:

```json
{
  "kLa_10": 3448.06, "kLa_25": 1762.12, "kLa_50": 552.84,
  "kLa_inst_10": 3928.23, "kLa_inst_25": 1396.43, "kLa_inst_50": 1076.45,
  "dtmix_0.50": 0.210, "dtmix_0.75": 0.263, "dtmix_0.95": 5.469,
  "vor_mean": 1.763, "vel_rms_qss": 0.772, "kla_fit_rmse_25": 0.0077,
  "tau_95_qss": 0.00283, "tau_98_qss": 0.00327, "tau_100_qss": 0.00394,
  "tau_95_max": 0.00405, "tau_98_max": 0.00462, "tau_100_max": 0.00575,
  "tau_mean_max": 0.00194
}
```

## 5. See it

Fidelity 3 is too coarse to look at — 8×8 cells barely resolves the
interface. Bumping to fidelity 5 (32×32) and using `BioReactor-video`
instead makes the sloshing actually visible, at the cost of ~1 minute
instead of ~15 seconds:

```bash
make build-video
mkdir -p runs/tutorial_video_demo
cat > runs/tutorial_video_demo/params.json <<'EOF'
{
  "run_id": "tutorial_video_demo",
  "fidelity": 5,
  "omega_b": 3.93,
  "n_harmonics": 1,
  "theta_max": [7.0, 0.0, 0.0],
  "phi_angular": [0.0, 0.0, 0.0],
  "omega_h": 0.0,
  "amplitude_h": [0.0, 0.0, 0.0],
  "phi_horizontal": [0.0, 0.0, 0.0],
  "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
  "fill_level": 0.5,
  "n_mix_cycles": 8,
  "t_end": 20.0
}
EOF
build/BioReactor-video runs/tutorial_video_demo/params.json
uv run python scripts/render_videos.py runs/tutorial_video_demo
```

`BioReactor-video` itself only dumps raw binary frames to
`runs/tutorial_video_demo/frames/` — `render_videos.py` is the separate step
that actually renders and encodes them (needs `ffmpeg` on `PATH`; `module
load ffmpeg` on OSCAR), producing `volume_fraction.mp4` (body frame, rocking
with the bag) and `volume_fraction_lab.mp4` (lab frame, fixed camera):

![Body-frame volume-fraction animation from a real fidelity-5 run: the liquid (dark) sloshing as the bag rocks back and forth.](../assets/img/first-simulation-fidelity5.gif)

This is the same VOF field that `vol_frac_interf.dat` records numerically —
the video is just that field rendered frame by frame, nothing the solver
computes differently.

!!! note "There's a second, separate video pipeline"
    `config/slurm_video_template.sh` + `scripts/submit_video_run.py` render
    videos a different way — directly from Basilisk's own view/`bview`
    output via `ppm2mp4`, producing `vorticity3.mp4`/`oxygen3.mp4`/`tracer*.mp4`
    instead of `volume_fraction*.mp4`. That path hasn't been exercised while
    writing this page; the steps above are the ones actually verified here.

## What you just exercised

`BioReactor` read `params.json`, ran the two-phase VOF solver with Henry's-law
oxygen transport, wrote its state to `.dat` files as it went, dumped a
`checkpoint.dump` at the end (this matters once you get to
[checkpoint restart](../explanation/checkpoint-restart.md)), and
`postprocess.py` reduced those raw files down to the KPIs in
`results.json`. Every other workflow in this project — sweeps, batch
sampling, the BO loop — is this same run → postprocess step, automated and
repeated.

## Next

- [Your first sweep](first-sweep.md) — chain several of these together with checkpoint restart
- [Your first optimization loop](first-optimization-loop.md) — no Basilisk build required
