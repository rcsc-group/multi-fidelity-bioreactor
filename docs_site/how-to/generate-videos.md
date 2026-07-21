# How to generate videos

## The actual pipeline (verified)

`BioReactor-video` does **not** produce an MP4 by itself — it dumps raw
binary frames (grid + VOF field per timestep) to `run_dir/frames/`.
`render_videos.py` is the separate step that actually renders and encodes
them:

```bash
make build-video
build/BioReactor-video runs/my_run/params.json
uv run python scripts/render_videos.py runs/my_run
```

Needs `ffmpeg` on `PATH` (`module load ffmpeg` on OSCAR) and Basilisk's own
`gl` headless-rendering libs to have linked when `build-video` was compiled
(`libgl1-mesa-dev libglu1-mesa-dev` on Debian/Ubuntu — see [Setup](../setup.md)).

Produces, in `runs/my_run/`:

| File | Content |
|---|---|
| `volume_fraction.mp4` | VOF interface, body frame (rotates with the bag) |
| `volume_fraction_lab.mp4` | Same field, lab frame (fixed camera) |

`render_videos.py` cleans up `frames/` once it's done. See
[Your first simulation](../tutorials/first-simulation.md) for a real
example, including the actual rendered output.

## When this runs automatically

`chain.py` automates the whole build→run→render sequence for every segment
when a chain config sets `videos: true` — it switches to
`config/slurm_video_template.sh`, which calls `render_videos.py` itself at
the end of the SLURM job. For a one-off local run (as above), you run the
two steps yourself.

## The other video pathway (not verified here)

`scripts/submit_video_run.py` + `config/slurm_video_template.sh` exist as a
second, standalone way to submit a single video run via SLURM directly.
It calls the exact same `render_videos.py` step, so it should produce the
same `volume_fraction*.mp4` files — but that specific script hasn't been
run end-to-end while writing this page, unlike the direct-invocation path
above.
