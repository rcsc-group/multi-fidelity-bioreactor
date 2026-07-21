# How to generate videos

Videos are produced automatically during every simulation run that uses
`BioReactor-video` — there's no separate rendering step to run.

## Build the video-capable binary

```bash
make build-video    # compiles build/BioReactor-video (only needed if stale)
```

`ffmpeg` and Basilisk's `ppm2mp4` helper are loaded automatically by the
SLURM template (`module load ffmpeg`) — no additional setup required.

## Use it

Any workflow that submits via SLURM (chain, sweep, sample, or the BO loop)
uses `BioReactor-video` already if you set `"videos": true` in a chain
config, or if the params include the video flags the sweep runner sets. For
a one-off local run, point directly at the binary:

```bash
build/BioReactor-video runs/my_run/params.json
```

## Output

Four MP4s land in `runs/<run_id>/` alongside the data files:

| File | Content |
|------|---------|
| `vorticity3.mp4` | Vorticity field (body frame) |
| `volume_fraction3.mp4` | VOF interface (body frame) |
| `oxygen3.mp4` | Dissolved oxygen concentration |
| `tracer.mp4` | Tracer mixing (top-half injection) |

See [Output files reference](../reference/output-files.md) for everything else written to `runs/<run_id>/`.
