# Video generation

Videos are generated automatically during every simulation run — no separate step
is required. All SLURM jobs use `BioReactor-video`, which renders frames inline
and encodes them to MP4 at the end of each segment.

## Output files produced per run

| File | Content |
|------|---------|
| `vorticity3.mp4` | Vorticity field (body frame) |
| `volume_fraction3.mp4` | VOF interface (body frame) |
| `oxygen3.mp4` | Dissolved oxygen concentration |
| `tracer.mp4` | Tracer mixing (top-half injection) |

All MP4s land in `runs/<run_id>/` alongside the data files.

## Build requirements

```bash
make build-video    # compiles BioReactor-video (only needed if binary is stale)
```

`ffmpeg` and Basilisk's `ppm2mp4` helper are loaded automatically by the SLURM
template (`module load ffmpeg`). No additional setup is required.
