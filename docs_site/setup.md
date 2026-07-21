# Setup

## 1. Install Basilisk

The Makefile needs `qcc` (Basilisk's C preprocessor) somewhere on `PATH`, or
`$BASILISK` pointing at its `src/` directory. Basilisk has no version tags —
this builds whatever the current tarball is:

```bash
wget https://basilisk.fr/basilisk/basilisk.tar.gz
tar xzf basilisk.tar.gz
cd basilisk/src
ln -s config.gcc config
make -k
make
export BASILISK="$PWD"
export PATH="$PATH:$BASILISK"
```

Only requirement per Basilisk's own docs: a C99-compliant compiler and a
GNU-make-compatible `make`. `gawk` specifically (not `mawk`, Debian/Ubuntu's
default) is also needed — Basilisk's build scripts are awk-based.

!!! note "Building `build-video` also needs headless GL libraries"
    `make build-video` links against Basilisk's `gl/` utilities, which need
    `libGL`/`libGLU` headers to link even for headless (no-display) frame
    dumps. On Debian/Ubuntu: `apt install libgl1-mesa-dev libglu1-mesa-dev`.

!!! note "Basilisk's grid/gpu, grid/cuda, grid/hip, grid/opencl may fail to build — that's expected"
    On a machine with no CUDA/HIP/OpenCL toolchain, those four subdirectories'
    `Makefile.deps` regeneration fails and makes `make -k`/`make` exit
    non-zero, even though everything this CPU-only solver actually needs
    (`qcc`, `libglutils.a`, `libfb_tiny.a`) builds fine regardless. Verify
    success by checking those specific files exist under `basilisk/src/` and
    `basilisk/src/gl/`, not by trusting `make`'s exit code.

### If you're on Brown's OSCAR cluster specifically

The project's own SLURM jobs use a Basilisk build under persistent storage,
not `~/scratch` — `/oscar/scratch` purges files unmodified for ~30 days, and
since Basilisk is only ever read from, never modified, it can silently look
stale and get swept. The OSCAR spack module (`module load basilisk`) has a
broken header path and shouldn't be used either way. If you're extending
this project's own production sweeps (not just building it yourself),
coordinate with whoever maintains that cluster-specific install rather than
building a second, separate copy.

## 2. Build the simulation binary

```bash
make build           # standard kLa-only binary  → build/BioReactor
make build-video     # + frame dumps for videos  → build/BioReactor-video
make build-health    # + Poisson diagnostics      → build/BioReactor-health
```

## 3. Install Python dependencies

```bash
pip install uv      # one-time, if not already installed
uv sync
```

Run all scripts via `uv run python scripts/foo.py ...` or `uv run pytest`.
`uv` manages the virtual environment automatically — no manual activation needed.

## Building this documentation site

```bash
uv run --group docs mkdocs serve   # live-reloading local preview
uv run --group docs mkdocs build   # static site → site/
```
