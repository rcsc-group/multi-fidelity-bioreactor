# Setup

## Build the simulation binary

```bash
cd multi-fidelity-bioreactor
make build           # standard kLa-only binary  → build/BioReactor
make build-video     # + frame dumps for videos  → build/BioReactor-video
make build-health    # + Poisson diagnostics      → build/BioReactor-health
```

The Makefile calls `/oscar/data/dharri15/eaguerov/basilisk/src/qcc` (built under
the persistent data allocation, not scratch — `/oscar/scratch` purges files
unmodified for ~30 days, and Basilisk is only ever read from, never modified,
so it silently looked stale and got swept there in July 2026. The OSCAR spack
module has a broken header path — never use `module load basilisk`).

## Install Python dependencies

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
