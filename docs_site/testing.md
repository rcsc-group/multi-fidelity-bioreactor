# Test suite

`pyproject.toml` sets `addopts = -m 'not medium and not hpc'`, so the default
invocation runs only the fast, hermetic unit tests. `medium` and `hpc` tests are
opt-in — pass an explicit `-m` to include them (this overrides addopts, it does
not combine with it).

```bash
# Fast unit tests only (no binary, no SLURM — ~20 s). This is the default.
uv run python -m pytest tests/

# Include numerical verification tests (runs the real binary at fidelity 3, ~2 min each)
uv run python -m pytest tests/ -m medium

# Full suite including SLURM integration tests (submits real jobs, needs cluster allocation)
uv run python -m pytest tests/ -m "medium or hpc"
```

!!! note "Run `medium`/`hpc` via CI, not locally"
    The `medium` and `hpc` marked tests are wired to run on GitHub Actions,
    not routinely on OSCAR — that's the whole point of the CI setup below.
    Push your change and let CI validate it; only run these locally when you're
    specifically debugging something that requires comparing local OSCAR
    behavior against CI (e.g. isolating a toolchain-specific difference).

GitHub Actions (`.github/workflows/ci.yml`) runs the fast suite on every push/PR,
and the `medium` numerical-verification suite (real fidelity-3 CFD runs) on
pushes to `main`. It builds Basilisk from the official source tarball into a
cached `basilisk/` directory rather than relying on any OSCAR-specific path —
`hpc`-marked tests (real SLURM submission) never run in CI, since they need an
actual cluster allocation.

**Known open issue:** the `medium` job's numerical-verification step is
currently non-blocking (`continue-on-error`). `test_interface_oscillates_at_rocking_frequency`
passes reproducibly on OSCAR — including against the exact same Basilisk
tarball CI uses, built fresh there — but fails on the CI runner's Ubuntu/gcc
toolchain (6.5% spectral power in the expected band vs. a 20% threshold,
typical correct value ~37%; not a borderline rounding difference). Basilisk
source version and simple run-to-run nondeterminism have both been ruled out
as the cause; not yet root-caused.
