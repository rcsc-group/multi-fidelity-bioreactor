# Fidelity guide

Runtimes are per segment (one SLURM job) with the sweep t_end of ~35–98 nondim units.

| Fidelity | Grid | 4 CPUs | 16 CPUs | Use for |
|----------|------|--------|---------|---------|
| 3 | 8×8 | seconds | — | Smoke tests, CI |
| 4 | 16×16 | ~2 min | — | Quick debugging |
| 5 | 32×32 | ~2 min | ~1 min | BO DoE, low-fidelity surrogate |
| 6 | 64×64 | ~15 min | ~5 min | Mid-fidelity check |
| 7 | 128×128 | ~2 h | ~30 min | Standard production sweeps |
| 8 | 256×256 | ~16 h (approx.) | ~4 h (approx.) | Grid-convergence check vs. fidelity 7 |
| 9 | 512×512 | ~days | ~10–11 h/condition (measured) | High-fidelity reference |
| 10 | 1024×1024 | — | ~2.7–6.6 days/condition (measured) | Publication (Kim et al. 2024/2025) |

Specify `"cpus": 16` in `_sweep` when using fidelity ≥ 7.

For which level to actually pick, see [How to choose a fidelity level](../how-to/choose-fidelity.md).
