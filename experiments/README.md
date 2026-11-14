# experiments/ — results of numerical experiments

Small result files (JSON, CSV, a few figures) that the diary and the figure scripts read. Run data
(snapshots, dumps, logs) is not here; it is under `runs/` on the cluster.

| Path | Contents |
|---|---|
| `multifidelity/` | Surrogate and optimisation results: one `test_*.json` (and `.png`) per numbered test in `diary.md`. `kp_parity/`: comparison with the Kim et al. reference. `guide/`: figures for the docs. |
| `kimetal2024/` | Reference data from Kim, Harris and Cimpeanu (2025): raw CSV, the paper source and figures, and our replicas of their figures (`figure_replicas/`). |
| `canonical_case/` | One reference run: parameters, results and two videos. |
| `sweep_fb_*`, `l9_l10_*`, `geometry_*`, `kim_one_case_*` | Early sweeps and level checks (theta, fill, rate). |
| `settling_*`, `dtheta_*`, `*_manifest.json` | Settling-time and restart studies: the manifest lists the runs, the results file the measured values. |
| `hypothesis_ledger.json`, `run_provenance.csv` | Index of tested hypotheses; which binary and parameters produced each run. |
| `docs/` | Older internal notes (MPI checkpoint development, problem statement). |

## How to find the source of a result

1. Search `diary.md` (and `diary_archive_to_2026-10-07.md`) for the file name. The entry gives the
   question, the command, the result and the conclusion.
2. The script that writes the file usually has the same stem: `experiments/multifidelity/test_t4_retro_bo.json`
   comes from a script named like `scripts/*t4*`. Search with `git grep -l test_t4_retro_bo scripts/`.
