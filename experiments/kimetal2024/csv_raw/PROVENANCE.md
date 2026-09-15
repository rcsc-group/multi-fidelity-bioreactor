# Provenance of the Kim et al. (2024) reference CSVs

Recorded 2026-09-15, after three separate incidents this week where a stored
artifact was trusted on its filename alone. Check this file before using any
CSV here as a reference series.

| file | source | status |
|---|---|---|
| `mixing_kla_vs_frequency.csv` | `dataset_kimetal2024.xlsx`, sheet `Data`, block labelled **"Fig. 9,11"** (rows 3-12) | **Author's own published values.** Verified bit-exact against the xlsx (max abs diff 0.0 across all 10 rows x 6 columns). Not digitized. |
| `mixing_kla_vs_angle.csv` | same xlsx, block labelled **"Fig. 10,12"** | Author's own published values. Column is `Angle_deg` (note: *not* `theta_deg`). |
| `shear_ediss_vs_frequency.csv` | digitized from `Figures/Fig_tau_Ediss_rpm_deg.pdf` panel (a) | Digitized. **Row 1 is a units/label row, not data** -- read with `skiprows=[1]`. |
| `shear_ediss_vs_angle.csv` | digitized from the same PDF, panel (b) | Digitized 2026-09-15; the digitiser reproduces the rpm CSV to 0.1% median error. No units row. Column is `theta_deg`. |

Figs 9-12 are **not** in the paper's figure set as digitizable curves we need
to reconstruct -- Kim published the underlying numbers. Fig 13 is not in the
xlsx, which is why it had to be digitized.

`replicate_plots.py` redraws these CSVs and reads **zero** of our runs. Its
output belongs in `../kim_redrawn/`, never in `../figure_replicas/`.

Units: `dtmix_strict_*` are **seconds**; `vor_meanabs_steady_streaming` is
1/s; `kLa_*` are 1/h.
