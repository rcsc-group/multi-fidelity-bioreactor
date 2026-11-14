# Meeting 2026-10-15

We rebuilt this folder on 2026-10-21 from the git history. The figures are the versions of 14 October.

| No. | File | What it showed | What the PIs said |
|---|---|---|---|
| 01 | `01_fig_numdiff_share.png` | Share of physical transport in the tracer and oxygen mixing, for each grid level | They did not accept it. They called it a black box. |
| 02 | `02_fig_numdiff_dscan.png` | Mixing time and kLa against the scalar diffusivity, levels L6 to L9 | Same. |
| 03 | `03_fig_variance_budget.png` | Tracer variance budget for each level | Same. |
| 04 | `04_fig_h3_rpm.png` | Shear and dissipation against rocking speed. Chain L8, L9, L10. | They asked for more L10 runs and a real test set. They found the uncertainty too small. They were right. |
| 05 | `05_fig_h3_parity.png` | Predicted against simulated values | We dropped it. |

Open questions from the PIs:
- Does our postprocessing match Kim's? (Yes. See the 2026-10-22 meeting, card 11.)
- What is the meaning of the variance ratio r?
- Numerical diffusion can be useful. Is it a problem?
- kLa looks random at L10.

We withdrew the numerical diffusion claim after this meeting.
