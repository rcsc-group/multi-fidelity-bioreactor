# Meeting 2026-10-15 (reconstructed on 2026-10-21)

The figure files are the versions that existed on 2026-10-14/15, restored from git history (the later redraws in `experiments/multifidelity/` are different).

| # | File | What it showed | How computed | Outcome of the meeting |
|---|---|---|---|---|
| 01 | `01_fig_numdiff_share.png` | "physical share of transport" of the tracer and oxygen per grid level L6-L10 | share = 1/r with r = (numerical) variance loss over the loss the physical diffusivity alone gives, from the variance budget d sigma^2/dt = -2 D <grad c^2> (scripts/fig_numdiff_share.py; test V) | PIs did not accept it: black box, "diffusion only", advection dominates, numerical diffusion can be helpful. The claim was withdrawn as a headline (diary 2026-10-15). |
| 02 | `02_fig_numdiff_dscan.png` | mixing time and kLa against the scalar diffusivity D at levels L6-L9 | 22 runs, D x 1 ... 1e5 (test K, scripts/fig_numdiff_dscan.py) | same critique |
| 03 | `03_fig_variance_budget.png` | tracer variance budget per level | scripts/fig_variance_budget.py | same critique |
| 04 | `04_fig_h3_rpm.png` | 3x3: shear (tau_95, mean tau) and dissipation against rpm for theta 4, 7, 9: L8 fit, L9 fit, L10 runs, three-level chain with leave-one-out bars | chain L8 -> L9 -> L10, log of the quantity (scripts/fig_h3_rpm.py, scripts/test_h3_mfbml.py) | PI asked for more L10 runs and a real test set; the uncertainty was "surprisingly low" (it was: it ignored the L9 uncertainty, fixed 2026-10-19) |
| 05 | `05_fig_h3_parity.png` | predicted against simulated parity plot | scripts/fig_h3_parity.py | rejected by the user ("not interested in simulated against predicted") |

Questions the PIs raised: parity with Kim's postprocessing (answered 2026-10-15/16: same jumps in Kim's code), variance equation "diffusion only", the meaning of r, variance loss at the interface, numerical diffusion can be helpful, kLa looks random at L10.
