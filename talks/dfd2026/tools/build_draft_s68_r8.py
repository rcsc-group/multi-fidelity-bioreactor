"""Build slides_s68draft.md = slides.md + slides 6, 7, 8, 9 inserted before '# Results 1'. Round 3.

Slide 6: heading fixed (not clipped). Slide 7: unchanged from round 2. Slide 8: three methods, no L9 ratio curves.
Slide 9: 2D campaign heatmaps (template, baseline data). Reads slides.md (unchanged), writes only slides_s68draft.md.
"""
import shutil
from pathlib import Path

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
src = (TALK / "slides.md").read_text()
anchor = "-->\n\n---\n\n# Results 1"
assert src.count(anchor) == 1, src.count(anchor)

for f in ["s6_profile.png", "s7a_kim.png", "s7b_predict.png", "s8_cost.png", "s9_heatmaps.png"]:
    shutil.copy(TALK / "figs" / f, TALK / "public" / f)

S6 = """---
transition: fade
---

<div class="s6-head">Extending the multi-fidelity algorithm to three levels</div>

<div class="s6-left">

$$f_\\ell(x) = \\beta_0^{(\\ell)} + \\rho_\\ell\\, f_{\\ell-1}(x) + \\delta_\\ell(x)$$

$$\\ell = 9,\\ 10$$

$$w_g \\propto p(\\lambda_g)\\, p(\\eta_g)\\, p(y \\mid \\lambda_g, \\eta_g)$$

</div>

<img src="/s6_profile.png" class="s6-img" />

<style>
.s6-head { position: absolute; left: 32px; top: 14px; font-size: 26px; line-height: 1.5; padding-bottom: 6px; color: var(--ink); }
.s6-left { position: absolute; left: 28px; top: 200px; width: 470px; font-size: 24px; color: var(--ink); }
.s6-left .katex { font-size: 1.0em; }
.s6-left .katex-display { margin: 26px 0; }
.s6-img { position: absolute; left: 520px; top: 80px; width: 470px; }
</style>

<!--
The recursion has three levels. Level 8 is a kernel ridge fit. Level 9 and level 10 are Gaussian processes on the level below. Each level's mean is a linear map of the level below. The weights are the prior times the marginal likelihood over a grid of length scale and noise.

Level 8 is a kernel ridge fit. Its length scale and noise come from the leave-one-out objective. The objective is the mean of the squared closed-form leave-one-out residual, [K^-1 y]_i / [K^-1]_ii. The code runs 10 restarts (test_h1_mfbml.py lines 101 to 107; fast_loo_krr.py lines 14 to 18).

For levels 9 and 10, the basis is F = [1, standardised f_below] (krr_lr_gpr_fullbayes.py lines 99 to 106). The level below is a plug-in mean (h3_band_check.py lines 23 to 27); its uncertainty is not carried into the fit.

The grid over length scale lambda and noise ratio eta has a log-normal prior (krr_lr_gpr_fullbayes.py line 79). The kernel is K_g = k_lambda + eta^2 I (line 126). The posterior shape parameters are a_n = a0 + (n - p)/2 and b_n = b0 + r^T K^-1 r / 2 (line 137). The log marginal likelihood is at line 140, and the grid weights are at line 142. The prediction is the weighted mean (lines 169 to 171).

Right: the speed along the vertical line through the crop centre, from the bottom wall. The centre is L8 cell (row 125, col 33). Measured separation, relative RMS of L8 minus L10 over L10: 0.78 for the vertical line and 0.74 for the horizontal line. The vertical line was plotted. Most of the separation is within about 0.02 domain units of the wall, where the boundary layer is.

Cost labels at the curve ends are the chain-set medians: L8 50 core-h, L9 550 core-h, L10 1,139 core-h per run. The L9 figure is 10 cold runs (median 722) and 6 warm runs (median 122), both with their source. The L10 figure is the warm run plus its source.
-->
"""

S7 = Path(TALK / "tools" / "slide7_block.md").read_text()

S8 = """---
transition: fade
---

<div class="s8-head">accuracy against cost</div>
<img src="/s8_cost.png" class="s8-img" />

<style>
.s8-head { position: absolute; left: 32px; top: 14px; font-size: 26px; line-height: 1.5; padding-bottom: 6px; color: var(--ink); }
.s8-img { position: absolute; left: 42px; top: 46px; width: 940px; }
.s8-cap { position: absolute; left: 32px; right: 32px; top: 470px; font-size: 18px; line-height: 1.3; color: var(--ink); }
</style>

<!--
Each curve shows error against training cost, in core-hours. The axis is logarithmic. The points are the same random subsets for each method.

Same held-out subsets for each method. Cost is the summed accounting core-hours of each training set.

Three methods are shown. The chain uses L8, L9 and L10. The L10-only model uses L10 runs alone. The two-level model H1 uses L8 and L10 directly.

Cost means the accounting cost of each run. It is the elapsed time times the number of cores, summed over all jobs of that run. Warm starts include the run they start from.

The chain needs all 40 L8 runs and all 16 L9 runs. That base is about 11,200 core-hours. Each L10 run adds about 670 core-hours.

H1 needs all 40 L8 runs. That base is about 2,000 core-hours. Each L10 run adds about 1,200 core-hours, because each L10 run starts from an L9 or L8 run.

Reading the curves. For tau 95, the chain reaches 10 % at about 13,200 core-hours. H1 reaches 10 % only at 12 runs, about 16,700 core-hours. For EDR, neither reaches 10 % within 12 runs.
-->
"""

S9 = """---
transition: fade
---

<img src="/s9_heatmaps.png" class="s9-img" />

<style>
.s9-img { position: absolute; left: 32px; top: 8px; width: 960px; }
</style>

<!--
This slide shows the 2D campaign as a template. The data are the baseline runs only. No new run is included.

Rows are the three gated quantities: tau 95, the energy dissipation rate, and kLa at 10 percent oxygen. Columns are the predicted median on the left, and the relative uncertainty on the right.

The uncertainty is two sigma over the mean. The black line marks 0.30, the gate for a point to count as known.

The coloured dots are the runs. Blue is L8, green is L9 and pink is L10.

The energy dissipation rows use the time mean as a stand-in. The 95th percentile (edr_95) is not yet saved for the existing runs. Its backfill jobs are queued and not finished. The final version will use edr_95.
-->
"""

text = src.replace(anchor, "-->\n\n" + S6 + "\n" + S7 + "\n" + S8 + "\n" + S9 + "\n---\n\n# Results 1", 1)
(TALK / "slides_s68draft.md").write_text(text)
print("wrote", TALK / "slides_s68draft.md")
