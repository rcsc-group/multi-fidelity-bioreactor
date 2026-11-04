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

for f in ["s6_zoom.png", "s7a_kim.png", "s7b_predict.png", "s8_cost.png", "s9_heatmaps.png"]:
    shutil.copy(TALK / "figs" / f, TALK / "public" / f)

S6 = """---
transition: fade
---

<div class="s6-head">Extending the multi-fidelity algorithm to three levels</div>

<div class="s6-left">

<div class="s6-lab">level 8: kernel ridge, leave-one-out objective</div>

$$f_8 = \\mathrm{KRR}(\\mathcal{D}_8), \\qquad \\hat\\theta_8 = \\arg\\min_\\theta\\ \\frac{1}{n}\\sum_{i=1}^{n}\\left(\\frac{[K_\\theta^{-1} y]_i}{[K_\\theta^{-1}]_{ii}}\\right)^{2}$$

<div class="s6-lab">levels 9 and 10: recursive GP, weights over (λ, η)</div>

$$f_\\ell(x) = \\beta_0^{(\\ell)} + \\rho_\\ell\\, f_{\\ell-1}(x) + \\delta_\\ell(x), \\qquad \\ell = 9,\\, 10$$

$$\\delta_\\ell \\sim \\mathcal{GP}\\!\\left(0,\\ \\sigma_\\ell^2\\, k_{\\lambda}\\right), \\qquad k_\\lambda(x,x') = \\exp\\!\\left(-\\frac{\\lVert x - x'\\rVert^2}{2\\lambda^2}\\right)$$

$$\\mu_\\ell(x) = \\sum_g w_g\\, \\mu_{\\ell,g}(x), \\qquad w_g \\propto p(\\lambda_g)\\, p(\\eta_g)\\, |K_g|^{-1/2}\\, |F^\\top K_g^{-1} F|^{-1/2}\\, b_n^{-a_n}$$

$$K_g = k_{\\lambda_g} + \\eta_g^2 I, \\qquad a_n = a_0 + \\tfrac{n-p}{2}, \\qquad b_n = b_0 + \\tfrac12\\, r^\\top K_g^{-1} r$$

</div>

<img src="/s6_zoom.png" class="s6-img" />

<style>
.s6-head { position: absolute; left: 32px; top: 14px; font-size: 26px; line-height: 1.5; padding-bottom: 6px; color: var(--ink); }
.s6-left { position: absolute; left: 24px; top: 70px; width: 470px; font-size: 14px; color: var(--ink); }
.s6-left .katex { font-size: 0.92em; }
.s6-lab { font-size: 14px; color: var(--mute); margin-top: 6px; }
.s6-img { position: absolute; left: 520px; top: 66px; width: 470px; }
</style>

<!--
The recursion has three levels. Level 8 is a kernel ridge fit. Its length scale and noise come from the leave-one-out objective. The code calls this in CachedLooKRR, with 10 restarts (test_h1_mfbml.py lines 101 to 107). The closed form of the leave-one-out residual is in fast_loo_krr.py lines 14 to 18.

Levels 9 and 10 are Gaussian processes on top of the level below. The mean of the lower level is the input basis: F = [1, standardised f_below] (krr_lr_gpr_fullbayes.py lines 99 to 106). The level below enters as a plug-in mean (h3_band_check.py lines 23 to 27). Its uncertainty is not carried into the fit.

The code does not minimise a likelihood over the Gaussian process hyperparameters. It evaluates a grid over the length scale lambda and the noise ratio eta. The grid prior is at krr_lr_gpr_fullbayes.py line 79. Kernel and noise are at line 126. The posterior shape parameters a_n and b_n are at line 137. The log marginal likelihood is at line 140. The grid weights are at line 142. The prediction is the weighted mean, at lines 169 to 171.

Right: the same window of the flow at the three levels. The window is 32 L8 cells wide. It sits where the L8 and L10 fields differ most. The L8 texture is blocky. The L10 field is finer. Cell edges are drawn, so the grid size can be seen.
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
