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
.s6-left { position: absolute; left: 32px; top: 90px; width: 460px; height: 360px; display: flex; flex-direction: column; justify-content: center; font-size: 31px; color: var(--ink); }
.s6-left .katex { font-size: 1.0em; }
.s6-left .katex-display { margin: 14px 0; }
.s6-img { position: absolute; left: 520px; top: 90px; width: 470px; }
</style>

<!--
The recursion has three levels. Level 8 is a kernel ridge fit. Levels 9 and 10 are Gaussian processes. Each level's mean is a linear map of the level below. The weights are the prior times the marginal likelihood over a grid of length scale lambda and noise eta.

Full forms, for reference. Level 8: the length scale and noise minimise the mean squared closed-form leave-one-out residual, the mean over i of ([K^-1 y]_i / [K^-1]_ii)^2 (fast_loo_krr.py lines 14 to 18; test_h1_mfbml.py lines 101 to 107, 10 restarts). Levels 9 and 10: the kernel is k_lambda(x, x') = exp(-||x - x'||^2 / (2 lambda^2)). The grid kernel is K_g = k_lambda + eta^2 I (krr_lr_gpr_fullbayes.py line 126). The posterior shape parameters are a_n = a0 + (n - p)/2 and b_n = b0 + r^T K^-1 r / 2 (line 137). The basis is F = [1, standardised f_below] (lines 99 to 106). The below-level mean enters as a plug-in (h3_band_check.py lines 23 to 27).

The log-normal prior on lambda is at krr_lr_gpr_fullbayes.py line 79. The log marginal likelihood is at line 140. The grid weights are at line 142. The prediction is the weighted mean (lines 169 to 171).

Right: the speed along one vertical line, from the bottom wall up. The line is at x = -0.369 (domain units, origin at the centre). It runs from the wall at y = -0.143 up to y = -0.010. Distance is measured from the first liquid cell. The line goes through the same physical point at all three levels. The curves end at distance 0.13 because the line stops at the crop centre. The separation is largest within about 0.02 of the wall.

Cost labels use the slide 8 rule, applied to the full chain set. Each run is charged its own cost, plus its source cost when the source is not in the set. The labels show the mean per run. L8: 51 core-h (40 runs, cold). L9: 575 core-h (16 runs; warm sources are in the set, so only their own cost counts). L10: 667 core-h (13 runs; each source is in the set).
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
Each curve shows error against training cost, in core-hours. The axis is logarithmic. Every curve uses one accounting rule. Each run is charged its own accounting cost, plus its warm-start source when that source is not already in the method's set.

The approved chain uses 40 L8 runs and 16 L9 runs. Each L10 run then adds about 670 core-hours.

The cheaper chain variants use 10 or 20 L8 runs and 6 warm L9 runs. Their L10 runs start from cold L9 runs that are not in the set. Each L10 run then costs about 1,220 core-hours. Each variant point is one subset size, and the solid line is the lower envelope of all variants.

The L10 GP (ours) costs 1,244 core-hours per run, including the source. The H1 model costs 2,043 core-hours for the 40 L8 runs, then 1,220 core-hours per L10 run. The Kim code is assumed to cost 120 core-hours per run, with no sources. Its error values are the same as our L10 GP.

Crossings at 10 % error, relative RMS, from the cost table of this round.
- tau 95: the chain envelope reaches 10 % at 4,582 core-hours (6 warm L9 runs, 10 L8 runs, 3 L10 runs). The approved chain reaches it at 13,246. H1 reaches it at 16,700, with 12 L10 runs. The L10 GP does not reach it by 12 runs.
- EDR: only the chain reaches 10 %, at 13,246 core-hours (approved set, 3 L10 runs). The cheap variants do not.

Crossings at 20 % error.
- tau 95: the chain envelope reaches 20 % at 4,582 core-hours. H1 reaches it at 5,707. The L10 GP reaches it at 12,442. The Kim code, at the assumed 120 core-hours per run, reaches it at 1,200 core-hours (10 L10 runs). So the Kim code is cheaper than the chain at 20 %. This is an assumed cost, not a measured one.
- EDR: only the chain reaches 20 %, at 13,246 core-hours.
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
