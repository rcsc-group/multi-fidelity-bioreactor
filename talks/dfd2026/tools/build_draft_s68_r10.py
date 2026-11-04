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

for f in ["s6_row1.png", "s7a_kim.png", "s7b_predict.png", "s8_cost.png", "s9_heatmaps.png"]:
    shutil.copy(TALK / "figs" / f, TALK / "public" / f)

S6 = """---
transition: fade
---

<div class="s6-head">Extending the multi-fidelity algorithm to three levels</div>

<img src="/s6_row1.png" class="s6-row1" />

<div class="s6-c s6-c1">

$$f_8(x) = \\mathrm{KRR}(\\mathcal{D}_8)$$

</div>

<div class="s6-c s6-c2">

$$\\begin{aligned} f_9(x) &= \\beta_9 + \\rho_9\\, f_8(x) \\\\ &\\quad + \\delta_9(x) \\end{aligned}$$

</div>

<div class="s6-c s6-c3">

$$\\begin{aligned} f_{10}(x) &= \\beta_{10} + \\rho_{10}\\, f_9(x) \\\\ &\\quad + \\delta_{10}(x) \\end{aligned}$$

</div>

<div class="s6-r3 s6-r3a">40 runs &middot; 0.22 core-h/cycle</div>
<div class="s6-r3 s6-r3b">16 runs &middot; 4.10 core-h/cycle</div>
<div class="s6-r3 s6-r3c">13 runs &middot; 86.29 core-h/cycle</div>

<div class="s6-row4">

$$\\delta_\\ell \\sim \\mathcal{GP}, \\qquad w_g \\propto p(\\lambda_g)\\, p(\\eta_g)\\, p(y \\mid \\lambda_g, \\eta_g)$$

</div>

<style>
.s6-head { position: absolute; left: 32px; top: 10px; font-size: 26px; line-height: 1.5; color: var(--ink); }
.s6-row1 { position: absolute; left: 32px; top: 50px; width: 960px; }
.s6-c { position: absolute; top: 296px; width: 400px; font-size: 34px; color: var(--ink); text-align: center; }
.s6-c .katex { font-size: 1.0em; }
.s6-c1 { left: 37px; width: 283px; } .s6-c2 { left: 314px; } .s6-c3 { left: 650px; }
.s6-r3 { position: absolute; top: 468px; width: 283px; font-size: 17px; color: var(--mute); text-align: center; }
.s6-r3a { left: 37px; } .s6-r3b { left: 373px; } .s6-r3c { left: 709px; }
.s6-row4 { position: absolute; left: 0; right: 0; top: 508px; text-align: center; font-size: 24px; color: var(--ink); }
</style>

<!--
ROW3 notes: run counts and cost per cycle. Core-h per cycle = sacct cost of the run divided by its simulated cycles (T per run). Median over the runs the chain uses: L8 40 runs, median 0.22 (IQR 0.16 to 0.26); L9 16 runs (10 cold fig9 at theta 7, 6 warm), median 4.10 (IQR 2.76 to 5.83); L10 13 runs (warm, about 6 cycles each), median 86.3 (IQR 65.9 to 94.4). The L10 per-cycle value is high because each L10 run is short and includes the start-up. Source: data/row3_r9.json.

Row 2 notes. f_8 is a kernel ridge fit on the L8 runs. f_9 and f_10 use the level below as their mean basis (beta, rho), with a Gaussian process delta on top. The code: krr_lr_gpr_fullbayes.py lines 99 to 106 (basis), 116 to 142 (weights, marginal likelihood), 169 to 171 (prediction). The level below enters as a plug-in mean (h3_band_check.py lines 23 to 27).

Row 4: the grid over length scale lambda and noise eta. Weights are the prior times the marginal likelihood p(y given lambda, eta). Full forms: K_g = k_lambda + eta^2 I (line 126); a_n = a0 + (n - p)/2 and b_n = b0 + r^T K^-1 r / 2 (line 137); the log marginal likelihood at line 140.

Left-to-right: L8, then L9 (warm-started from L8), then L10 (warm-started from L9). The zoom inset shows the same window at each level (32 L8 cells across), with cell edges drawn. The window is where |u_L8 - u_L10| is largest.
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
x axis: training cost in L10-run equivalents. One unit is the measured cost of one L10 run of our code, 666.7 core-h. Each curve is core-h divided by that unit. The Kim code is divided by its assumed 120 core-h per run.

Warm-start accounting. Every chain level is warm-started from the level below. The ten cold fig9 theta-7 L9 runs are costed as warm L9 runs. Their cost is estimated as the warm per-cycle median (4.65 core-h/cycle, 45 warm L9 runs) times the post-release cycles plus 12 warm-up cycles. The estimate is for 10 runs, 7,353 core-h in total. The 6 measured warm x8_l9 runs are unchanged.

Per-run rule (same for all curves): each run is charged its own cost plus its warm-start source when the source is not already in the method's set.
- Approved chain (full set): 40 L8 runs, 16 L9 runs, 13 L10 runs. Base 9,841 core-h, then 667 core-h per L10 run.
- Cheap chain variants (10 or 20 L8 runs, 6 warm L9 runs): base 918 to 2,488 core-h, then 1,091 core-h per L10 run. The L10 runs bring their L9 sources, which are estimated.
- H1: 40 L8 runs (2,043 core-h), then 1,091 core-h per L10 run.
- L10 GP, ours: 1,113 core-h per L10 run, including sources.
- Kim code: 120 core-h per run, no sources (assumed).

Crossings at 10 % error (relative RMS), in L10-run equivalents (core-h in brackets):
- tau 95: the chain envelope reaches 10 % at 6.3 (4,189 core-h), with 3 L10 runs and the cheap variant (warm9_L8_10). The approved chain reaches it at 17.8 (11,841). H1 reaches it at 22.7 (15,129). The L10 GP does not reach it by 12 runs.
- EDR: only the approved chain reaches 10 %, at 17.8 (11,841 core-h).

Crossings at 20 % error:
- tau 95: the chain envelope at 6.3 (4,189 core-h). H1 at 8.0 (5,314). L10 GP at 16.7 (11,133). The Kim code at 10.0 (1,200 core-h, 10 runs). So the Kim code is cheaper than the chain at this error level. This comparison uses the assumed 120 core-h per Kim run, not a measured cost.
- EDR: only the chain reaches 20 %, at 17.8 (11,841 core-h).

Cheap variants: they reach the 10 % tau level at lower cost, but they fail on EDR (errors 24 to 1,080 %). The envelope keeps them, because they are cheaper on tau_95.
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
