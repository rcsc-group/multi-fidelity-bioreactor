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

for f in ["s6_rows.png", "s7a_kim.png", "s7b_predict.png", "s8_cost.png", "s9_heatmaps.png"]:
    shutil.copy(TALK / "figs" / f, TALK / "public" / f)

S6 = """---
transition: fade
---

<div class="s6-head">our addition: a chain of three grids</div>
<img src="/s6_rows.png" class="s6-img" />

<div class="s6-eq">

$$f_{10} = \\rho_{10}\\,f_9 + \\delta_{10}, \\qquad f_9 = \\rho_9\\,f_8 + \\delta_9$$

</div>

<style>
.s6-head { position: absolute; left: 32px; top: 14px; font-size: 26px; line-height: 1.5; padding-bottom: 6px; color: var(--ink); }
.s6-img { position: absolute; left: 32px; top: 46px; width: 960px; }
.s6-eq { position: absolute; left: 0; right: 0; top: 522px; text-align: center; font-size: 24px; color: var(--ink); }
.s6-eq .katex { font-size: 1.05em; }
</style>

<!--
Our model adds one level to the multi-fidelity model. We run the same flow on three grids: L8, L9 and L10.

Each finer grid starts from the coarser one. The three fields look alike at the same rocking phase.

The cost under each level is the cost inside the chain. Costs are summed from the accounting records, per run.

L8: median 50 core-hours per run, 40 runs, all cold.

L9: median 550 core-hours per run. Ten runs are cold, at a median of 722 core-hours. Six runs are warm, at a median of 122 core-hours. The warm cost includes the L8 run it starts from.

L10: median 1,139 core-hours per run. Each L10 run is warm. Its cost includes its own run and the L9 or L8 run it starts from.

Colour is the speed of the liquid. The L8 median is 50 core-hours over the 40 runs of the chain set. The accounting audit gives 37 for all 232 cold L8 runs, which includes other sweeps and theta 2 runs.

The warm-start source of each x8_l9 run is the L8 run at the same rpm and theta (t6_l8, or t5_l8 where that is the only one). Evidence: the launch script submit_l10_from_l8.py with --level 9 resolves the source with source_run(); the restart time t_checkpoint of x8_l9_rpm25_th4 (129.04) equals the last recorded time of t6_l8_rpm25_th4 (shear_stress.dat). No t5 run exists for that point.

The model has two steps. f_9 is rho_9 times f_8, plus a Gaussian process delta_9. f_10 is rho_10 times f_9, plus delta_10. The Gaussian process carries what the coarse grid cannot show.
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
