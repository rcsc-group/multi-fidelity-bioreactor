"""Build slides_s68draft.md = slides.md + new slides 6-8 inserted before '# Results 1'. Round 2.

Reads slides.md (unchanged), copies figs/*.png to public/, writes slides_s68draft.md only.
Numbers on slides: none (figures carry the numbers, speaker notes give the cost definitions).
"""
import shutil
from pathlib import Path

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
src = (TALK / "slides.md").read_text()
anchor = "-->\n\n---\n\n# Results 1"
assert src.count(anchor) == 1, src.count(anchor)

for f in ["s6_chain.png", "s7a_kim.png", "s7b_predict.png", "s8_cost.png"]:
    shutil.copy(TALK / "figs" / f, TALK / "public" / f)

S6 = """---
transition: fade
---

<div class="s6-head">our addition: a chain of three grids</div>
<img src="/s6_chain.png" class="s6-img" />

<div class="s6-eq">

$$f_{10} = \\rho_{10}\\,f_9 + \\delta_{10}, \\qquad f_9 = \\rho_9\\,f_8 + \\delta_9$$

</div>

<style>
.s6-head { position: absolute; left: 32px; top: 22px; font-size: 26px; color: var(--ink); }
.s6-img { position: absolute; left: 32px; top: 40px; width: 960px; }
.s6-eq { position: absolute; left: 0; right: 0; top: 400px; text-align: center; font-size: 28px; color: var(--ink); }
.s6-eq .katex { font-size: 1.05em; }
</style>

<!--
Our model adds one level to the multi-fidelity model. We run the same flow on three grids: L8, L9 and L10.

Each finer grid starts from the coarser one. The three fields look alike at the same rocking phase.

The cost under each level is the cost inside the chain. Costs are summed from the accounting records, per run.

L8: median 50 core-hours per run, 40 runs, all cold.

L9: median 550 core-hours per run. Ten runs are cold, at a median of 723 core-hours. Six runs are warm, at a median of 122 core-hours. The warm cost includes the L8 run it starts from.

L10: median 1,139 core-hours per run. Each L10 run is warm. Its cost includes its own run and the L9 or L8 run it starts from.

The model has two steps. f_9 is rho_9 times f_8, plus a Gaussian process delta_9. f_10 is rho_10 times f_9, plus delta_10. The Gaussian process carries what the coarse grid cannot show.
-->
"""

S7 = """---
transition: fade
---

<img src="/s7a_kim.png" class="s7-a" />
<img src="/s7b_predict.png" class="s7-b" />

<style>
.s7-a { position: absolute; left: 47px; top: 56px; width: 930px; }
.s7-b { position: absolute; left: 47px; top: 318px; width: 930px; }
</style>

<!--
Top: our L10 runs and the digitised Kim et al. 2025 curves at 7 degrees. Both statistics follow Kim. The 22.5 rpm bump is in both.

At 25 rpm there are two L10 runs. The filled point is the training run. The open point is the second run. They differ by about 8 % in the EDR. We show both, so the spread is visible.

Bottom: the model predicts each run before it sees that run. The bars show the 95 % band.

Left: the 13 training runs, each left out once.

Right: six new runs at 5.5 degrees. Their predictions use no refit.

The EDR is the maximum over time, as Kim defines it. At the new angle, the time mean is more accurate than this maximum.
-->
"""

S8 = """---
transition: fade
---

<img src="/s8_cost.png" class="s8-img" />
<div class="s8-cap">Cost is the summed accounting cost of each training set. The L8 and L9 grids set a floor near 11,000 core-h.<br/>The L9 log-linear ratio reaches 5 % at lower cost than the chain for tau 95. For EDR, only the log-linear ratio reaches 5 %.</div>

<style>
.s8-img { position: absolute; left: 42px; top: 74px; width: 940px; }
.s8-cap { position: absolute; left: 32px; right: 32px; top: 462px; font-size: 18px; line-height: 1.3; color: var(--ink); }
</style>

<!--
Each curve shows error against training cost, in core-hours. The axis is logarithmic.

Cost means the accounting cost of each run. It is the elapsed time times the number of cores, summed over all jobs of that run.

The chain needs all 40 L8 runs and all 16 L9 runs, plus their sources. That base is about 11,200 core-hours. Each L10 run adds about 670 core-hours on average.

The L9 ratios need the 16 L9 runs and their L8 sources. That base is about 9,500 core-hours. The log-linear ratio has the same cost as the constant ratio.

The L10-only model needs no L9 grid. Each L10 run brings its own warm-start source. Each run then costs about 1,240 core-hours on average.

Reading the curves: at the 10 % level, the constant L9 ratio is cheaper than the chain for both quantities. At 5 %, the log-linear L9 ratio is cheaper than the chain for tau 95. For EDR, only the log-linear ratio reaches 5 % within 12 runs.
-->
"""

text = src.replace(anchor, "-->\n\n" + S6 + "\n" + S7 + "\n" + S8 + "\n---\n\n# Results 1", 1)
(TALK / "slides_s68draft.md").write_text(text)
print("wrote", TALK / "slides_s68draft.md")
