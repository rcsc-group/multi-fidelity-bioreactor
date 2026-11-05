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

for f in ["s6_slide_r11.png", "s7a_kim.png", "s7b_predict.png", "s8_cost.png", "s9_heatmaps.png"]:
    shutil.copy(TALK / "figs" / f, TALK / "public" / f)

S6 = """---
transition: fade
---

<img src="/s6_slide_r11.png" class="s6-bg" />

<div class="s6-head">Extending the multi-fidelity algorithm to three levels</div>

<div class="s6-eq" style="top:70px">

$$f_8(x) = \\mathrm{KRR}(\\mathcal{D}_8)$$

</div>
<div class="s6-sub" style="top:150px">40 runs &middot; 0.22 core-h/cycle</div>

<div class="s6-eq" style="top:215px">

$$f_9(x) = \\beta_9 + \\rho_9\\, f_8(x) + \\delta_9(x)$$

</div>
<div class="s6-sub" style="top:295px">16 runs &middot; 5.76 core-h/cycle</div>

<div class="s6-eq" style="top:360px">

$$f_{10}(x) = \\beta_{10} + \\rho_{10}\\, f_9(x) + \\delta_{10}(x)$$

</div>
<div class="s6-sub" style="top:440px">13 runs &middot; 64.9 core-h/cycle</div>

<div class="s6-strip">

$$\\delta_\\ell \\sim \\mathcal{GP}, \\qquad w_g \\propto p(\\lambda_g)\\, p(\\eta_g)\\, p(y \\mid \\lambda_g, \\eta_g)$$

</div>

<style>
.s6-head { position: absolute; left: 32px; top: 12px; font-size: 26px; line-height: 1.5; color: var(--ink); }
.s6-bg { position: absolute; left: 0; top: 0; width: 1024px; height: 576px; }
.s6-eq { position: absolute; left: 32px; width: 450px; font-size: 30px; color: var(--ink); text-align: left; }
.s6-eq .katex { font-size: 1.0em; }
.s6-sub { position: absolute; left: 32px; width: 450px; font-size: 22px; color: var(--ink); }
.s6-strip { position: absolute; left: 0; right: 0; top: 508px; text-align: center; font-size: 26px; color: var(--ink); }
</style>

<!--
Rows, top to bottom: L8, L9, L10. Each row has the equation for that level on the left, and the cropped tank on the right.
Row 1 (L8): f_8 is a kernel ridge fit on 40 cold L8 runs. Median cost 0.22 core-h per simulated cycle (IQR 0.16 to 0.26).
Row 2 (L9): the chain trains on 16 L9 runs. Six are measured warm runs (x8_l9) at 5.76 core-h per cycle (wall-clock logs; 11.8 simulated cycles each). The other 10 are cold fig9 theta-7 runs, costed as warm estimates: 5.76 core-h per cycle times 11.8 cycles, about 68 core-h per point, 682 core-h in total. The 16-run figure includes these 10.
Row 3 (L10): 13 warm runs. Median 64.9 core-h per cycle, steady rate excluding the first two cycles. Each run simulates about 6 cycles (x8 about 12). The logs do not separate start-up from steady running, so the rate is not corrected for it.
Bottom: the level below is the mean basis (beta, rho) of the level above, with a Gaussian process delta. The weights over the length scale lambda and noise eta are the prior times the marginal likelihood. Code: krr_lr_gpr_fullbayes.py lines 99 to 106 (basis), 126 (K_g), 137 (a_n, b_n), 140 (log marginal likelihood), 142 (weights), 169 to 171 (prediction). The level below enters as a plug-in mean (h3_band_check.py lines 23 to 27).
Tanks: cropped to the liquid plus a thin gas band, equal aspect. One shared colour scale: the 10th and 90th percentiles of |u| over the pooled liquid cells of the three levels, each level weighted equally (vmin 0.0228, vmax 0.2026, nondimensional speed). Saturated liquid cells above vmax: L8 8.4 %, L9 12.0 %, L10 7.2 %. Below vmin: L8 16.6 %, L9 5.1 %, L10 9.5 %. The upper limit moved from the 85th to the 90th percentile because the 85th left 19 % of L9 saturated.
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
x axis: training cost in core-hours, log scale.
Plotted: the full-set chain only (solid line, filled markers), with all runs warm-started. Its cost is 40 L8 runs and 16 L9 runs (10 estimated warm, 6 measured warm) as the base, 3,170 core-h, then 667 core-h per L10 run. Each point is one number of L10 runs, n = 3, 4, 6, 8, 10, 12.
Also shown: H1 (2,043 core-h base, then 703 per L10 run) and L10 GP only (726 core-h per run).

Variant table (cheap chain subsets, not plotted). Errors are relative RMS, held-out. Errors at n = 3, 4, 6, 8, 10, 12 L10 runs.
| variant | n_L8 | n_L9 | base core-h | tau_95 error (%) | EDR error (%) |
|---|---|---|---|---|---|
| full set | 40 | 16 | 3,170 | 5.4, 5.8, 5.4, 5.2, 4.3, 2.3 | 8.7, 6.8, 8.9, 8.1, 6.6, 5.3 |
| warm9 | 40 | 6 | 2,488 | 9.0, 9.7, 6.4, 5.7, 5.1, 2.7 | 23.9, 43.5, 33.1, 23.9, 17.8, 11.9 |
| warm9, 20 L8 | 20 | 6 | 1,455 | 6.9, 6.9, 5.1, 5.1, 4.7, 3.4 | 48.9, 504, 58.9, 42.3, 28.4, 15.5 |
| warm9, 10 L8 | 10 | 6 | 918 | 5.6, 5.2, 4.9, 5.1, 4.0, 3.3 | 53.6, 1,080, 64.6, 43.5, 31.3, 18.9 |
Each cheap variant adds 703 core-h per L10 run.

Crossings, full-set chain only (core-h):
- tau_95 at 10 %: 5,170 core-h (n = 3, error 5.4 %).
- tau_95 at 20 %: 5,170 core-h (n = 3).
- tau_95 at 5 %: 9,837 core-h (n = 10, error 4.3 %). Not reached at n = 6 or 8.
- EDR at 10 %: 5,170 core-h (n = 3, error 8.7 %).
- EDR at 20 %: 5,170 core-h (n = 3).
- EDR at 5 %: not reached by n = 12 (5.3 % at n = 12).

Kim's code, measured on Oscar (job 8363941, N=1024, 64 ranks, steady t=0.1–0.4): ≈102 core-h per rocking cycle, against 26.3 for ours (v_l10_rpm32.5_th7); at N=256, 0.98 against 0.22.
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
