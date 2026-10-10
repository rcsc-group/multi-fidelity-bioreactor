# Backlog

Ideas and loose ends not yet scheduled. Newest first within a section. When an
item starts, give it a diary entry; when it lands, delete it here.

## APS DFD 2026-11-25 (pivot 2026-10-08: vanilla mfbml, HF target)

- Canonical case: choose after the prior-art check (docs/lit_review_mf_bioreactor.md) and test A.
- L9 dtmix has no noise estimate: rep_l9_rpm25_rel82 was set up but never ran.
- The time-horizon fidelity (truncated chi curve) does not carry the late rate (test B); only
  revisit with a physical late-rate model, not a prior on r.

## Methods

- **kLa as a Floquet exponent of the oxygen operator ("hydrodynamics once,
  oxygen replayed").** Oxygen is passive, so on a settled periodic flow C*
  obeys a *linear, periodic* advection-diffusion problem and kLa is its
  leading decay rate -- exact, no closure. Record one settled period of
  (u, f); solve only the oxygen equation on it, looped, and get the decay rate
  by power iteration over periods instead of waiting 30-90 cycles for C* to
  hit 25%. Allows oxygen on a finer mesh than the flow (Sc ~ 500), the
  likely cause of L9's kLa scatter. Feeds Yi et al. cleanly: LF = replay kLa,
  HF = coupled L10 kLa, same QoI, so the linear transfer is justified.
  *First falsifier:* L6 serial, replay vs coupled kLa at the same level.
  *Before anything:* read how the interfacial oxygen BC is implemented.
- **Dense L8 grid for the MF demo.** Yi et al. assume abundant LF (200*d);
  `plot_mf_tau_max.py` has 9 L8 points. L8 is ~1.5 min/cycle at 8 ranks, so
  ~15-20 points over 15-40 rpm is about a day of compute and puts the LF in
  the regime the method is built for.
- **Heteroscedastic noise in KRR-LR-GPR.** L10 cycle-to-cycle spread of
  tau_max grows ~35x from 17.5 to 37.5 rpm; the GPR noise is homoscedastic,
  so the band is too wide at low rpm and the pooled value is set by 37.5.
- **KRR-LR-GPR band collapses at small n_HF.** The likelihood matches Yi eq 23
  (checked 2026-10-08), but predict() scales the variance by the plug-in
  sigma2 = r'K^-1 r / n; with n_HF = 3 and a linear transfer it goes to ~0
  (test A coverage 0-57%). Fix candidates: prior on the residual amplitude,
  or n_HF >= 5.

## Figures / runs

- Bare `return;` in VIDEOS-only events movies_output / movies_output_tau (src/BioReactor.c ~2512, ~2623): the CLAUDE.md hazard (can stop the run). Lean production builds are unaffected (no VIDEOS). Fix before any VIDEOS run.

- **Fig 8(b) tau tails ~13x Kim's** after the transient cleared. Test for
  grid-scale noise with a volume metric (EDR field spectrum / |grad u|^2
  distribution) at L8/L9/L10 before touching anything near the wall.
- `submit_fig8_hist.py` still restarts from zero forcing: set every `*_prev`
  to the current condition (as submit_fig11_l10_kla.py does). Cost of not
  doing it: ~15 cycles of L10 (~15 h) per re-record.
- Fig 9 L10 anchor (32.5 rpm): segment 2 running (job 6763888); ~4
  segments total. Add the chain to plot_fig9/plot_fig11 when it lands.
- MF tau_rpm: replace the 6-cycle L10 HF points at 35/37.5 with the ~30-40
  L10 cycles logged by fig11_l10_rpm35/37.5 (pre-release window).

- L6->L7->L8 kLa ladder at one bad L9 point (27.5 rpm or theta=4): does it
  converge toward Kim like 32.5 rpm, or stay anomalous?
- `make replicas` target wrapping the README's regenerate list.
- Unsubmitted: `submit_figset.py` (Figs 2-7), `submit_fig14.py` (14/15),
  `submit_figA18.py`. Parked by request: A.16(c).
- Reproducibility of production kLa: at L4, OpenMP round-off moves kLa_50 by 2x (diag_omp_determinism.py). Test whether L7 kLa changes with MPI rank count (8 vs 16). If it does, kLa error bars need a round-off/chaos term.
- t_end period rounding: FIXED in source (2026-09-30). Production scratch binaries still use the old rule; picked up at the next production rebuild.
- Fig 13 MF figures still use the old buggy 3-grid u_bar (0.24/0.30). Faithful E&H gives 0.5-1.5 (diag_eh_ubar.py). Decide: apply it, or get L11 grids.
- kLa: whole-period estimator adopted in Figs 11/12 (2026-10-01). Still open: re-check mixing time (dtmix) the same way; the CI restart test still uses the 5-point kLa.
