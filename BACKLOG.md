# Backlog

Ideas and loose ends not yet scheduled. Newest first within a section. When an
item starts, give it a diary entry; when it lands, delete it here.

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
- **Upstream likelihood check.** `mfbml` KRR-LR-GPR `_logLikelihood` uses
  `-0.5*n*sigma2` (not `log sigma2`) when the noise is a free hyperparameter.
  Upstream's own code, not project drift -- confirm against the paper's
  Appendix B before relying on free-noise fits.

## Figures / runs

- L6->L7->L8 kLa ladder at one bad L9 point (27.5 rpm or theta=4): does it
  converge toward Kim like 32.5 rpm, or stay anomalous?
- Figs 10/12: regenerate now that theta=2 has results (6/6 angles).
- `make replicas` target wrapping the README's regenerate list.
- Unsubmitted: `submit_figset.py` (Figs 2-7), `submit_fig14.py` (14/15),
  `submit_figA18.py`. Parked by request: A.16(c).
