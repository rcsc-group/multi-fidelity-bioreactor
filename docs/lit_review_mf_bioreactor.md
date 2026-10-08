# Literature review: multi-fidelity methods for bioreactor flows

Date of search: 2026-10-08. The reviewer reports evidence only. It takes no decisions.

Method note. "Opened" means the reviewer read the text of the paper (arXiv full text). "Abstract only" means the reviewer did not read the methods. Claims from memory are marked "unverified".

## 1. APS DFD 2026 abstract (Aguero Vera)

Result: NOT FOUND. No abstract text is reported here.

What was searched:

- Web search: `APS DFD 2026 "Aguero Vera" Brown University abstract`. Returned only DFD 2024 and 2025 pages and unrelated items.
- Web search: `meetings.aps.org DFD26 "Elvis Alexander Aguero Vera"`. Returned nothing relevant.
- Web search: `"DFD26" Bulletin ... rocking bioreactor multi-fidelity`. Returned older DFD talks (2021, 2022, 2023) on rocking bioreactors. None is by the user.
- Direct fetch of `meetings.aps.org/Meeting/DFD26/Content/Search?q=Aguero`. Redirects (302) to `meetings-archive.aps.org/dfd/2026`. That URL returns 404. `archive.aps.org/dfd/2026` returns 404. `meetings.aps.org/Meeting/DFD26/Session/A01` redirects to `meetings-archive.aps.org/dfd/2026/a01`, which returns 404.

Interpretation (not proven): the tool cannot reach the DFD26 program pages. The abstract may exist but is not reachable from here. The user should paste the abstract text or the session URL.

## 2. Prior art: multi-fidelity methods and bioreactors

### 2.1 Query log

Tools: WebSearch (a general web search, not Google Scholar itself), arXiv search tool, Semantic Scholar API, Crossref API.
Semantic Scholar API: HTTP 429 (rate limit), no data. Crossref: 1 query only (`Aguero Vera sloshing bioreactor`), no relevant result. Google Scholar: not accessible. So the log below is arXiv plus general web. Publisher-side journals are covered only by what the web search surfaced.

| # | Query (abridged) | Where | Relevant hits | Hits |
|---|---|---|---|---|
| 1 | multi-fidelity surrogate bioreactor CFD | Web | 2 | 2311.05776 (via snippet); Manchester/RSC chapter (single-fidelity GP) |
| 2 | co-kriging stirred tank mixing time multi-fidelity | Web | 0 | generic co-kriging papers only (1210.0686, 1112.5389) |
| 3 | (multi-fidelity OR multifidelity OR co-kriging) AND (bioreactor OR stirred tank OR cell culture) | arXiv | 3 | 2311.05776; 2211.14493; 2507.11640 (single-fidelity PINN) |
| 4 | abs "multi-fidelity" AND "mixing time" | arXiv | 0 | none |
| 5 | rocking/wave bioreactor + surrogate/optimization/Bayesian | arXiv | 0 | only 2504.05421 (single-fidelity CFD of rocking bioreactor) |
| 6 | multifidelity Bayesian optimization bioprocess CFD kLa shear stress | Web | 2 | 2311.05776; 2508.10970 |
| 7 | variable-fidelity surrogate wave / rocking / single-use bioreactor | Web | 0 | commercial pages only |
| 8 | multi-level Monte Carlo bioreactor CFD uncertainty stirred tank | Web | 0 | MLMC for RANS (1811.00872, not a bioreactor); CMA Monte Carlo (not MLMC) |
| 9 | multi-fidelity Gaussian process stirred tank mixing CFD surrogate chem. eng. | Web | 2 | Savage 2305.00710; hybrid-fidelity CSTR surrogate (ScienceDirect S1876107026000623, abstract only, not opened) |
| 10 | multi-fidelity bubble column / stirred tank NN kLa / dissipation | Web | 0 | single-fidelity ML for kLa only |
| 11 | GP / Bayesian optimization CFD bioreactor impeller shear stress | Web | 0 MF (3 single-fidelity) | see 2.3 |
| 12 | hierarchical surrogate shaken flask / orbital shaken CFD multi-fidelity | Web | 0 | generic MF surveys only |
| 13 | (multifidelity OR multi-fidelity OR variable-fidelity) AND (mixer OR mixing OR agitated OR bioprocess) | arXiv | 5 | 2508.10970; 2211.14493; 2609.17440; 2305.00710; 2311.05776 |
| 14 | (Gaussian process OR surrogate OR Bayesian optimization) AND bioreactor AND CFD | arXiv | 2 | 2311.05776; 2507.11640 |
| 15 | Bayesian optimization hollow orbitally shaken bioreactors | Web | 1 (single-fidelity) | J. Chem. Eng. Japan 2023 (DOAJ page returned 403; abstract seen in search snippet only) |
| 16 | Surrogate and multiscale modelling for (bio)reactor scale-up (chapter) | Web | 1 (single-fidelity) | Cho 2023 RSC chapter (abstract only) |
| 17 | multi-fidelity OR multifidelity photobioreactor / wave bioreactor / cultivated meat CFD surrogate | Web | 0 MF | same as 16, plus 2504.05421 |
| 18 | BO (shaken OR orbitally shaken OR rocking OR sloshing) bioreactor | arXiv | 0 | physics of orbital shaking only; 2504.05421 |
| 19 | Crossref: Aguero Vera sloshing bioreactor | Crossref | 0 | irrelevant |
| 20 | Semantic Scholar: multi-fidelity bioreactor | S2 API | n/a | HTTP 429, no data |

### 2.2 Candidates that could be prior art (opened and read)

1. **Eskandari, Puiman, Zeitler, arXiv 2311.05776** (2023). Opened; read problem and solution sections.
   - Fidelities: high = 3D CFD (RANS turbulence, Eulerian multiphase, 370,000 cells, about 5 days per run). Low = ideal-mixing 0D model (under 2 minutes).
   - Method: multi-fidelity Bayesian optimisation, discrete misoKG (Poloczek et al. 2017). Outputs: CO and H2 uptake rate q_i, gas conversion rate.
   - System: industrial-scale syngas fermentation bioreactor (gas-liquid). Not a cell-culture, rocking or shaken system. QoI is conversion rate, not shear, kLa-as-QoI or mixing time. The two fidelities are two different physics models, not grid levels.
   - Status: short workshop-style paper. It describes the approach and concerns. The part I read shows no results table. I did not read the last third (discussion, from character 12,000).
   - Verdict: multi-fidelity applied to bioreactor CFD, but only partly in scope.

2. **Sun et al., arXiv 2211.14493** (2022). Opened; read introduction and contributions only.
   - Fidelities: data from smaller-scale bioreactors (low) and larger-scale bioreactors (high). Second case: different cell lines as fidelities. Methods: Kennedy-O'Hagan linear MFGP and a nonlinear MFGP.
   - Data: real-world process data (the paper states "real-world datasets"). No CFD. I did not open the data section.
   - Verdict: multi-fidelity on bioreactors, but not on flow or CFD.

3. **Martens et al., arXiv 2508.10970** (2025). Opened; read introduction and the start of Methods.
   - Fidelities: experimental scales (microtiter plate, mini-bioreactor, pilot). Method: multi-fidelity batch BO with mixed variables. Benchmark: a custom ODE simulation of a CHO process (Monod kinetics).
   - No CFD. No flow QoI.
   - Verdict: multi-fidelity in bioprocess development, not flow.

4. **Triantafyllou et al., arXiv 2609.17440** (2026-09). Abstract only.
   - Fidelity-augmented GP BO on process flowsheet simulators (plasmid DNA bioprocess in SuperPro Designer; Aspen HYSYS). Flowsheet level, not flow.
   - Verdict: not flow CFD.

5. **Savage et al., arXiv 2305.00710** (2023). Opened; read abstract, introduction, and the fidelity section.
   - Fidelities: continuous mesh fidelities (cell count) of 3D CFD, plus others. Method: multi-fidelity BO (DARTS). System: helical-tube chemical reactor, pulsed flow, mixing QoI (tanks-in-series from tracer). Experimentally checked by 3D printing.
   - Not a bioreactor. It mentions "biologically-derived chemicals" only as motivation.
   - Verdict: closest chemical-engineering mixing precedent. It uses grid-based fidelities.

### 2.3 Near neighbours: single-fidelity surrogates / BO on bioreactor CFD

- Cho (Manchester), RSC 2023 chapter "Surrogate and multiscale modelling for (bio)reactor scale-up". GP surrogate of a multiscale CFD model inside BO; flat-plate photobioreactor with static mixers. Abstract only (page seen via search). Single fidelity.
- Hollow orbitally shaken bioreactors, Bayesian optimisation of diameter ratio, shaking speed and shaking radius; QoI includes mass transfer, suspension, shear stress. J. Chem. Eng. Japan, Dec 2023 (DOAJ). Abstract only; page returned 403. Single fidelity.
- Trávníková, von Lieres, Behr, arXiv 2507.11640. PINN surrogate for 2D stirred-tank flow field, Re 50 to 5000. Opened (abstract, introduction, methods start). It asks whether "lower fidelity or approximate data" can train the PINN, but it is not a multi-fidelity model in the Kennedy-O'Hagan sense.
- Kim, Harris, Cimpeanu, arXiv 2504.05421. CFD framework for a rocking bioreactor (mixing, oxygen transfer, shear). Single fidelity; no surrogate. Abstract only.
- Photobioreactor deep-learning surrogate (Imperial spiral record). Abstract seen via search only.

### 2.4 Near neighbours: multi-fidelity in other fluid/mixing settings (for context)

- Savage et al. 2305.00710 (above): chemical reactor mixing, mesh fidelities.
- Lima et al., arXiv 2511.23140: MF-BO of a burner, mesh element size as continuous fidelity. Abstract only. Not a bioreactor.
- Eaheart and Radaideh, arXiv 2603.14143: multi-fidelity GP and NNs on coarsened-mesh CFD of a gas reactor. Abstract only.
- ScienceDirect S1876107026000623: a physics-informed hybrid-fidelity surrogate of a stirred reactor; the search snippet says it was compared with a multi-fidelity GP. Not opened.

### 2.5 Conclusion (task 2)

**Borderline.** Not "no prior art".

- Prior art exists for multi-fidelity ML on bioreactors in general (2211.14493, 2508.10970, 2609.17440) and for multi-fidelity BO using CFD versus 0D models in a bioreactor (2311.05776).
- Not found in my searches: any multi-fidelity work on rocking/wave bioreactors, single-use bags, shaken flasks, or cell-culture bioreactors with CFD, with QoIs such as shear stress, kLa or mixing time, and grid-convergence or resolution-based fidelities.
- Limits of this claim: Google Scholar and Semantic Scholar were not usable; publisher paywalled pages were seen in snippets only; conference proceedings (AIChE, Mixing 18) were not searched in depth. Absence in these searches does not prove absence in the literature.

## 3. Canonical flow candidates for a cheap 2D multi-fidelity demonstration

Basilisk source checked locally at `/oscar/data/dharri15/eaguerov/basilisk/src` and the basilisk.fr index pages. Reference papers named from memory are marked "unverified".

### 3.1 Flow past a cylinder (Cd, St vs Re)

- (a) Basilisk: `examples/karman.c`, Re = 160, embedded boundary, centred NS solver, https://basilisk.fr/src/examples/karman.c . The file gives no Cd or St reference data. Also `test/cylinders.c` (Stokes flow, periodic cylinder array): https://basilisk.fr/src/test/cylinders.c .
- (b) Reference data: `cylinders.c` compares with Sangani and Acrivos 1982 (Int. J. Multiphase Flow 8, 193; verified in file). For Cd and St vs Re, the usual sources are Williamson 1989 and Henderson 1995 (unverified).
- (c) Multi-fidelity ML papers on this case: none found in my searches (arXiv queries 8 results; none used a cylinder MF benchmark). Not verified further.
- (d) Link to bioreactors: weak (no free surface). Applies to impeller/baffle wake only.

### 3.2 Lid-driven cavity

- (a) Basilisk: `test/lid.c`, Re = 1000, https://basilisk.fr/src/test/lid.c . Compares to Ghia et al. (file plots `yprof.ghia`; verified).
- (b) Ghia, Ghia, Shin 1982, J. Comput. Phys. 48 (author list from the file label "Ghia et al."; full citation unverified).
- (c) MF ML papers: none found with the cavity as the benchmark. Search returned only a PINN (2605.13892) and a viscoelastic cavity paper.
- (d) Link: weak (single phase, steady, closed box).

### 3.3 Sloshing tank / standing wave

- (a) Basilisk: no sloshing-tank example in `src/examples` or `src/test`. Closest: `test/capwave.c` (capillary standing wave, Popinet and Zaleski 1999; Navier-Stokes + VOF + surface tension; https://basilisk.fr/src/test/capwave.c ). `test/oscillation.c` is also present (not read). Basilisk is used for 2D and 3D sloshing and orbital sloshing at EPFL (page read: validation on simple 2D cases against Bouvard et al. 2017, Phys. Rev. Fluids 2).
- (b) Reference: `capwave.c` uses the Prosperetti analytic solution (`prosperetti.h` in the Basilisk tree; verified in file). Gravity-wave sloshing: Faltinsen and Timokha (unverified).
- (c) MF ML papers: van Halder, Sanderse, Koren, arXiv 2106.03491 (opened intro): intrusive deconvolutional NN trained on low- and high-fidelity data (grid resolution) for a 3D PIC/FLIP sloshing simulation in a rectangular tank. Jagtap, Mitsotakis, Karniadakis, arXiv 2202.02899: multi-fidelity data for inverse water-wave problems (abstract only; not a tank). A web-search snippet that credited a ship-hull MF paper to sloshing was checked and is wrong (the PDF is a hull-form paper, not sloshing).
- (d) Link to bioreactors: strong (free surface, rocking/orbital motion). Same physics class as the rocking bioreactor.

### 3.4 Rayleigh-Benard convection (Nu vs Ra)

- (a) Basilisk: no RB example in the local `src/examples` or `src/test` (only a mention in `ginzburg-landau.c`). Search snippets say Basilisk was used for RB in arXiv 1909.02270 and 2301.10833 (not opened).
- (b) Reference: Nu vs Ra scaling, e.g. Ahlers, Grossmann, Lohse, arXiv 0811.0471 (review; abstract only).
- (c) MF ML papers: none found with Nu(Ra) in my searches. Ravi et al. arXiv 2404.11965 (MF GP "for regression problems in physics") appeared; I did not open it, so I cannot say whether it uses RB.
- (d) Link: weak (heat transfer; no free surface). Long statistics make the run cost high.

### 3.5 Ranking (the reviewer's reading of the evidence)

1. Sloshing / standing wave. Reason: strongest link to rocking bioreactors; Basilisk has a related validated test (`capwave.c`, analytic reference); one MF-ML precedent exists (3D PIC/FLIP), so a 2D Basilisk grid-fidelity study is not a duplicate of it. Caveat: no ready Basilisk gravity-sloshing example; the user must build one.
2. Cylinder, Cd/St vs Re. Reason: Basilisk example exists (`karman.c`), cheap, QoI varies with Re, and published data exists (unverified). Weak bioreactor link. Caveat: St needs a long transient.
3. Lid-driven cavity. Reason: Basilisk test with reference data, very cheap. But a steady, smooth QoI gives a small and mostly constant-factor grid error, and it has no link to bioreactors.
4. Rayleigh-Benard. Reason: QoI Nu(Ra) is non-trivial, but no local Basilisk example, statistics are costly, and the link to bioreactors is weak.

## 4. Things not accessed

- APS DFD 2026 program pages (all routes returned 404 or redirect).
- Semantic Scholar API (429). Crossref used once.
- Google Scholar.
- Paywalled pages (ScienceDirect, RSC, J. Chem. Eng. Japan, MDPI returned 403): only abstracts or snippets seen.
- Methods sections of 2609.17440, 2511.23140, 2603.14143, 2404.11965, 2202.02899 and the data sections of 2211.14493, 2508.10970 were not read.
