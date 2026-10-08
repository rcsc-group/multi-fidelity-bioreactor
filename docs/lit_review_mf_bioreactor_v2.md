# Literature review v2: multi-fidelity methods and bioreactor flows

Date: 2026-10-08. The reviewer reports evidence. It takes no decisions.
This file extends `lit_review_mf_bioreactor.md` (v1). Raw data: `/oscar/scratch/eaguerov/tmp/ss_raw/`.

Evidence labels. "Opened" = the reviewer read full text (arXiv tools). "Abstract only" = only title and abstract seen. "Unverified" = not checked in a source.

## 1. Task 1: prior art with the Semantic Scholar Graph API

### 1.1 What was run

- Scripts: `run_queries.py`, `run_retry_and_citations.py`, `screen.py`, `analyze.py`, `show.py` (all in `ss_raw/`).
- Free tier, no key. 3 s between calls. On HTTP 429: wait 30 s, retry up to 5 times. Then 3 extra retry rounds over the failed queries.
- The free tier was saturated. 429 was very common. 12 of 35 query files still have no data after all retries. A key is needed to finish them. No key was sought.
- Raw JSON: `ss_raw/search_*.json` (100 hits per query, with abstracts) and `ss_raw/cites_*.json`.
- Note: `t.json` is a 5-hit probe of `multi-fidelity bioreactor` (HTTP 200, total 33,899). It confirms the API worked at that moment. It is not a full result.
- Screen rule: a hit counts as "relevant" if title+abstract match a multi-fidelity regex (multi-fidelity, co-kriging, variable fidelity, multi-level, bi-fidelity, low-fidelity, transfer learning) AND a bioreactor regex (bioreactor, cell culture, cultivated meat, bioprocess, CHO, shaken, rocking, single-use, etc.). "+flow" adds a flow regex (CFD, mixing, kLa, shear, flow, etc.).
- The regex is loose. Every "relevant" hit was then read by title and judged by hand. Most are false hits (see 1.3).

### 1.2 Query log

total = S2 estimate of matches. MF&bio = regex hits among the 100 returned. FAIL = HTTP 429 after all retries, no data.

| # | Query | Total | Returned | MF&bio | +flow |
|---|---|---|---|---|---|
| 1 | multi-fidelity wave bioreactor | 1055 | 100 | 2 | 2 |
| 2 | multi-fidelity rocking bioreactor | 149 | 100 | 2 | 2 |
| 3 | multi-fidelity stirred tank bioreactor | 105 | 100 | 0 | 0 |
| 4 | multi-fidelity cell culture CFD | 203 | 100 | 0 | 0 |
| 5 | multi-fidelity cultivated meat | 87 | 87 | 0 | 0 |
| 6 | multi-fidelity shear stress cells | 1263 | 100 | 1 | 1 |
| 7 | multifidelity shaken bioreactor | 388 | 100 | 0 | 0 |
| 8 | co-kriging bioprocess CFD | 41 | 41 | 0 | 0 |
| 9 | co-kriging mixing time | 10476 | 100 | 0 | 0 |
| 10 | co-kriging shear stress cells | 1105 | 100 | 0 | 0 |
| 11 | variable fidelity wave bioreactor | 477 | 100 | 1 | 1 |
| 12 | variable fidelity rocking bioreactor | 71 | 71 | 1 | 1 |
| 13 | multi-level cell culture CFD | 3894 | 100 | 5 | 1 |
| 14 | multi-level cultivated meat | 3353 | 100 | 1 | 1 |
| 15 | multi-level shaken bioreactor | 502 | 100 | 10 | 4 |
| 16 | low-fidelity high-fidelity bioprocess CFD | 35777 | 100 | 1 | 1 |
| 17 | low-fidelity high-fidelity kLa | 32535 | 100 | 0 | 0 |
| 18 | low-fidelity high-fidelity mixing time | 88018 | 100 | 0 | 0 |
| 19 | transfer learning surrogate bioreactor | 1620 | 100 | 0 | 0 |
| 20 | transfer learning surrogate wave bioreactor | 48 | 48 | 0 | 0 |
| 21 | transfer learning surrogate rocking bioreactor | 11 | 11 | 0 | 0 |
| 22 | transfer learning surrogate stirred tank bioreactor | 840 | 100 | 0 | 0 |
| 23 | Bayesian optimization bioreactor CFD | 404 | 100 | 1 | 1 |
| F1 | multi-fidelity bioreactor | FAIL (probe: 33,899) | | | |
| F2 | multi-fidelity kLa | FAIL | | | |
| F3 | multi-fidelity sloshing | FAIL | | | |
| F4 | multi-fidelity free surface two-phase flow Bayesian optimization | FAIL | | | |
| F5 | multifidelity cell culture CFD | FAIL | | | |
| F6 | multifidelity cultivated meat | FAIL | | | |
| F7 | multifidelity single-use bioreactor | FAIL | | | |
| F8 | multi-level single-use bioreactor | FAIL | | | |
| F9 | co-kriging kLa | FAIL | | | |
| F10 | low-fidelity high-fidelity shear stress cells | FAIL | | | |
| F11 | variable fidelity bioreactor | FAIL | | | |
| F12 | variable fidelity stirred tank bioreactor | FAIL | | | |

Result: 23 queries with data, 12 failed. 1,865 unique papers were screened. Spelling variants (multi-fidelity vs multifidelity) were mostly covered by only one of the two spellings, because S2 search is fuzzy. The failed rows leave these gaps: single-use bioreactor, kLa for the multi-fidelity spelling, shear stress for the low-fidelity spelling, sloshing.

### 1.3 Candidates opened or screened

"Fidelities" = what the low and high levels are. "Bioreactor?" = is the system a bioreactor.

| Paper | Fidelities combined | Bioreactor? | Evidence level |
|---|---|---|---|
| Eskandari et al., arXiv 2311.05776 (2023) | 3D CFD (RANS, Eulerian multiphase) vs ideal-mixing 0D model. Discrete misoKG. | Yes: industrial syngas fermentation (gas-liquid, not cell-culture, not rocking) | Opened in full this time. The text shows a BoTorch Hartmann6D test study and a plan ("next step ... MFBO campaign"). It shows no result table for the bioreactor simulator. Fidelities differ in physics model, not in grid. |
| Sun et al., arXiv 2211.14493 (2022) | Small-scale vs large-scale bioreactor process data; cell lines | Yes, but process data, no CFD | Opened in v1 (intro). |
| Martens et al., arXiv 2508.10970 (2025) | Experimental scales (plate, mini, pilot). Benchmark is an ODE CHO model. | Yes, process level, no CFD | Opened in v1 (intro, methods start). |
| Triantafyllou et al., arXiv 2609.17440 (2026) | Flowsheet simulator fidelities | Bioprocess flowsheet, no flow | Abstract only (v1). |
| Savage et al., arXiv 2305.00710 and 2210.17213 | Mesh-level CFD fidelities (cells count) | No: helical-tube chemical reactor, tracer mixing | 2305.00710 opened in v1. 2210.17213 abstract only. |
| Savage et al., arXiv 2308.08841 | CFD + multi-fidelity BO, coiled reactors | No: chemical flow reactor | Abstract only. |
| van Halder et al., arXiv 2004.13128 | Grid levels of a 3D PIC/FLIP sloshing tank (32x16x16 doubled per level). Uncertain inputs: tank motion amplitude and wavelength. QoI: particle counts per cell. Section "Surrogate Modelling for 3D Fluid Sloshing". | No: single-phase free-surface sloshing tank (computer-graphics solver) | Passages of the sloshing section read this time. |
| van Halder et al., arXiv 2106.03491 | Low and high grid resolution PIC/FLIP, intrusive CNN | No (sloshing tank) | Opened in v1 (intro and method start). |
| Setinek et al., arXiv 2511.01830 | Two RANS fidelities: wall-resolved (y+<1) vs wall-function (y+ 30-300) | No: airfoils | Opened in full. |
| Out-of-the-loop MF-BO, arXiv 2608.04113 | Generic MF-BO; chemistry, hyperparameters | No | Abstract only. Cites Martens 2508.10970. |
| Human-in-the-loop BO, arXiv 2606.19230 | Single-fidelity GP BO, CHO simulator | Yes, process level. Not multi-fidelity | Abstract only. |
| Lafzi and Dabiri, arXiv 2103.11303 | Multi-fidelity GP, droplet in oscillatory microchannel flow | No | Title/abstract snippet only. |
| "Bayesian Optimization in Bioprocess Engineering - Where Do We Stand Today?" (2025, cites 2311.05776) | Review | Bioprocess | Not opened. Title only. |
| "Machine learning and physics-driven modelling of multiphase systems" (Int. J. Multiphase Flow 2024, cites 2305.00710) | Review | No | Not opened. Title only. |

False positives from the regex (not relevant): "multi-level" metabolic engineering (K. phaffii, S. albus), "multi-level perspective" on cultured meat politics, multi-level climate governance, multilevel laminar-shear perfusion bioreactor design (a device with several levels, not a model), microfluidic multi-level molds.

### 1.4 Citation chasing

| Seed | Citing papers returned by S2 | Multi-fidelity + bioreactor + flow among them |
|---|---|---|
| 2311.05776 (Eskandari) | 1 | none (one review of BO in bioprocess, not opened) |
| 2305.00710 (Savage) | 20 | 2511.23140 (burner), 2308.08841 (flow reactor), 2026 comp. chem. eng. (burner, DOI 10.1016/j.compchemeng.2026.109787): none is a bioreactor |
| 2210.17213 (Savage DGP) | 5 | none new |
| 2504.05421 (Kim, rocking bioreactor CFD) | 3 | none: aquaculture tank sloshing, cultured-meat supply-chain review, interfacial particles |
| 2211.14493 (Sun) | 6 | none with CFD (knowledge transfer across cell lines/scales; kinetic models) |
| 2508.10970 (Martens) | 3 | 2608.04113, 2606.19230, one hybrid-model paper: none uses flow CFD |
| 2004.13128 (van Halder) | 2 | 1 multifidelity encoder/decoder CFD paper (AIAA 2022-0803): not a bioreactor |
| 2507.11640 (PINN stirred tank) | 0 | none |

S2 citation lists are often shorter than real ones. Treat as lower bounds.

### 1.5 Verdict

Statement X (precise): "A multi-fidelity surrogate or Bayesian-optimisation study where the fidelities are CFD grid/resolution levels of a two-phase Navier-Stokes solver, applied to a rocking (or wave, or shaken, or single-use) cell-culture bioreactor, with oxygen transfer (kLa) and cell shear stress as objectives."

- No prior art found for X in S2 (23 queries with data, 1,865 papers, 8 seeds of citation chasing) plus the arXiv searches of v1.
- Prior art exists for these neighbours:
  - Multi-fidelity BO on a bioreactor CFD vs 0D model: Eskandari 2311.05776 (syngas, gas-liquid, industrial scale). Different physics fidelities, not grid levels. No bioreactor result shown in the text.
  - Multi-fidelity ML on bioreactor process data (no CFD): 2211.14493, 2508.10970, 2609.17440.
  - Multi-fidelity BO with mesh-level CFD fidelities on a mixing reactor: Savage 2305.00710, 2210.17213, 2308.08841 (chemical reactor, not a bioreactor).
  - Multi-fidelity (grid-level) neural surrogates on a free-surface sloshing tank with prescribed tank motion: van Halder 2004.13128 (single-phase PIC/FLIP, not a bioreactor, no BO).
  - Single-fidelity CFD of a rocking bioreactor: Kim et al. 2504.05421. Single-fidelity BO on shaken bioreactors and photobioreactors (v1, abstract only).
- Limits: 12 queries failed (429). Not searched: Web of Science, Scopus, Google Scholar, publisher full-text, AIChE/Mixing conference proceedings, theses. S2 coverage of non-arXiv engineering journals is incomplete. Absence in these sources does not prove absence.

## 2. Task 2: canonical flows for a multi-fidelity demonstration

Criteria in order: (1) wow for an APS DFD audience, (2) zero or low compute, ideally an open multi-fidelity dataset, (3) published reference values, (4) free-surface or mixing link. Basilisk not required.

URL status was checked with `curl` on 2026-10-08 (script `ss_raw/check_urls.sh`).

### 2.1 Ranking (the reviewer's reading of the evidence)

**Rank 1. SPH dam-break against a wall, four resolutions (Lipari and Vuik, TU Delft / 4TU).**
- Solver and fidelities: DualSPHysics 5.0, weakly compressible SPH. H/dx = 800, 1600, 3200, 6400. 1M, 5M, 20M, 82M fluid particles. Re_eff up to 256,000. Same flow, four levels. Cost is already paid.
- Data: collection https://doi.org/10.4121/c.5353691 (resolves, HTTP 200, redirects to data.4tu.nl/collections/fd47d945-5e9d-478b-8322-69b7f14b6e4f). General page https://data.4tu.nl/datasets/f7659339-633e-4dbc-a6df-3afaaee7af57/2 resolves. CC BY 4.0. About 735 GB, 1,650 files (from a search snippet). The page states file downloads were "currently unavailable due to storage maintenance" when read. Previews only. This is the main risk.
- Reference values: the data are the reference. The `Commentary.pdf` in the collection describes the QoIs (not read).
- Multi-fidelity ML papers on it: none found (searches in S2 and arXiv). Unverified beyond that.
- Why impressive: free surface, violent impact, four resolutions up to 82M particles, zero compute for us. Weakness: SPH, not a grid solver. Fidelity is particle count. The data are large.

**Rank 2. Rising bubble benchmark of Hysing et al. (2009), 2D, two cases.**
- Solvers and fidelities: three codes (TU Dortmund TP2D, EPFL FreeLIFE, Magdeburg MooNMD) at up to 5 grid-refinement levels (files `c#g#l#.txt`: case, group, level). Listed levels: c1g1 l4-l7, c1g2 l1-l3, c1g3 l1-l4, c2g1 l4-l8, c2g2 l1-l3, c2g3 l2-l4.
- Data: ASCII time series (time, centre of mass, circularity, rise velocity). All-in-one zip, 1.7 MB, HTTP 200: http://www.featflow.de/media/bubble/data_bench_quantities.zip . Index page: http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble/bubble_reference.html (HTTPS handshake fails from this host; plain HTTP works).
- Reference values: yes. Test 1: Re=35, Eo=10, density and viscosity ratio 10. Test 2: Re=35, Eo=125, ratios 1000 and 100. Paper: Hysing et al., Int. J. Numer. Meth. Fluids 60, 1259 (2009), DOI 10.1002/fld.1934. Verified on the Featflow page.
- Multi-fidelity ML papers on it: none found (arXiv search "Hysing rising bubble ... multi-fidelity" returned only generic MF papers).
- Basilisk: `examples/bubble.c` exists (https://basilisk.fr/src/examples/bubble.c, HTTP 200). The reviewer did not confirm it uses the Hysing setup (the text found uses unit gravity). Do not assume.
- Why impressive: two-phase interface flow with surface tension, grid levels and published QoIs, free data, tiny. It looks like bioreactor physics (bubbles, interface). Weakness: classical, 2D, low Re.

**Rank 3. SPHERIC Test 10, sloshing wave impact (experiment).**
- Data: roll-motion tank, pressure time histories for the first impacts, water and sunflower oil, repeatability records, videos. ZIP of 236 MB (from a search result). Page resolves: https://www.spheric-sph.org/tests/test-10 (HTTP 200).
- Fidelities: the experiment is the high fidelity. Low fidelities would be our own cheap runs (2D, coarse grids). So this one needs compute.
- Multi-fidelity ML on it: none found.
- Why impressive: closest physics to a rocking bioreactor (periodic roll of a partly filled tank, free surface). Experimental ground truth. Weakness: impact pressure is a spiky QoI, not smooth.

**Rank 4. Airfoil RANS two-fidelity data (Setinek, Galletti, Brandstetter, arXiv 2511.01830, built on AirfRANS).**
- Fidelities: high = wall-resolved (first cell about 2 micrometre, y+<1, 13.4 core-hours, 180K nodes). Low = wall-function (first cell 1,200 micrometre, y+ 30-300, 4.8 core-hours, 96K nodes). 611 matched pairs. OpenFOAM simpleFoam, k-omega SST, Re 2e6 to 6e6.
- Data: the paper (opened in full) does not state a release link. Dataset availability: unverified. AirfRANS itself (single fidelity) is open: https://github.com/Extrality/AirfRANS and https://airfrans.readthedocs.io (both resolve).
- Finding relevant to our shear objective (from the opened text): pressure and velocity show positive transfer from low to high fidelity. Wall shear stress shows none (nMAE of LF vs HF: 0.405 and 0.796 for the two WSS components, 0.043 for surface pressure). The authors say wall shear is sensitive to boundary-layer resolution.
- Why impressive: current ML-for-CFD topic, compute-budget scaling laws. Weakness: not free-surface, not mixing, and a workshop paper.

**Rank 5. RANS vs DNS/LES curated dataset (McConkey, Yee, Lien, arXiv 2103.11515).**
- Fidelities: four RANS models (k-eps, k-eps-phi-f, k-omega, k-omega SST) with matching DNS/LES for periodic hills, square duct, parametric bumps, converging-diverging channel, curved backward-facing step. 29 cases per model, 895,640 points.
- Data: https://doi.org/10.34740/kaggle/dsv/2044393 (resolves to Kaggle, HTTP 200). Reference values: the DNS/LES data themselves (Xiao et al. periodic hills Re=5600, etc., as listed in the paper).
- Multi-fidelity ML on it: the dataset is for closure modelling. A multi-fidelity use was not found.
- Why impressive: DNS-quality truth with several RANS levels, parametric sweeps. Weakness: no free surface, single phase.

**Rank 6. Turbulence databases with several resolutions (JHTDB, BLASTNet 2.0).**
- JHTDB https://turbulence.pha.jhu.edu/ redirects to https://turbulence.idies.jhu.edu/home (HTTP 200): DNS of isotropic turbulence and channel flow at several Reynolds numbers, queryable by web service. BLASTNet 2.0 https://blastnet.github.io (HTTP 200): 2.2 TB, 744 samples from 34 DNS (compressible, reacting and non-reacting), used for super-resolution (NeurIPS 2023 D&B, arXiv 2309.13457).
- These give high-fidelity DNS and the user derives coarse versions by filtering. That is super-resolution, not independent LF solvers. No MF-BO paper found.
- Why impressive: huge, canonical. Weakness: filtered data are not true low-fidelity simulations. No free surface.

### 2.2 Other items seen (not ranked)

- Cylinder at Re=200 and 3900, Bayesian optimisation with gradient and/or multi-fidelity information, drag reduction (thesis "Variations on Bayesian optimization applied to numerical flow simulations", Univ. Cadiz, https://rodin.uca.es/handle/10498/29454). Abstract-level evidence only (search snippet). It shows that cylinder MF-BO already exists.
- "Multi-fidelity Rayleigh-Benard" and "multi-fidelity cavity": no hits in arXiv searches in this session (consistent with v1).
- Other open multi-fidelity CFD datasets named in web search snippets (not opened): RANS-LES jet database (92 RANS, 18 LES), double-delta wing (VLM and RANS, 2448 snapshots), AhmedML (500 geometries, hybrid RANS-LES; opened abstract page: single fidelity, CC-BY-SA, OpenFOAM case included).
- Basilisk: v1 found `test/capwave.c` (capillary wave with Prosperetti analytic reference) as the nearest free-surface test.

### 2.3 What was not done

- Dataset downloads were not tested. Only landing pages and a 1.7 MB zip header were checked.
- The 4TU `Commentary.pdf` was not read.
- Multi-fidelity ML on each dataset was searched with arXiv and S2 queries only (not Google Scholar).
