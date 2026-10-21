# Talking points, 2026-10-22 (order of presentation)

1. Definitions first (00, 11): how kLa and mixing time are computed; Kim's own code gives the same numbers and the same jumps. Not a postprocessing bug.
2. The two 5x3 grids (01, 02): what the models predict at L10 and what data are missing (L9 kLa only at 7 deg; L10 mixing time at 3 runs; the 6 warm L9 runs at 4 and 9 deg are running).
3. The headline (04): one converged flow period repeated predicts the normal run: kLa(0.10) to 8 %, the 50 % mixing time to 10 %, at 1/25 of the cost; the scatter is a property of each condition's flow, not noise. mfbml itself is not better than a constant here.
4. The C* question (03): the lower threshold carries the information; a constant reaches the run-to-run floor from 4 runs; stopping at C* = 0.25 saves about 4x.
5. Leads (05, 06, 07, 09): kLa ~ dissipation^0.6 at L9/L10 (blind test queued); the 22.5 rpm window seen only at fine levels; warm starting saves 7-13x in cycles; L9 + 3-4 L10 runs gives ~5 % for shear and dissipation.
6. Experiment only, if asked (08): same flow, finer scalar mesh changes the values by about +30 % per level. We make no claim about which value is physical.

Say plainly: negative results: mfbml adds no accuracy over a ratio fitted on a cheaper level; L8 does not predict L10 kLa; the L8 -> L10 ratio varies with rpm and angle.
Still running: 6 warm L9 runs (4 and 9 deg), 3 protocol checks (7 deg), 7 resumed L10 runs (L10 mixing time, ~2 days), 5 L10 oxygen runs at 5.5 deg (blind test of the dissipation law).
