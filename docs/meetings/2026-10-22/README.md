# Meeting 2026-10-22

The figures are copies from 21 October. The scripts and the rules for each test are in `diary.md` (entries from 2026-10-19 to 2026-10-21).

Colours: L8 blue, L9 green, L10 pink. Marker shape shows the quantity.

You asked for four things: (1) a 5x3 chain L8, L9, L10. (2) a 5x3 two-level L9, L10. (3) the C* extrapolation. (4) a model where the cheap level is one repeated flow period. Cards 05 to 11 are extra results.

## 00 How we compute kLa and mixing time
File: `00_definitions.png`. One L8 run (25 rpm, 7 deg).
- Left: oxygen in the liquid, C*(t). Markers at C* = 0.10, 0.25, 0.50.
- Right: tracer mixing, chi(t). Markers at chi = 0.50, 0.75, 0.95.
- kLa at C* = x is the slope of ln(1 - C*) over one rocking period around the first time C* reaches x.
- The oxygen curve is not one exponential. So kLa depends on x. The median ratio kLa(0.50) / kLa(0.10) is 0.51.
- Code: `scripts/postprocess.py`.

## 01 Chain L8, L9, L10 (5 rows, 3 columns)
File: `01_chain_L8_L9_L10_5x3.png`.
- Rows: 95th-percentile shear, mean shear, dissipation (EDR), kLa at C* = 0.10, mixing time at chi = 0.50. Columns: angle 4, 7, 9 deg.
- Dashed lines: the fits of L8 and L9. Grey line and band: L10 prediction with the 95 % range.
- The model uses the log of each quantity. It trains on all L10 runs of the row. The L9 uncertainty goes into the L10 band.
- All curves span 14.5 to 38.5 rpm. The L10 runs span 17.5 to 32.5 rpm. So you see the extrapolation.
- Shear and dissipation: leave-one-out error 4.6 %, 3.0 %, 4.1 %. Test on 5 new runs at 5.5 deg: 2.0 %, 3.1 %, 7.3 %.
- Training set: all 13 L10 runs (7 at 7 deg, 3 at 4 deg, 3 at 9 deg) train each shear row, the dissipation row and the kLa row. Each panel shows only the runs of its own angle.
- The dissipation, shear, kLa and mixing-time axes are logarithmic.
- Gaps: L9 kLa and mixing time now exist at 4, 7 and 9 deg (three L9 runs each at 4 and 9 deg). L10 mixing time has 3 runs: no prediction.

## 02 Two-level model L9 to L10 (5 rows, 3 columns)
File: `02_two_level_L9_L10_5x3.png`.
- Same layout. The cheap level is a fit of the L9 runs. The L10 runs correct it. L8 is not used.
- A ratio fitted on L9 is as good as this model. The model adds an uncertainty band.

## 03 C* extrapolation
File: `03_threshold_extrapolation.png`.
- Idea of your PI: find kLa at a higher C* from a lower C*.
- Cheap level: kLa at the lower C*. A run can stop early. Expensive level: kLa at the higher C*.
- Plot: error on unseen runs against the number of expensive runs. Data: 79 L8 runs.
- Methods: mfbml. GPR with a linear prior. GPR with the lower-C* value as prior. The lower value times one constant.
- Result: the lower C* carries the information. One constant gives an error of 0.20 to 0.27 from 4 runs. mfbml is not better.
- Practical: a run that stops at C* = 0.25 saves a factor 4 in time. It predicts kLa(0.50) to about 20 %.
- Pair 0.25 to 0.75 is hard: error 0.4.

## 04 One repeated flow period as the cheap level
File: `04_repeated_period_mfbml.png`.
- Cheap level: record one converged flow period. Repeat it 22 cycles. Solve oxygen and tracer on it. Cost: 1.2 core-hours.
- Expensive level: the normal L8 run. Cost: about 32 core-hours.
- Data: 39 conditions, angles 4, 5.5, 7, 9 deg.
- Top: cheap against expensive value. Bottom: error against the number of normal runs.
- The two agree. Correlation of the logs: 0.98 (kLa at 0.10), 0.91 (kLa at 0.25), 0.99 (mixing time).
- The scatter around the trend agrees too. So the scatter belongs to the flow of each condition. It is not noise.
- With 4 normal runs and one constant, the error is 0.08, 0.19, 0.11. A GPR on the normal runs alone gives 0.57, 0.71, 0.36.
- mfbml is not the best method here. A constant is as good or better.
- Limit: the cheap level starts from the converged state of the normal run. A new condition needs 12 warm-start cycles first.
- Check at 32.5 rpm: the cheap level differs from the cold run by 2.7 % and 6.7 %. A coupled control with the same restart agrees to 0.55 %. The gap comes from the release time.
- Not tested at L9 and L10.

## 05 A law for kLa from the dissipation rate
File: `05_kla_law_dissipation.png`.
- Left: kLa against EDR at L9 and L10. Middle: mixing time against EDR at L8 and L9. Right: error of each method for kLa at L10.
- kLa follows EDR to the power 0.6 at L10 (13 runs, correlation 0.97) and 0.5 at L9. L8 shows no relation.
- The interface area is the same at all levels. So the law is for the transfer velocity.
- Small-eddy theory gives 0.25. We do not explain our value.
- Blind test: we froze the predictions for five runs at 5.5 deg. The runs are queued.

## 06 Warm start
File: `06_warm_start.png`.
- A: change of condition at one level. B: change of level. C: is the training flow settled? D: cycles before the release.
- Cold start: 80 cycles. Start from L8: 12 cycles. Start from L9: 6 cycles.
- Resume at the same condition: exact at L8.
- A change of rpm or angle is the least reliable.

## 07 Cost against error, shear and dissipation
File: `07_cost_vs_error.png`.
- L9 at every design point plus 3 or 4 L10 runs gives about 5 % at about 4,000 core-hours.
- L10 runs alone give 13 % to 190 % at that cost.
- A ratio fitted on L9 is as accurate as the chain.

## 08 Experiment: same flow, finer mesh for the scalars
File: `08_fixed_flow_experiment.png`.
- One L8 flow. We solve oxygen and tracer on the L8, L9 and L10 mesh.
- Mixing time: 11.8 s, 19.0 s, 24.7 s. kLa at 0.25: 28.6, 38.7, 51.6 per hour.
- This is an experiment. We make no claim about the physical value.
- Limits: one flow. The finer-mesh flux is interpolated (about 5 % error).

## 09 Runs 0.2 rpm apart
File: `09_run_pairs.png`.
- At 30 rpm, L9 pairs differ by 2 % to 3 %. L8 pairs differ by 12 % to 37 %.
- At 22.5 rpm, both levels differ a lot.
- At 22.4 and 22.6 rpm, the L9 dissipation differs by 30 %. The late mixing rate differs by a factor 3.6.
- A narrow flow regime near 22.5 rpm exists at L9 and L10. L8 does not show it.

## 10 One repeated period, first test
File: `10_one_period_repeated.png`.
- L8, 32.5 rpm. kLa at 0.10 and mixing time agree to 0.6 % and 0.0 %. kLa at 0.25 differs by 15 % to 51 %.
- Card 04 replaces it.

## 11 Is it a bug?
File: `11_is_it_a_bug.png`.
- Left: the same run again. Difference: 0.01 % to 1 %.
- Middle: Kim's code on our fields. His default mixing time is 8 % (L8) and 3 % (L9) slower. He counts air-filled surface cells. Without them the difference is below 2 %. His kLa fits agree to 1.6 % at 0.25 and 0.50.
- Right: another estimator keeps the scatter.
- The jumps between runs 0.2 rpm apart are in the simulated flow. Kim's code sees them too.

## Questions for the PIs
1. Which kLa definition do the experiments use?
2. Is there any measured kLa of this bag?
3. Which mechanism gives an exponent near 0.5?
4. What convergence criterion would you accept for a scalar?

## Feedback after the meeting
(empty)
