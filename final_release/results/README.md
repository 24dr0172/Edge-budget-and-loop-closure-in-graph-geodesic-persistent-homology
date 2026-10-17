# Evidence for Paper 1

Closure means retaining every consecutive pair in the known clean cyclic order. It is specific to each graph order. Connectivity, closure, a nonempty diagram and low distance error are separate events. `R2_global` is radius/global-chord order; `OR_kNN` and `mutual_kNN` use minimum and maximum endpoint ranks. `R1_local` is CkNN and occurs only in the order-certificate example.

## Main results and their evidence

1. **Connected is not closed.** The circle parameter-choice theorem is supported by `clean_closure/event_table.csv` and `rank_landmarks.csv`. All 60 exact-count first-connected graphs per native order have empty loop diagrams. The first connected integer-k graph is empty in 38/60 OR and 40/60 mutual cases; it closes exactly when the two largest consecutive rank keys tie. These counts concern positive clean reference loops and the released tie convention.
2. **Closure switches the loop on under the conditions.** The pre-closure theorem and separate post-closure theorem are supported by last-pre-closure and closure rows for radius, OR and mutual in the 11-row event table. The per-order kNN checks are in `knn_event_checks.csv`; it retains only connectivity, first connected integer k, last pre-closure and closure. `own_closure_cells.csv` records each order's own closure; it is not a split at radius closure for all methods.
3. **After closure, extra edges do not improve the certified quantities.** `covered_monotonicity.csv` retains 180 native method–sample series with no saved covered-grid decreases in maximum distance error, RMS distance error or diagram error. `first_closed_errors.csv` reports all three errors for nine profile–order groups, without a redundant reference-selected Best column. The distance statement is restricted to covered budgets. For diagrams, the useful-diagram theorem guarantees the capped error; ordinary error needs its small-gap condition. Every retained clean sample satisfies g <= 0.089558 < 2/pi - 1/2. Figure 2 shows the saved clean grid only; Figure 1 is the conceptual closure schematic.
4. **A new order needs its short-hops-first check.** `order_checks/certificates.csv` has 60/60 passing radius, OR and mutual cases, and 59/60 passing CkNN cases. The failed CkNN certificate is preserved in `failure_witness.json`. Failure is inconclusive, and these counts are not a general CkNN theorem. CkNN has no error-ranking table or performance curve in this release. `order_checks/shuffled_control.csv` adds the placement control: 0/60 short-hops-first passes, six loops at first connectivity before closure, and 16/60 first-nonempty/first-closed grid matches. It contains diagnostics only.
5. **Under noise, inspect lifetime rather than exact emptiness alone.** `matched_noise/` contains the matched first-connected comparison for radius and OR, plus the before/at-closure comparison at delta=0.01. Short positive bars remove exact emptiness without recovering a comparably long loop. This is an empirical result, with no universal noise boundary or noisy-data selector.

## Saved-data extra-k diagnostic

`extra_k/samples.csv` has 120 paired rows from the existing rank audit: 20 seed streams x 3 profiles x 2 kNN rules, all n=512. The statistic is the individual difference k_close - k_conn. Its median is not computed by subtracting two marginal medians.

OR's observed maximum extra k is 3 in each profile; its medians are 1 (uniform), 1 (cosine), and 0.5 (trough). Mutual's medians are 0.5, 1, and 7, with maxima 4, 3, and 26. The sparse trough therefore illustrates how far mutual closure can lie beyond connectivity. These are descriptive saved-sample results, not a fixed recommendation to add a certain number of neighbors. The appendix ECDF displays the full distributions. The 120-row panel describes n=512; the rank-only sweep below measures sample-size effects.

## Rank-only sample-size sweep

`extra_k/candidate_audit/` adds 30,780 OR/mutual records from 15,390 profile–seed samples at n=128 through 1,000,000. There are 5,130 seed–size combinations shared across the three profiles. The exact local identity is checked against 18,980 full chord-rank values and 480 actual connectivity events. The rank-only study computes no persistence diagrams. All its gaps are below one third, so the existing circle theorem identifies ring closure with loop retention. Its README records the model calculation, targeted literature check, uncertainty and why the universal law and multi-loop table are excluded. The size sweep is distinct from the paired n=512 extraction.

## Supporting checks

`useful_diagram/saved_prefixes.csv` retains 3,556 closed graphs from 60 saved clean samples; `saved_transitions.csv` retains 3,376 transitions. The budgets are the distinct saved closed budgets plus each order's exact closure and the complete graph, not every intervening edge count. Each stored diagram has one positive bar with birth at the longest consecutive chord, within the recorded tolerance. The `.npz` checkpoints contain coordinates, reference diagrams and estimated diagrams. `useful_diagram/execution_summary.json` reports 3,259 overlapping saved-score comparisons agreeing to about 1e-16. This extraction checks the stored matching again; it does not rerun Ripser for those graphs.

`useful_diagram/rational_counterexample.json` retains the seven-point witness showing why the theorem caps the error. Radius M=13 to M=14 lowers ordinary error from about 0.04406 to the empty-diagram baseline 0.03006; the capped error stays fixed. The independent finite reduction is rerun by the verifier. The stored rational coordinates are evaluated with floating-point chord geometry.

`ellipse_checks/` retains the analytic locality and chord–arc constants, 60 reused coordinate sets, and 180 native order cases. The 539 distinct saved diagrams cover 540 connectivity, last-pre-closure and closure event checks; one connectivity event coincides with the last pre-closure event. The lower chord–arc factor is 1/pi, the local chord threshold is 0.18, and B < 1/(3pi) supplies the post-closure interval. The ellipse supports the smooth-curve scope of the closure results. Theorem 3 establishes the single-bar and capped-error conclusions for arc-supported graphs, including the local closed ellipse prefixes. This does not establish the circle numerical gap bound or ordinary-error monotonicity for ellipses. Both the locality argument and its analytic constants must remain in the paper's appendix.

## Sampling and units

The clean data use seeds 20310101–20310120 and the original uniform, cosine and quarter-circle trough densities. The three profiles share seed streams and the ellipse reuses the latent samples; these are not independent new replications. `clouds.csv` stores the reference longest lifetime divided by diameter in `oracle_longest_lifetime`. In the noise `cells.csv` that field is the raw reference lifetime; `longest_over_reference` is dimensionless. Checkpoints contain unnormalized diagram endpoints.

The matched-noise study uses 20 separate streams 20320101–20320120, n=512, and normal amplitudes 0, 0.0005, 0.001, 0.002, 0.005 and 0.01. The same latent points and signed displacement vector are reused across amplitudes. Its 240 rows are 20 matched streams, not 240 independent samples. At first connectivity the empty counts per method are 20, 9, 0, 0, 0 and 0. No row passes the half-reference-lifetime screen. At delta*n=5.12, median lifetime ratios are about 0.0253 for radius and 0.0291 for OR. `closure_comparison.csv`, `closure_summary.csv` and `closure_diagrams.npz` add the 80 before/at-closure observations at delta=0.01 from those same 20 seeds. Before closure, 20/20 graphs per order pass the half-lifetime screen; every closure diagram is nonempty. Gap crossing remains outside the package.

The label-free selection rule is conditional: the first nonempty diagram identifies closure only after the clean geometric and edge-order conditions are established. The controlled extra-k analysis itself uses the known cyclic sample order and is not a label-free tuning algorithm.

`matched_noise/coordinates.npz` stores all 20 seeds' latent arc coordinates, clean Cartesian points, unit normals, signed uniform draws, and observed Cartesian points at each of the six amplitudes. The delta=0.01 arrays are shared by the first-connectivity and before/at-closure observations. Array keys and units are documented in `matched_noise/config.json`. The archive was reconstructed from the existing seeds on 2026-10-15; it adds no independent sample or persistence calculation. `code/export_noise_coordinates.py` reproduces it, and `code/verify_release.py` checks every array against the recorded seed convention.

## Arc-support extension validation

`arc_support_verification.json` records 80 deterministic small samples and 5,767 closed graph budgets, including 3,796 supported budgets, 307 empty supported diagrams, 974 small-gap cases and 90 selected odd-field Betti checks. The new checker recomputes support on all 180 saved ellipse closure graphs; all pass and all their saved closure diagrams contain one bar. The theorem proof covers arbitrary fields. These validation samples are separate from the original 60 clean samples and are not additional independent replications of the main empirical study.

## Noise grid and selectors

`noise/` contains all two-model, six-amplitude, three-method, 20-seed grid measurements and both selectors. Clean observations are shared across the two models: there are 60 distinct clean method/seed cases, represented twice in the model-indexed tables. Both selectors choose the first closed tested graph in all 60. At delta=0.01, longest-lifetime selection gives median normalized diagram error 0.0085--0.029 across the six groups, compared with about 0.32 for first nonemptiness. The noise streams 20320101--20320120 are separate from clean-circle streams 20310101--20310120. Table C.1 reports the three lifetime-screen thresholds.

## Ellipse sweep and reconstruction

`ellipse/` has 3,516 graph evaluations in 180 method/sample series from 60 samples. All 3,099 strictly local graphs and all 180 exact first nonlocal events pass arc support and have one positive bar. Beyond the proved range, 60 graphs fail support and 40 have multiple bars, up to eight. No ordinary or capped error decrease occurs in the measured series at tolerance 1e-6. The cycle-specific Euclidean death is computed, and 59/60 samples pass both sufficient inequalities in Theorem 3(c); trough seed 20310103 fails the death margin. NN-crust independently reconstructs the ring on all 60 permuted coordinate sets.

## Exact example and presentation audit

`arc_support_example/exact_calculation.md` and `separation_verification.json` provide the six-point exact calculation and exhaustive orientation checks cited in Example 1 and Appendix D. `rank_scaling_audit.json` checks the closure fractions and Wilson intervals against every raw rank record. Profiles share underlying uniform draws and phase, explaining identical paired closure indicators for uniform/cosine and trough OR at n=1,000,000; these are paired observations.

## Phase 1 computations

`phase1/README.md` describes the completed Phase 1 evidence. Its `bridgeless/` folder contains the original-sample replay, the complete 30,000-record method sweep, Wilson intervals and timings. `euclidean/` holds all 220 distinct full diagrams and paired selector comparisons. Sample-level grid oracles, all-amplitude summaries, complete/displayed rank counts, literature checks, verification and the completion record are directly in `phase1/`. All paths and keys preserve the original pairing.

`separation_verification.json` now includes an exact radical metric certificate and every admissible cyclic order at all distinct thresholds. This strengthens the curve-order calculation: no single order works throughout the filtration, although each graph separately admits an order.
