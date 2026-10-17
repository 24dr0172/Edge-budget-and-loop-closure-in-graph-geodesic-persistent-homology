# Results guide

The paper separates four events: a graph becomes connected, it contains every consecutive pair around the clean curve, its persistence diagram becomes nonempty, and its diagram error becomes small. Here, **closure** means the second event. It is measured separately for each graph family.

File paths in this guide are relative to `results/` unless they begin with `code/`. Saved method names are `R2_global` (radius), `OR_kNN` and `mutual_kNN`. `R1_local` denotes CkNN in the order checks.

## Clean-circle results

`clean_closure/` contains 60 samples from three chosen density profiles, with 20 seeds per profile and n=512. The three graph families give 180 method/sample series.

The first connected graph at an exact edge count has an empty H1 diagram in all 60 samples per family. At the first connected whole integer-k layer, OR misses the loop in 38/60 samples and mutual in 40/60. These counts describe the chosen profiles; the larger sweep gives the profile-specific rates reported in the introduction and Table 3.

`clean_closure/event_table.csv` has 11 event rows. `clean_closure/knn_event_checks.csv` supplies each kNN family's own connectivity and closure events. `clean_closure/own_closure_cells.csv` marks closure separately for each family.

After closure, the saved series have no decreases in maximum distance error, RMS distance error or diagram error within the recorded tolerance. The distance result holds for graphs containing the ring. The diagram theorem always concerns capped error; ordinary error needs the additional conditions stated in the paper. Every retained clean sample satisfies the circle small-gap condition. Figure 2 shows the common budget grid. Table D.1 gives the first-closed-grid errors and exact-closure values.

## Edge-order checks

`order_checks/certificates.csv` reports 60/60 passes for radius, OR and mutual, and 59/60 for CkNN. The failing CkNN insertion is recorded in `order_checks/failure_witness.json`. `order_checks/precision_check.json` confirms a genuine violation of condition 2 by about 5.87e-6 using exact ordering and rigorous cost bounds. A failed check alone does not prove that a loop is lost or that extra bars occur.

The shuffled control has 0/60 passes, six loops at first connectivity before closure, and 16/60 first-nonempty/first-closed grid matches. CkNN is used for these condition checks, without a performance-ranking table or error curve.

## Additional neighbors and sample size

`extra_k/samples.csv` contains 120 paired measurements: 20 seeds, three profiles and two kNN rules. Each extra-k value is the sample's `k_close - k_conn`. Its median is computed from those differences.

For OR, the profile medians are 1, 1 and 0.5, with maximum 3 in each group. For mutual, the medians are 0.5, 1 and 7, with maxima 4, 3 and 26. Figure D.1 shows these distributions.

`extra_k/candidate_audit/` contains the rank-only sweep: 30,780 method records from 15,390 profile/seed samples and 5,130 shared seed/size combinations. It computes ranks without persistence diagrams. Table D.2 displays four sizes through n=32,768; the release also includes n=262,144 and n=1,000,000.

`extra_k/candidate_audit/main_text_summary.csv` contains all six Table 3 rows, including the connected bridge-free rates and pointwise Wilson intervals. The [rank-study guide](extra_k/candidate_audit/README.md) explains the count formula, sampling and probability model. Neighboring gap keys on one circle are dependent, so the independent-key model is a heuristic.

`rank_scaling_audit.json` checks the raw rates and intervals. `rank_pairing_replay.json` explains equal rates across paired profiles. At n=32,768, two opposite sample-level changes cancel in the uniform/trough OR closure rate. At the two larger sizes those OR indicators agree sample by sample; uniform and cosine also agree at n=1,000,000.

## Diagrams after closure

`useful_diagram/saved_prefixes.csv` contains 3,556 closed graphs and `saved_transitions.csv` contains 3,376 transitions. The budgets include saved closed grid points, exact closure and the complete graph. They do not include every intervening edge count. Each saved diagram has one positive bar, born at the longest consecutive chord within tolerance. The checkpoints store coordinates, reference diagrams and estimated diagrams.

`useful_diagram/rational_counterexample.json` gives the seven-point example in Remark 2 and Appendix D.2. One added edge lowers ordinary error from about 0.04406 to 0.03006, while capped error stays fixed. The verifier checks its bars and both bottleneck matching options independently.

## Noise and selection rules

`noise/` contains two models, six amplitudes, three graph methods and 20 seeds on the prescribed 24-budget grid. Extra connectivity and latent-closure budgets are diagnostics and are excluded from both selection rules and the grid oracle.

The models share zero-noise observations: there are 60 distinct clean method/seed cases. Both rules choose the first closed grid point in every clean case. At delta=0.01 both rules choose before latent closure in all 120 model/method/seed cases. Longest-lifetime selection has median normalized diagram error 0.0085–0.029 across the six groups, versus about 0.32 for first-nonempty selection. These are empirical comparisons; no general noisy-data guarantee is established.

`matched_noise/` provides radial checks at first connectivity and immediately before/at latent closure. The first-connectivity empty counts per method are 20, 9, 0, 0, 0 and 0 across the six amplitudes. None passes the half-reference-lifetime screen. At delta=0.01, every before/at-closure diagram is nonempty; all 20 pre-closure graphs per method pass the half-lifetime screen. Thus noisy graphs can have substantial loops before latent closure.

`matched_noise/coordinates.npz` stores latent coordinates, clean points, normals, signed draws and all observed coordinates. Its configuration gives keys and units. The same displacement vectors are scaled across amplitudes. The same delta=0.01 coordinates are used for both event comparisons.

## Ellipse and arc support

`ellipse_checks/` stores 60 coordinate sets and 180 method/sample closure checks. Its 539 distinct diagrams cover 540 events because one pair of events coincides. Appendix B gives the analytic constants and the separate loop-survival argument.

`ellipse/` contains 3,516 graph evaluations in 180 series. All 3,099 strictly local graphs and all 180 exact first nonlocal events pass arc support and have one positive bar. Farther beyond the proved range, 60 graphs fail support and 40 have multiple bars, with a maximum of eight.

The cycle's Euclidean death and the two sufficient inequalities in Theorem 3(c) are checked separately. They pass for 59/60 samples; trough seed 20310103 fails the death margin. No ordinary or capped error decrease occurs in the measured series at tolerance 1e-6. NN-crust recovers the ring on all 60 permuted samples.

`arc_support_verification.json` records 80 deterministic small samples and 5,767 closed graph budgets: 3,796 pass arc support, including 307 with empty diagrams. It also includes 974 small-gap cases and 90 selected odd-field checks. These tests check finite examples; the proof supplies the arbitrary-field result.

## Exact example and comparisons

The [exact Example 1 calculation](arc_support_example/exact_calculation.md) and `separation_verification.json` check all 120 cyclic orders across 14 threshold graphs. No single order works throughout the filtration, although every threshold graph has an order separately.

The [comparison guide](phase1/README.md) covers the 30,000-record connected bridge-free sweep, 220 distinct Euclidean diagrams, 720 sample-specific grid-oracle choices and all-amplitude summaries. The oracle minimizes error separately for each sample using the clean reference; it is a benchmark.

## Seeds, units and precision

Clean and ellipse studies use seeds 20310101–20310120. Noise uses separate seeds 20320101–20320120. Profiles share seed streams, and amplitudes share displacement vectors, so these paired observations should not be pooled as independent repetitions.

`clean_closure/clouds.csv` stores normalized reference lifetime in `oracle_longest_lifetime`. In `matched_noise/cells.csv` the same field is raw lifetime. `longest_over_reference` is dimensionless. Checkpoint diagram endpoints are unnormalized.

The [code guide](../code/README.md) gives precision, comparison tolerances, the separate witness-helper tolerance, package versions and runtimes. The release README gives the full reproduction sequence.
