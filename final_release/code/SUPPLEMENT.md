# Supplement: computations and supporting results

Shivam Singh, Annesha Sen and S. P. Tiwari. *Edge budget and loop closure in graph-geodesic persistent homology*.

Paths are relative to this release. The manuscript and its source notes are kept outside the release. The mathematical one-gap model and its proofs are in Appendix D of the paper.

## S1. Data map

- `results/clean_closure/`: clean-circle coordinates, graph events, budget cells, error summaries and checkpoints. `cells.csv` supplies the post-closure extremes; `clouds.csv` supplies reference lifetimes and diameters.
- `results/order_checks/`: CkNN order checks, the permuted-scale control and its diagrams, plus `precision_check.json`.
- `results/useful_diagram/`: saved closed-prefix diagrams, transitions and `rational_counterexample.json`.
- `results/extra_k/`: the original paired rank study. `candidate_audit/` holds the full size sweep, its seeds, rank checks, probability tails, model predictions and fixed-k summary.
- `results/matched_noise/`: radial connectivity and closure diagnostics, including the observed coordinates.
- `results/noise/`: both noise models, all amplitudes, graph budgets, diagrams and selector outputs.
- `results/comparisons/`: bridge-free graphs, unrestricted Euclidean diagrams, paired noise increments and grid-oracle summaries.
- `results/ellipse_checks/`: ellipse coordinates and connectivity/closure checks. `results/ellipse/` holds its post-closure sweep and NN-crust comparison.
- `code/config.json`: the seeds, formulas, grids, numerical conventions and package versions. `manuscript_figures/` supplies PDF and PNG assets.

Across the clean-circle grid, the largest post-closure diagram error is 0.0439884. The largest within-series increase is 0.0431358, from 0.0008526 to 0.0439884 for cosine OR, seed 20310119. Across the 60 clean samples, the empty-diagram cost P* ranges from 0.244 to 0.324. `results/number_audit.json` records the extremes; each P* is half the normalized reference lifetime in `results/clean_closure/clouds.csv`.

## S2. Software, precision and checks

Shortest paths use SciPy Dijkstra; persistence uses Ripser.py. Scores and distances use float64; returned persistence endpoints use float32. A diagram is empty exactly when no finite returned H1 interval has death strictly greater than birth. Comparison tolerances do not remove short positive intervals. The shortest stored positive lifetime is 1.2107193e-8 in raw metric units; this observed minimum is not a general resolution guarantee.

Raw comparison tolerances are 5e-8 for the clean grid, 1e-7 for selected-budget verification, and 1e-6 for the ellipse. The longest-lifetime selector returns the earliest grid budget within 1e-7 of the maximum raw lifetime. Diagnostic event budgets are excluded from selection and from the grid oracle. The support-check tolerance is 1e-10. Versions and settings are in `code/config.json` and `code/requirements.txt`.

The adaptive cutoff proof is in the remark after the arc-support proof in Appendix A.5. `code/validate_cutoff.py` compares sorted float32 endpoints with unrestricted Ripser on 150 random positive weighted graphs and 126 saved n=512 diagrams. Endpoints agree exactly in `results/cutoff_validation.json`.

`code/verify_arc_support.py` checks deterministic small graphs and 180 saved ellipse closure cases. Its independent F2 reduction of the Rips 2-skeleton checks ring death through triangle-boundary spans; selected Betti numbers are checked over F3 and F5. Tests cover nonempty and empty arc-supported diagrams and the six-point shortcut example. The report is `results/arc_support_verification.json`. The proof supplies the all-field claim.

## S3. Runtime and procedure cost

For a given edge order sigma, K budgets, n vertices and at most m graph edges, a direct implementation makes K all-pairs shortest-path and persistence calls. Binary-heap Dijkstra costs O(K n (m+n) log n), in addition to persistence, and distance storage costs O(n^2). A naive 2-skeleton has O(n^3) simplices: n=512 permits 22,238,720 triangles per budget. Ripser avoids constructing this complete simplex list explicitly. First-nonempty selection can stop at the first positive diagram; longest-lifetime selection evaluates the full grid.

Full-diagram studies use at most n=512. The size sweep through one million points computes ranks and bridge events, without persistence diagrams. Representative serial blocks cover 24 budgets and all three methods; median runtimes are about 131-196 seconds across four cases. These are block times, not the total study time.

The machine is an AMD Ryzen 7 7445HS with Python 3.12.3. OPENBLAS, OMP and MKL thread limits are 1. Each of four cases is measured three times after a three-vertex Ripser warmup. Timed blocks include pair/rank orders, graph construction, all-pairs Dijkstra and finite H1 diagrams. They exclude sampling, reference construction, bottleneck errors, certificates, plotting and file I/O. `results/runtime_benchmark.json` records the platform, package versions and all measurements; `code/benchmark_diagrams.py` repeats the benchmark. Timings depend on hardware and load.

## S4. The CkNN certificate

The exception is the trough sample with seed 20310105. Using zero-based indices, the outer chord from 190 to 201 costs 0.01730770437513113; the consecutive route from 190 to the interior vertex 200 costs 0.01731357612075681. The opposite route costs about 0.0000022653.

Treating the saved float64 coordinates as exact dyadic numbers, the fourth power of a CkNN score is a rational function of squared chords. Exact ordering puts the outer edge at insertion 5840 and the required direct edge from 190 to 200 at 5841. Integer-square-root bounds at 100 decimal places enclose the route excess in [5.87174562568017, 5.87174562568018] x 1e-6, strictly above zero. Independent 60- and 100-digit calculations on the ideal circle, using saved arc parameters, give about 5.87174562567631e-6. The violation is therefore genuine. `code/check_cknn_precision.py` writes `results/order_checks/precision_check.json`.

## S5. Exact examples

`code/verify_separation.py` verifies Example 1 in all 120 cyclic orders at all 14 distinct threshold graphs, keeping equal distances together. It enumerates orientations directly at the two proof thresholds, including all 1,024 orientations at 19/10. The per-threshold admissible orders and their empty intersection are in `results/separation_verification.json`.

The complete seven-point capped-error example and its two bottleneck matchings are at the end of Appendix A.6. The coordinates, bars and transition are in `results/useful_diagram/rational_counterexample.json`. The release verifier checks them with a separate F2 reduction and bottleneck matching.

## S6. Euclidean baseline and noise comparisons

`code/euclidean_baseline.py` runs unrestricted Ripser on the 220 distinct observed clouds and retains every finite positive H1 interval. Clean clouds are shared across the two noise models, giving 240 model/amplitude/seed rows. Diagrams, timings, input hashes and paired comparisons are in `results/comparisons/euclidean/`.

For each seed and pipeline, subtract its own clean error from its noisy error before computing medians, quartiles and paired counts. `paired_selectors.csv` contains the paired absolute errors; `paired_summary.csv` includes noise increments. The longest-lifetime rule selects its budget separately at each amplitude. The clean Euclidean median error of about 0.1154 against the intrinsic reference is therefore separated from the noise increment.

`code/comparison_summaries.py` selects only connected rows with `grid=True` in `results/noise/cells.csv` to compute the reference-based grid oracle. It also summarizes all six amplitudes. The same displacement vector is scaled across amplitudes within each seed and noise model. These amplitude comparisons are paired. All comparisons use the same reference, normalization and endpoint convention.

## S7. Rank sweep and bridge-free protocol

The primary sizes are 128, 512, 4,096 and 32,768, with 1,500 seeds per profile at the first three sizes and 500 at the last. This gives 5,000 seed/size combinations, 15,000 profile/seed samples and 30,000 method records. The full sweep adds 262,144 with 100 seeds per profile and 1,000,000 with 30, giving 5,130 seed/size combinations, 15,390 profile/seed samples and 30,780 method records. The full bridge-free sweep uses these same records at every size. Seed settings are in `results/extra_k/candidate_audit/config.json` and `code/config.json`.

Profiles share sorted uniform draws and phase for each seed and size; their inverse CDFs produce the three profiles. Equation (6) computes consecutive keys directly. Its 18,980 keys on 240 small samples agree with full chord ranks; 480 connectivity checks agree with the second-largest key. All sweep samples have largest gap below one third, so the circle theorem turns rank closure into a positive H1 diagram. The sweep does not recompute persistence.

Bridge-free checks rebuild full OR and mutual adjacency using `code/bridge_free.py` and `code/bridgeless_circle.cpp`. All edges in a whole integer-k layer enter before connectedness and bridges are tested separately. Consecutive rank keys do not determine bridges. Independent graph traversals check the optimized calculation on random small graphs and all 60 original samples using full Cartesian ranks. Measurements, validation, timings and pointwise Wilson intervals are in `results/comparisons/bridgeless/`.

Equal profile percentages can arise from paired draws and the piecewise-affine inverse trough CDF. A window contained in one affine branch preserves its two outer counts, including in the denser region. At n=32,768, two opposite seed-level changes cancel in the OR connection rate and all +3 indicators agree; most closing gaps lie outside the trough. `results/rank_pairing_replay.json` records the replay. At n=262,144 and n=1,000,000, the OR closure-at-connection and +3 indicators agree for every paired uniform/trough seed, as checked by `code/audit_rank_scaling.py`.

## S8. Probability calculations and small-study calibration

`code/check_knn_miss_candidate.py` evaluates the one-gap tails by exact rational sums before conversion and checks 12 values against numerical integration (maximum difference 1.11e-16). `results/extra_k/candidate_audit/poisson_tails.csv` and `iid_predictions.csv` hold the tails and independent-key predictions. Neighboring circle keys are dependent; the independent-key model is a heuristic for the sampled graphs.

The original 60 samples close at first connected integer k in 22 OR and 20 mutual cases. Among the larger sweep's 1,500 n=512 seed triples, the empirical probabilities for mutual closure in 0, 1, 2 or 3 profiles are (0.592667, 0.204, 0.183333, 0.020). Convolving this distribution 20 times gives mean 12.6133 and upper tail 0.04033 at 20. OR gives mean 19.76 and upper tail 0.3701 at 22. This descriptive calibration preserves the three profiles within each seed; it changes none of the observed counts. `results/comparisons/small_study_calibration.json` stores the full distributions. Regenerate it with `python3 tools/generate_supplement.py`.

## S9. Supporting tables and figures

**Table S1. Median normalized errors**, with 20 samples per profile. Exact closure uses M=M_close; the first closed grid point is the smallest tested budget containing every consecutive pair. The edge counts can differ between methods.

| Profile | Order | Exact e_inf | Exact e_2 | Exact e_PH | Grid e_inf | Grid e_2 | Grid e_PH |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Uniform | Radius | 0.00022 | 0.00013 | 0.00015 | 0.00026 | 0.00014 | 0.00017 |
| Uniform | OR kNN | 0.00017 | 0.00008 | 0.00010 | 0.00021 | 0.00010 | 0.00012 |
| Uniform | Mutual kNN | 0.00017 | 0.00009 | 0.00010 | 0.00019 | 0.00010 | 0.00011 |
| Cosine | Radius | 0.00151 | 0.00074 | 0.00090 | 0.00164 | 0.00081 | 0.00102 |
| Cosine | OR kNN | 0.00113 | 0.00031 | 0.00039 | 0.00137 | 0.00037 | 0.00048 |
| Cosine | Mutual kNN | 0.00102 | 0.00028 | 0.00036 | 0.00115 | 0.00030 | 0.00045 |
| Trough | Radius | 0.00257 | 0.00135 | 0.00164 | 0.00317 | 0.00164 | 0.00190 |
| Trough | OR kNN | 0.00216 | 0.00075 | 0.00064 | 0.00234 | 0.00083 | 0.00079 |
| Trough | Mutual kNN | 0.00287 | 0.00102 | 0.00102 | 0.00306 | 0.00106 | 0.00112 |

**Table S2. Full size sweep.** Median integer k values and percent reaching closure at connectivity, at the first connected bridge-free layer, and at k_conn+3. Brackets are pointwise 95% Wilson intervals. Counts are per profile and rule; profiles and rules are paired. The machine-readable table is `results/comparisons/size_sweep_summary.csv`.

| Profile | n | Order | Samples | Median k_conn | Median k_br | Median k_close | At connection (%) | At k_br (%) | +3 (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Uniform | 128 | OR kNN | 1500 | 5 | 6 | 6 | 34.3 [32.0, 36.8] | 74.5 [72.2, 76.6] | 95.1 [93.9, 96.1] |
| Uniform | 128 | Mutual kNN | 1500 | 8 | 9 | 9 | 28.7 [26.4, 31.0] | 70.0 [67.6, 72.3] | 90.1 [88.4, 91.5] |
| Uniform | 512 | OR kNN | 1500 | 7 | 7 | 8 | 32.9 [30.5, 35.3] | 74.9 [72.7, 77.1] | 94.7 [93.5, 95.8] |
| Uniform | 512 | Mutual kNN | 1500 | 10 | 11 | 11 | 28.9 [26.7, 31.3] | 71.7 [69.4, 74.0] | 91.1 [89.5, 92.4] |
| Uniform | 4,096 | OR kNN | 1500 | 9 | 10 | 10 | 32.2 [29.9, 34.6] | 74.8 [72.5, 76.9] | 95.5 [94.3, 96.4] |
| Uniform | 4,096 | Mutual kNN | 1500 | 13 | 14 | 14 | 28.7 [26.5, 31.1] | 71.8 [69.5, 74.0] | 91.0 [89.4, 92.3] |
| Uniform | 32,768 | OR kNN | 500 | 12 | 12 | 13 | 31.2 [27.3, 35.4] | 74.6 [70.6, 78.2] | 96.6 [94.6, 97.9] |
| Uniform | 32,768 | Mutual kNN | 500 | 16 | 17 | 17 | 29.4 [25.6, 33.5] | 67.8 [63.6, 71.7] | 90.8 [87.9, 93.0] |
| Uniform | 262,144 | OR kNN | 100 | 14 | 15 | 15 | 27.0 [19.3, 36.4] | 81.0 [72.2, 87.5] | 93.0 [86.3, 96.6] |
| Uniform | 262,144 | Mutual kNN | 100 | 19 | 20 | 20 | 26.0 [18.4, 35.4] | 71.0 [61.5, 79.0] | 88.0 [80.2, 93.0] |
| Uniform | 1,000,000 | OR kNN | 30 | 16 | 16 | 17 | 46.7 [30.2, 63.9] | 76.7 [59.1, 88.2] | 96.7 [83.3, 99.4] |
| Uniform | 1,000,000 | Mutual kNN | 30 | 21 | 22 | 22 | 23.3 [11.8, 40.9] | 66.7 [48.8, 80.8] | 86.7 [70.3, 94.7] |
| Cosine | 128 | OR kNN | 1500 | 5 | 6 | 6 | 33.7 [31.4, 36.2] | 72.9 [70.6, 75.1] | 92.1 [90.6, 93.3] |
| Cosine | 128 | Mutual kNN | 1500 | 8 | 9 | 10 | 23.2 [21.1, 25.4] | 64.7 [62.3, 67.1] | 78.7 [76.6, 80.7] |
| Cosine | 512 | OR kNN | 1500 | 7 | 7 | 8 | 33.2 [30.9, 35.6] | 74.1 [71.8, 76.2] | 94.1 [92.8, 95.2] |
| Cosine | 512 | Mutual kNN | 1500 | 10 | 11 | 11 | 27.0 [24.8, 29.3] | 71.4 [69.1, 73.6] | 89.2 [87.5, 90.7] |
| Cosine | 4,096 | OR kNN | 1500 | 9 | 10 | 10 | 31.8 [29.5, 34.2] | 75.2 [73.0, 77.3] | 95.4 [94.2, 96.3] |
| Cosine | 4,096 | Mutual kNN | 1500 | 13 | 14 | 14 | 28.2 [26.0, 30.5] | 72.1 [69.8, 74.3] | 91.3 [89.8, 92.7] |
| Cosine | 32,768 | OR kNN | 500 | 12 | 12 | 13 | 30.8 [26.9, 35.0] | 73.6 [69.6, 77.3] | 96.6 [94.6, 97.9] |
| Cosine | 32,768 | Mutual kNN | 500 | 16 | 17 | 17 | 29.4 [25.6, 33.5] | 68.0 [63.8, 71.9] | 90.6 [87.7, 92.9] |
| Cosine | 262,144 | OR kNN | 100 | 14 | 15 | 15 | 26.0 [18.4, 35.4] | 81.0 [72.2, 87.5] | 93.0 [86.3, 96.6] |
| Cosine | 262,144 | Mutual kNN | 100 | 19 | 20 | 20 | 26.0 [18.4, 35.4] | 71.0 [61.5, 79.0] | 88.0 [80.2, 93.0] |
| Cosine | 1,000,000 | OR kNN | 30 | 16 | 16 | 17 | 46.7 [30.2, 63.9] | 76.7 [59.1, 88.2] | 96.7 [83.3, 99.4] |
| Cosine | 1,000,000 | Mutual kNN | 30 | 21 | 22 | 22 | 23.3 [11.8, 40.9] | 66.7 [48.8, 80.8] | 86.7 [70.3, 94.7] |
| Trough | 128 | OR kNN | 1500 | 5 | 6 | 7 | 24.2 [22.1, 26.4] | 55.9 [53.4, 58.4] | 71.1 [68.8, 73.4] |
| Trough | 128 | Mutual kNN | 1500 | 10 | 12 | 18 | 6.8 [5.6, 8.2] | 37.7 [35.3, 40.2] | 36.5 [34.1, 39.0] |
| Trough | 512 | OR kNN | 1500 | 7 | 7 | 8 | 32.7 [30.4, 35.1] | 74.1 [71.8, 76.2] | 94.1 [92.8, 95.2] |
| Trough | 512 | Mutual kNN | 1500 | 12 | 15 | 19 | 7.1 [5.9, 8.5] | 41.2 [38.7, 43.7] | 38.1 [35.6, 40.6] |
| Trough | 4,096 | OR kNN | 1500 | 9 | 10 | 10 | 32.3 [29.9, 34.7] | 74.6 [72.3, 76.7] | 95.5 [94.3, 96.4] |
| Trough | 4,096 | Mutual kNN | 1500 | 14 | 16 | 20 | 11.5 [10.0, 13.2] | 43.9 [41.4, 46.5] | 44.6 [42.1, 47.1] |
| Trough | 32,768 | OR kNN | 500 | 12 | 12 | 13 | 31.2 [27.3, 35.4] | 74.6 [70.6, 78.2] | 96.6 [94.6, 97.9] |
| Trough | 32,768 | Mutual kNN | 500 | 17 | 18.5 | 21 | 14.8 [12.0, 18.2] | 47.4 [43.1, 51.8] | 52.8 [48.4, 57.1] |
| Trough | 262,144 | OR kNN | 100 | 14 | 15 | 15 | 27.0 [19.3, 36.4] | 81.0 [72.2, 87.5] | 93.0 [86.3, 96.6] |
| Trough | 262,144 | Mutual kNN | 100 | 19 | 21 | 23 | 12.0 [7.0, 19.8] | 47.0 [37.5, 56.7] | 61.0 [51.2, 70.0] |
| Trough | 1,000,000 | OR kNN | 30 | 16 | 16 | 17 | 46.7 [30.2, 63.9] | 76.7 [59.1, 88.2] | 96.7 [83.3, 99.4] |
| Trough | 1,000,000 | Mutual kNN | 30 | 22 | 23 | 24 | 16.7 [7.3, 33.6] | 50.0 [33.2, 66.8] | 60.0 [42.3, 75.4] |

**Figure S1. Paired additional-neighbor distributions.** Each curve uses 20 samples per profile at n=512. A jump at zero gives closure at first connectivity. Horizontal scales differ to show the mutual trough tail.

![Figure S1](manuscript_figures/FigureS1.png)

[PDF](manuscript_figures/FigureS1.pdf) | [PNG](manuscript_figures/FigureS1.png). Generated from `results/extra_k/samples.csv` by `code/build_figures.py`.

**Figure S2. Radial first-connectivity diagnostic.** Empty-diagram fraction and longest lifetime relative to the clean reference for radius and OR graphs.

![Figure S2](manuscript_figures/FigureS2.png)

[PDF](manuscript_figures/FigureS2.pdf) | [PNG](manuscript_figures/FigureS2.png). Generated from `results/matched_noise/summary.csv` by `code/build_figures.py`.

**Figure S3. Noise error curves.** Diagram, worst-pair and RMS errors for both noise models at amplitudes 0, 0.001 and 0.01.

![Figure S3](manuscript_figures/FigureS3.png)

[PDF](manuscript_figures/FigureS3.pdf) | [PNG](manuscript_figures/FigureS3.png). Generated from `results/noise/curves.csv` by `code/summarize_studies.py`.

**Figure S4. Unsupported-shortcut illustration.** A separate rendering of the six-point construction discussed in Remark 1. The paper's Figure 1 is drawn in its LaTeX source.

![Figure S4](manuscript_figures/FigureS4.png)

[PDF](manuscript_figures/FigureS4.pdf) | [PNG](manuscript_figures/FigureS4.png). Generated by `code/summarize_studies.py`.
