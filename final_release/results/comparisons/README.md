# Study comparisons

This folder holds the bridge-free study, the Euclidean baseline and the grid oracle. The paper cites two bridge-free proportions at n=512 (Section 7.2). The other comparisons are additional studies that the paper does not rely on. The exact certificate for Example 1 is linked below for convenience. Links are relative to this folder.

## Evidence

- [Exact Example 1 certificate](../separation_verification.json): all 120 cyclic orders at 14 distinct threshold graphs. Each graph admits an order separately; no order works throughout the filtration. The certificate records distances and admissible orders, including the thresholds 3/2 and 19/10 inside the positive interval. Table 3 is exported separately to [rank_table.csv](../rank_table.csv); it is not the bridge-free summary.
- [Clean bridge events](bridgeless/clean_samples.csv): 60 samples and two kNN rules, checked against full Cartesian ranks and an independent bridge search. The [summary](bridgeless/clean_summary.csv) gives rates and Wilson intervals.
- [Bridge sweep](bridgeless/sweep_samples.csv): 15,390 profile/seed/size samples and 30,780 method records through n=1,000,000. See the [summary](bridgeless/sweep_summary.csv), [validation](bridgeless/validation.json) and [timings](bridgeless/sweep_timings.csv).
- [Euclidean baseline](euclidean/samples.csv): 220 distinct observed clouds. Reusing clean clouds across models gives 240 model/amplitude/seed rows. The baseline is also shared across graph methods. See [error quartiles](euclidean/summary.csv) and [paired comparisons](euclidean/paired_summary.csv).
- [Grid oracle](grid_oracle_samples.csv): the lowest error among connected points on the prescribed 24-budget grid, chosen separately for each of 720 cases. Extra diagnostic budgets are excluded. See [group summaries](grid_oracle_summary.csv) and [selector gaps](selected_vs_oracle.csv). This benchmark uses the clean reference.
- [All amplitudes](all_amplitudes.csv): both rules, all six amplitudes, error medians/quartiles and median normalized budgets. Samples and displacements remain paired.
- [Sample counts](rank_display_counts.json): the primary display has 5,000 shared seed/size combinations, 15,000 profile/seed samples and 30,000 method records. The full release has 30,780 method records.
- [Literature sources](literature_checks.json): source details for the connectivity, spectral-distance and cyclic-graph comparisons in the paper.

## Connected bridge-free graphs

The test inserts a whole integer-k layer, then checks that the graph is connected and has no bridge. A bridge is an edge whose removal disconnects the graph. Every record satisfies k_conn ≤ k_br ≤ k_close. The bridge-free test can still stop before closure.

At n=512, each profile has 1,500 paired seeds. Closure rates at connection and at the first connected bridge-free layer are:

- Uniform, OR kNN: 32.9% at connection; 74.9% at the bridge-free layer.
- Uniform, mutual kNN: 28.9% at connection; 71.7% at the bridge-free layer.
- Cosine, OR kNN: 33.2% at connection; 74.1% at the bridge-free layer.
- Cosine, mutual kNN: 27.0% at connection; 71.4% at the bridge-free layer.
- Trough, OR kNN: 32.7% at connection; 74.1% at the bridge-free layer.
- Trough, mutual kNN: 7.1% at connection; 41.2% at the bridge-free layer.

The CSV summaries give pointwise 95% Wilson intervals. Profiles and rules reuse seeds, so pooled rows are not independent repetitions.

## Noise comparison at delta=0.01

The clean Euclidean median error is about 0.1154 against the intrinsic reference. To compare noise sensitivity, subtract each pipeline’s own clean error from its noisy error for each seed before computing medians and paired counts. The longest-lifetime rule selects its budget separately at each amplitude. Euclidean increments are shared across graph methods.

- Isotropic, OR kNN: median noise increment 0.02857 for longest-lifetime selection and 0.02731 for Euclidean; selection has the smaller increment in 9/20 paired samples. Grid-oracle absolute error: 0.00522.
- Isotropic, radius: median noise increment 0.00829 for longest-lifetime selection and 0.02731 for Euclidean; selection has the smaller increment in 19/20 paired samples. Grid-oracle absolute error: 0.00603.
- Isotropic, mutual kNN: median noise increment 0.01827 for longest-lifetime selection and 0.02731 for Euclidean; selection has the smaller increment in 15/20 paired samples. Grid-oracle absolute error: 0.00545.
- Radial, OR kNN: median noise increment 0.02282 for longest-lifetime selection and 0.02958 for Euclidean; selection has the smaller increment in 11/20 paired samples. Grid-oracle absolute error: 0.00632.
- Radial, radius: median noise increment 0.01001 for longest-lifetime selection and 0.02958 for Euclidean; selection has the smaller increment in 20/20 paired samples. Grid-oracle absolute error: 0.00805.
- Radial, mutual kNN: median noise increment 0.01454 for longest-lifetime selection and 0.02958 for Euclidean; selection has the smaller increment in 15/20 paired samples. Grid-oracle absolute error: 0.00645.

These results describe the tested models and grid. The selection rules have no general guarantee under noise. The grid oracle uses the clean reference and serves as a benchmark.

## Reproduction

Run commands from the release folder with the environment and complete sequence in the [code guide](../../code/README.md). To rebuild these summaries and this guide from saved measurements, run:

```sh
python3 code/reproduce_results.py
python3 code/comparison_report.py --data-only
python3 code/verify_comparisons.py
python3 tools/record_current.py --data-only
python3 code/verify_release.py --data-only
```

To rerun the bridge comparison and Euclidean baseline, use:

```sh
python3 code/bridge_free.py --validate --clean --sweep
python3 code/euclidean_baseline.py
python3 code/reproduce_results.py
python3 code/comparison_report.py --data-only
python3 code/verify_comparisons.py
python3 tools/record_current.py --data-only
python3 code/verify_release.py --data-only
```

The bridge calculation requires g++ with C++17 support. The helper is compiled in a temporary folder. Verification uses saved evidence and needs no C++ compiler. The [code guide](../../code/README.md) explains checkpoint reuse and fresh runs.

## Timing definitions

Individual Euclidean times include distances and persistence, excluding checkpoint file I/O. The per-checkpoint times cover all 220 diagrams.

Bridge times include sampling, helper startup, graph reconstruction and event calculations, excluding saved-row writing. Child RSS is the cumulative maximum across helper processes. The pilot also includes compiler memory when compilation was needed.

Timing records describe the runs that produced the saved measurements and depend on the computer and its load. Rebuilding this guide does not measure time again.
