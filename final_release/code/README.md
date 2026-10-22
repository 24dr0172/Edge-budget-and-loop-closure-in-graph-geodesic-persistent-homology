# Code

This folder holds the scripts that generate, summarize and check the data. The root [README](../README.md) maps every table and figure of the paper to its data and script, and defines the terms used below. Run all commands from the release folder (the folder that contains `code/` and `tools/`), with Python 3.12 and `code/requirements.txt` installed.

## Method and selector names in the data

The data files use short code names. They map to the paper as follows.

| Code name | Meaning in the paper |
|---|---|
| `R2_global` | radius order (pairs sorted by straight-line distance) |
| `OR_kNN` | OR-kNN order (a pair enters when either point ranks the other among its k nearest neighbors) |
| `mutual_kNN` | mutual-kNN order (both points must rank each other) |
| `R1_local` | the locally scaled CkNN order (additional study) |
| `R3_permuted` | the same local scales shuffled among vertices, an additional-study control |
| `first_nonempty` | first-nonempty selector (additional study) |
| `longest_lifetime` | earliest-longest-bar selector (additional study) |

The prefixes `R1_`, `R2_` and `R3_` in the method names are labels only; they are unrelated to the selection rules (R1) and (R2).

## Check the saved data

```sh
python3 code/verify_release.py --paper-only --data-only --report /tmp/edge-budget-paper-check.json
```

This runs the checks behind the paper's claims and checks the hash of every release file. The additional studies are reported as not run. Write this report outside the release. To compare with the manuscript as well, replace `--data-only` with `--manuscript-dir /path/to/paper`.

```sh
python3 code/verify_release.py --data-only
python3 code/verify_release.py --manuscript-dir /path/to/paper
```

These run every check, including the additional studies, all CSV column lists and every file hash. They replay seeded inputs and check saved diagrams, errors, selection rules, rank summaries, curve-order records and the exact examples. They do not regenerate the large persistence studies.

## Rebuild paper outputs from saved data

```sh
python3 code/reproduce_results.py --paper-table-only
python3 code/verify_separation.py
python3 code/verify_natural_instance.py
python3 code/verify_small_examples.py
python3 code/check_practical_claims.py
python3 code/build_retention_table.py
python3 code/plot_post_closure_error.py
python3 code/retention_extra_layers.py
python3 code/build_practical_choices.py
python3 code/peanut_neck_study.py --plot-only
python3 code/check_seven_point_example.py
python3 code/verify_retention_cost_witness.py --write
python3 code/verify_extra_layers_and_peanut.py
```

What each script does:

- `reproduce_results.py --paper-table-only` writes `results/rank_table.csv`, the 24 cells of Table 3 with counts, unrounded percentages and Wilson intervals. `rank_table.py` holds its checks: group sizes, shared seeds, gaps below 1/3, saved fractions and Wilson intervals. Without the option, `reproduce_results.py` rebuilds every summary table, including those of the additional studies.
- `verify_separation.py` certifies Example 1 with exact arithmetic on square roots, at both thresholds 3/2 and 19/10.
- `verify_natural_instance.py` certifies Example 2 with exact rational interval arithmetic (Machin's formula for pi and Taylor bounds, no extra packages). It also checks that each strict comparison has a margin of at least 0.03 raw units, and how far each rounded decimal bound is from its true value.
- `verify_small_examples.py` checks Remarks 1 and 2, both separation examples and the tie cases of Appendix B. It uses its own finite boundary reduction and reads no study data.
- `check_practical_claims.py` checks the fixed-k, bridge-free, independent-key and diagram-error numbers of Sections 7 and 8 from the raw records. It writes `results/practical_claims_verification.json`.
- `build_retention_table.py` writes the compact Table 3 CSV, three-layer CSV and fixed-k cells to `results/retention_table/`.
- `build_practical_choices.py` writes an additional table, comparing retention at first connectivity and the first bridge-free layer, and an additional table, partitioning fixed-k outcomes.
- `peanut_neck_study.py --plot-only` rebuilds the additional varying-waist figure from saved records; without this option it generates the study from its 20 paired seeds.
- `plot_post_closure_error.py` writes the additional post-closure error figure and its checks to `results/additional_figures/`.
- `retention_extra_layers.py` writes Figure 3 and its counts for all four sample sizes.
- `check_seven_point_example.py` recomputes every number of Appendix C and writes `results/seven_point_check/seven_point_numbers.json`. It numbers the points from 0, so its point i is x_(i+1) in the paper.
- `verify_independent_checks.py` independently checks the numbers of Appendix C, the Section 8.2 independent-key values and the geometry of the additional two-panel and ellipse drawings. `verify_retention_cost_witness.py --write` saves its output to `results/independent_checks/`, then recomputes every Table 3 cell, fixed-k cell, the additional post-closure error figure summary value and Appendix C value from the raw records.
- `verify_extra_layers_and_peanut.py` checks the Figure 3 counts and the additional peanut/trough samples, including small independent diagram replays.
- `verify_reported_claims.py` compares the capped-error formulas with the saved circle and ellipse errors and checks the additional CkNN, fixed-order and noise numbers. It is run by complete verification.

None of these scripts changes a measurement. `build_retention_table.py` and `plot_post_closure_error.py` accept a `RELEASE_ROOT` environment variable to work on another copy of the release. Native drawing sources and the Overleaf project are kept outside this release.

`reproduce_results.py --rank-table-only` exports the six-row n=512 comparison at connection, at the first bridge-free layer and at k_conn+3. `--fixed-k-only` rebuilds the uniform fixed-k summary.

## Record and verify after a change

After changing code, documentation or any tracked output, rebuild in this order:

1. Rebuild the outputs:

   ```sh
   python3 code/reproduce_results.py --paper-table-only
   python3 code/verify_separation.py
   python3 code/verify_natural_instance.py
   python3 code/verify_small_examples.py
   python3 code/check_practical_claims.py
   python3 code/build_retention_table.py
   python3 code/plot_post_closure_error.py
   python3 code/retention_extra_layers.py
   python3 code/build_practical_choices.py
   python3 code/peanut_neck_study.py --plot-only
   python3 code/check_seven_point_example.py
   python3 code/verify_retention_cost_witness.py --write
   python3 code/verify_extra_layers_and_peanut.py
   python3 code/comparison_report.py --data-only
   ```

2. Rewrite the guides and check the study comparisons:

   ```sh
   python3 tools/write_guides.py
   python3 code/verify_comparisons.py
   ```

3. Record the file hashes, then verify:

   ```sh
   python3 tools/record_current.py --data-only
   python3 code/verify_release.py --data-only
   ```

To include the manuscript, replace `--data-only` in steps 1 and 3 with `--manuscript-dir /path/to/paper`. With `--data-only`, the saved manuscript records are kept unchanged. To rebuild the summaries and figures of the additional studies, also run `python3 code/reproduce_results.py` and `python3 code/build_figures.py` before step 3.

`results/file_manifest.json` lists the SHA-256 hash of every release file except itself, `results/verification_report.json` and `results/run_logs/`. Do not change any listed file after recording. Keep the manuscript outside the release folder.

## Full generation from seeds

First check that every script, drawing, pinned package and the C++ compiler are present:

```sh
python3 code/reproduce.py --check
```

This reads no saved measurements and changes nothing. Then run:

```sh
python3 code/reproduce.py --output /path/to/new-run --workers 6
```

This copies the code, tools, documents and the TikZ drawings into the empty folder. The new `results/` folder starts empty. The pipeline generates every study, including the two-curve locality check, then builds the documents, records the file hashes and runs the complete verification. Add `--manuscript-dir /path/to/paper` to include the manuscript comparison.

Allow several hours and several GB of memory; lower `--workers` if needed. Continue an interrupted run with `--resume` (the code must be unchanged). Logs for each stage are in `results/run_logs/`. Timings depend on the computer and are not reproduction targets.

`python3 code/smoke_test.py` checks one seeded circle with n=128 points. It reads and writes no release files.

## Generators of the additional studies

- `generate_clean.py`, `generate_measurements.py`: the clean-circle study at n=512 (events, errors, CkNN and the shuffled-scale control), with full diagrams.
- `rank_sweep.py`: the rank sweep over six sample sizes and the probability checks. The paper uses the first four sizes.
- `bridge_free.py` with `bridgeless_circle.cpp`: connected bridge-free layers up to n=1,000,000. Needs a C++17 compiler.
- `run_studies.py`, `generate_matched_noise.py`, `euclidean_baseline.py`: noise studies, selection rules, Euclidean baselines and the grid oracle.
- `comparison_summaries.py` (called by `reproduce_results.py`): grid-oracle and selector summaries, rank display counts, the 36-row size sweep summary (`results/comparisons/size_sweep_summary.csv`: median k_conn, k_br and k_close, and closure at connection, at the first bridge-free layer and at k_conn+3, with Wilson intervals) and the paired calibration of the small study (`results/comparisons/small_study_calibration.json`).
- `generate_ellipse.py`: ellipse coordinates and event checks.
- `verify_arc_support.py`: small deterministic graphs and saved ellipses. Its support and homology functions are also used by the small examples; importing them runs nothing. Run `python3 code/verify_arc_support.py --report results/arc_support_verification.json` for a full report.
- `verify_curve_order.py`, `check_cknn_precision.py`, `check_rank_pairing.py`, `validate_cutoff.py`: further diagnostics, and the test of the persistence cutoff.
- `check_noncircular_locality.py`: the numerical check of the locality proposition of Appendix A on a 2:1 ellipse and a peanut curve, with a round circle as a control (about ten minutes on several cores). Its thresholds are grid values, not certified bounds; `beta` is the code name for the antipodal-chord bound.
- `check_peanut_trough.py`: 450 peanut clouds with trough density at n=128, 256 and 512, all three orders, full diagrams at exact closure. Settings and code hashes are in `results/peanut_trough/config.json`.
- `peanut_neck_study.py`: noise-free varying-waist study at n=512, paired across widths and methods, with unrestricted diagrams, a circle control and an intrinsic reference. `verify_peanut_neck.py` checks all saved bars and errors, cone certificates, and nine independent event-diagram replays using Floyd-Warshall distances. It also recomputes all 240 immediately preceding peanut prefixes to check that a single cross-neck insertion changes one bar alive at threshold 1/4 into two bars that die below it.
- `benchmark_diagrams.py`: timing of representative full-diagram evaluations.
- `audit_numbers.py`, `audit_rank_scaling.py`: recount reported numbers and size-sweep fractions from the raw records.

Shared definitions are in `graph_methods.py`, `experiment_core.py` and `finite_diagram.py`.

## Numerical conventions

`config.json` records the seeds, profiles, grids, tolerances and package versions. `result_schemas.json` lists the expected CSV columns. In the rank sweep, the three profiles share uniform draws and phase. Noise amplitudes share displacement draws within a seed and model, and the two noise models share clean samples. Repeated clean rows are therefore not independent samples.

Distances and shortest paths use float64. Ripser returns float32 endpoints. Every finite H1 interval with death greater than birth is kept. The small-example reducer uses a 1e-12 tolerance for positive lifetimes. Numerical checks support, but do not replace, the proofs in the paper.

Selection rules use connected grid graphs with a positive interval. If there is none, they return no budget and an infinite loss. Rule (R2) takes the earliest eligible budget whose longest lifetime is within 1e-7 of the maximum, in raw distance units. Extra diagnostic budgets are excluded from selection and from the grid oracle. Grid coefficients are rounded to four decimals before budgets are formed; at n=512 the coefficient 12 gives budget 7,302. Tolerances on normalized loss are 5e-8 for the clean grid, 1e-7 for selected budgets and 1e-6 for ellipse studies. Endpoint tolerances are in raw distance units.

**Persistence cutoff.** `finite_graph_diagram` in `run_studies.py` computes complete H1 diagrams of graph metrics with a cutoff, and keeps every positive bar it returns.

Why the cutoff is safe:

- Let w_max be the largest weight of a generating edge and K_t the Rips complex of the graph metric at threshold t. Take w_max <= t0 <= t1.
- For an edge xy of K_t1, choose a shortest graph path x=p_0,...,p_k=y. Each triangle x p_i p_(i+1) lies in K_t1, because D(x,p_i) <= D(x,y) <= t1 and D(p_i,p_(i+1)) <= w_max. So xy is homologous to that path.
- Every cycle of K_t1 is therefore homologous to a cycle of generating edges, which already lies in K_t0. Hence no H1 class is born above w_max.
- Once a cutoff above w_max leaves no H1 class alive, the diagram is complete.

The routine first tries the float32 successor of w_max. If a class is still alive, it raises the cutoff, up to the smallest vertex eccentricity (there one vertex is joined to all others and the complex is a cone).

`validate_cutoff.py` compares the sorted float32 endpoints with unrestricted Ripser on 150 random positive-weight graphs (a ring plus random chords, 5 to 25 vertices) and on 126 saved n=512 graphs from one noisy seed. These are test graph cases, not independent samples. `results/cutoff_validation.json` records exact agreement.

Full diagrams stop at n=512. The larger sweeps compute ranks and bridge events only.
