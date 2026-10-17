# Phase 1: completed computations and checks

All nine Phase 1 tasks are complete. The records below supply the computational evidence for the revision plan. Manuscript restructuring and the integration of these results are tracked separately in that plan.

- **1.1, Example 1:** all 120 cyclic orders were checked at all 14 distinct threshold graphs, including the empty graph. Exact radical arithmetic and outward integer bounds identify ties and certify distance comparisons. Each graph has an admissible order separately; the intersection over thresholds is empty. Direct enumeration independently checks all orientations for all 120 orders at each of the two proof graphs. See [the certificate](../separation_verification.json) and [the exact argument](../arc_support_example/exact_calculation.md).
- **1.2, original clean study:** all 60 samples, both graph rules, were replayed. Sparse event results agree with full Cartesian ranks and independent depth-first bridge traversal. See [raw events](bridgeless/clean_samples.csv) and [summary and Wilson intervals](bridgeless/clean_summary.csv).
- **1.3, full planned sweep:** 15,000 profile/seed/size samples and 30,000 method records through n = 32,768 are complete, using every corresponding released seed. See [raw events](bridgeless/sweep_samples.csv), [summary and intervals](bridgeless/sweep_summary.csv), [validation](bridgeless/validation.json), and [timings](bridgeless/sweep_timings.csv).
- **1.4, Euclidean baseline:** full unrestricted finite H1 diagrams were computed for all 220 distinct observed clouds. Reusing the shared clean cloud gives 240 model/amplitude/seed rows, with each baseline shared across graph methods. See [baseline errors](euclidean/samples.csv), [quantiles](euclidean/summary.csv), and [paired rule comparisons](euclidean/paired_summary.csv).
- **1.5, grid oracle:** 720 sample-specific minima are computed over connected points of the prescribed 24-budget grid; extra diagnostic budgets are excluded. Minimize within each sample before taking group medians or quartiles. See [sample choices](grid_oracle_samples.csv), [36 group summaries](grid_oracle_summary.csv), and [paired selector gaps](selected_vs_oracle.csv). This comparator uses the clean reference.
- **1.6, amplitudes:** [all_amplitudes.csv](all_amplitudes.csv) contains both rules at all six amplitudes, separately by noise model and method, including medians, IQRs and median normalized budgets. Shared clean samples and displacement vectors remain paired.
- **1.7, displayed sweep counts:** [rank_display_counts.json](rank_display_counts.json) confirms 5,000 seed/size combinations, 15,000 profile/seed samples and 30,000 method records in the primary display. The full released rank study retains all 30,780 method records, including the two larger sizes.
- **1.8, connectivity wording:** retain “a simple heuristic”. The bounded primary-source check does not establish that choosing the smallest connecting parameter is standard practice. The checked [Isomap documentation](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.Isomap.html) uses a fixed default; that does not support the stronger claim.
- **1.9, spectral graph:** the [Damrich et al. primary paper](https://papers.nips.cc/paper_files/paper/2024/file/4a32a646254d2e37fc74a38d65796552-Paper-Conference.pdf) specifies unweighted symmetric (OR) kNN graphs in the reported spectral comparisons. Keep the qualified robustness comparison and its parameter dependence. Source locations and the wording decision are in [literature_checks.json](literature_checks.json).

## Connected bridgeless events

Every complete integer-k layer is inserted before testing. The graph must be connected and have no bridges; minimum degree two is not substituted. For every record, k_conn ≤ k_br ≤ k_close. Closure at k_br is a necessary-condition comparison, not a loop guarantee.

At n = 512 in the 1,500-seed-per-profile sweep, the closure percentages at connection and at the first connected bridgeless layer are:

- Trough, OR: 32.7% at connection; 74.1% at the connected bridgeless layer.
- Trough, mutual: 7.1% at connection; 41.2% at the connected bridgeless layer.
- Cosine, OR: 33.2% at connection; 74.1% at the connected bridgeless layer.
- Cosine, mutual: 27.0% at connection; 71.4% at the connected bridgeless layer.
- Uniform, OR: 32.9% at connection; 74.9% at the connected bridgeless layer.
- Uniform, mutual: 28.9% at connection; 71.7% at the connected bridgeless layer.

All rates and pointwise 95% Wilson intervals are in the CSV summaries. Profiles and rules reuse seeds; pooled rows are not independent repetitions. The bridgeless rule leaves a substantial fraction of samples unclosed.

## Noise comparison at the largest amplitude

The following are medians across the same 20 seeds within each group. Euclidean errors are shared across graph methods. A difference of medians is not a median paired improvement; the latter and per-sample win counts are stored separately.

- isotropic, OR_kNN: Euclidean 0.14273; longest-lifetime selection 0.02867; grid oracle 0.00522.
- isotropic, R2_global: Euclidean 0.14273; longest-lifetime selection 0.00851; grid oracle 0.00603.
- isotropic, mutual_kNN: Euclidean 0.14273; longest-lifetime selection 0.01844; grid oracle 0.00545.
- radial, OR_kNN: Euclidean 0.14493; longest-lifetime selection 0.02293; grid oracle 0.00632.
- radial, R2_global: Euclidean 0.14493; longest-lifetime selection 0.01022; grid oracle 0.00805.
- radial, mutual_kNN: Euclidean 0.14493; longest-lifetime selection 0.01465; grid oracle 0.00645.

These comparisons describe the sampled noise models and grid. Neither a long bar nor the connected bridgeless check is a noisy-data theorem.

## Reproduction

Run from final_release with the documented Python environment. Recomputing the connected bridgeless study additionally requires a C++17 compiler (g++); the helper is compiled into the system temporary directory, and no binary is distributed.

```sh
python3 code/verify_separation.py
python3 code/phase1_bridgeless.py --validate --clean --pilot
python3 code/phase1_bridgeless.py --sweep
python3 code/phase1_noise.py
python3 code/phase1_summaries.py
python3 code/phase1_report.py
python3 code/verify_phase1.py
python3 code/record_current.py
python3 code/verify_release.py
```

The bridge and Euclidean scripts resume verified current checkpoints. To recompute an existing Euclidean diagram, remove its output checkpoint first; inputs remain in results/noise/checkpoints. Do not delete or regenerate the original study records. The regular verifier needs no C++ compiler and does not rerun the large persistence study.

Individual Euclidean timings include distances and persistence, excluding checkpoint I/O. Its pilot and full timing records describe separate invocations; per-checkpoint elapsed times cover all 220 diagrams. Bridge timings include sampling, helper startup, graph reconstruction and the event calculation; saved-row writing is excluded. The child RSS field is a cumulative maximum across helper processes, not a per-case allocation; the pilot invocation also includes compiler memory if compilation was needed. Bridge timings were measured while the independent serial Euclidean job was running, so they are feasibility measurements rather than isolated benchmark claims.

The existing Figure PDF/PNG pairs, manuscript and prior experimental data are preserved. No deposition, registration or manuscript release URL is part of this phase.
