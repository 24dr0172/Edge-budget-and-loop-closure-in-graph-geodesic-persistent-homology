# Results

Saved measurements, certificates and check records for the focused paper and additional studies. NPZ files can be loaded with NumPy using `allow_pickle=False`.

## Used by the focused paper

- `retention_table/retention_compact.csv`: the six-row numerical display for Table 3. `retention_counts.csv` contains all 24 primary groups and their counts and Wilson endpoints; `retention_three_layers.csv` contains all 72 rounded group/layer values, and `retention_fixed_k.csv` contains the 48 fixed-k cells. No LaTeX document source is distributed.
- `rank_table.csv`: the 24 primary connectivity-retention cells. `extra_k/rank_sweep/` contains raw consecutive-pair keys, sample sizes, profiles and seeds; the paper uses four sizes. `comparisons/bridgeless/` supplies the matching bridge-free records. `iid_predictions.csv` supplies the Section 8.2 heuristic values.
- `separation_verification.json`, `natural_instance_verification.json`, `small_examples_verification.json`: exact or rational certificates and independent finite checks for Examples 1-2, Remarks 1-2 and Appendix B. All proof arguments are in the paper.
- `clean_closure/`: the 60 independent n=512 clouds and their event diagrams and budget-grid records. These clouds are separate from the large rank sweep. `cutoff_validation.json` records the adaptive cutoff comparison with unrestricted persistence.
- `practical_claims_verification.json`: fixed-k, bridge-free, model and diagram-error values quoted in Sections 7 and 8.
- `retention_extra_layers/`: Figure 3 counts for every primary size, j=0 through 5 (144 rows), plus settings and the source hash.
- `noncircular_locality/`: the Appendix A two-curve study, numerical thresholds and first unsupported chords, with a circle control. Its continuum thresholds are grid estimates.
- `useful_diagram/rational_counterexample.json` and `seven_point_check/`: the Appendix C witness and values. `independent_checks/` records independent geometry and probability calculations. Its additional drawing and diagnostic calculations are retained too.

## Additional studies and shared checks

- `noise/`, `matched_noise/`: noisy coordinates, diagrams, first-nonempty and earliest-longest-bar selectors, summaries and timing records. These selectors are not prescribed by the focused paper.
- `peanut_neck/`: varying-waist study, including a circle control, full event barcodes and one-edge shortcut checks. `peanut_trough/`: 450 trough-density peanut clouds and 1,350 closure records. Both are additional studies.
- `ellipse/`, `ellipse_checks/`: ellipse sweeps, NN-crust comparison and geometry records. `order_checks/` contains CkNN and shuffled-scale controls. `curve_order/` checks the supplied curve order, rather than every possible order.
- `comparisons/`: bridge-free comparisons, Euclidean baselines, grid-oracle values and broader size sweeps; see [their guide](comparisons/README.md).
- `useful_diagram/saved_prefixes.csv`: all stored post-closure diagrams, including those behind the additional error plot. Shared checks verify that plot even when only the focused paper is compared.
- The remaining audit, replay, calibration, runtime and generation records document the preserved broader evidence. Full verification checks every CSV schema and every payload hash. It independently replays selected saved diagrams without regenerating sample cohorts.

## Release records

- `manuscript_record.json`: source, PDF and active-input hashes with paths relative to the external manuscript folder.
- `manuscript_structure.json`: page count, labels, citations and numbering.
- `generation_summary.json`: output counts and the current paper-to-evidence map.
- `file_manifest.json`: SHA-256 hashes of all payload files except itself, the mutable complete verification report and diagnostic run logs.
- `verification_report.json`: complete saved-data and optional manuscript verification. Focused reports must be written outside the release.

The [code guide](../code/README.md) gives generation and recording commands. No separate supplement is part of this release.
