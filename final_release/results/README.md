# Results

These are saved measurements and checks. Load NPZ files with NumPy using `allow_pickle=False`. The [release README](../README.md#tables-and-figures) maps the manuscript's tables and figures to their sources.

## Paper evidence

- `retention_table/` and `extra_k/rank_sweep/`: the 24 profile/method/size cells in Table 3, their counts, intervals and raw ranks.
- `clean_closure/` and `post_closure_rate/`: clean-circle diagrams and Figure 3's 320 records, including one omitted radius-doubling miss. `exact_optimum/` records the closure error and empty-diagram cost for 180 sample–order pairs.
- `peanut_neck/`: 100 clouds paired over five waist parameters and 300 method records, with 15 saved closure groups. Table 4 aggregates the 12 noncircular groups into four waist rows; Figure 4 also shows the circle control. `first_nonempty/` records the prefix search and radius-doubling outcomes on those same clouds.
- `radius_enlargement/`: 400,000 saved circle gap pairs and 12 rate rows. `radius_enlargement_extra/`: independent million-cloud checks at $\alpha=2$ for $n=4,096$ and $32,768$, with replayable batch counts. `smooth_curve_enlargement/`: the ellipse and peanut chord-event proxies, batch counts and full-diagram spot checks for Table 5. Proxy events are identified as such in the manuscript.
- `rational_certificates/`, `seven_point_check/`, `useful_diagram/rational_counterexample.json`, `natural_instance_verification.json`, and `noncircular_locality/`: finite examples, exact or interval checks, and numerical locality estimates.

## Additional studies

`winding_check/`, `retention_extra_layers/`, `certified_ellipse_window/`, `noise/`, `matched_noise/`, `peanut_trough/`, `ellipse/`, `ellipse_checks/`, `order_checks/` and [comparisons](comparisons/README.md) contain additional investigations and measurements. [additional_figures](additional_figures/README.md) contains the corresponding plots.

## Records

`manuscript_record.json` hashes the external manuscript source, PDF and its three figure files. `manuscript_structure.json` records pages, result labels, citations and environment counts. `generation_summary.json` gives saved-data counts and the evidence map. `file_manifest.json` hashes release files, excluding itself, the verification report and diagnostic logs. `verification_report.json` describes the most recent complete saved-data verification; a paper-only report is written outside the release. See the [tools guide](../tools/README.md).
