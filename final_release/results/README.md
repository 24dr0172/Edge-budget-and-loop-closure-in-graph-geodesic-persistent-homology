# Results

These files contain the scientific measurements, diagrams and summaries for the accompanying paper. The complete generation command is in [the code guide](../code/README.md). No other experiment folder is required.

- `clean_closure`: the three clean-circle profiles; coordinates, graph events, budget cells, errors, closure selections and sample checkpoints.
- `order_checks`: ordinary CkNN certificates, the shuffled-scale control, its diagrams and the exact precision check.
- `useful_diagram`: full diagrams of the tested closed graph prefixes, their transitions and the seven-point example.
- `extra_k`: additional integer neighbours needed for closure. `candidate_audit` contains the seed-generated size sweep and probability tables.
- `matched_noise`: the radial first-connectivity diagnostic and the before/at-closure comparison, including coordinates and diagrams.
- `noise`: the two noise models, the prescribed budget grid, exact event checks, both selectors and error summaries.
- `ellipse_checks`: seed-generated ellipse coordinates, analytic bounds and the exact connectivity/closure checks.
- `ellipse`: post-closure diagrams, arc-condition checks, nonlocal edge events and NN-crust results.
- `phase1`: the paper's bridge-free, Euclidean-baseline and grid-oracle comparisons. [The guide](phase1/README.md) explains the measurements.

CSV files contain measured rows or summaries. NPZ files contain arrays and persistence diagrams; load them with NumPy using `allow_pickle=False`. Checkpoint files let an interrupted calculation continue. The code guide defines pairing, numerical precision and selector conventions.

`generation_summary.json` lists the generated studies and their sizes. `verification_report.json` records the full release check. `file_manifest.json` lists SHA-256 hashes. When manuscript comparison is requested, the manuscript records identify the external source and PDF by their hashes. `generation_timing.json` and `runtime_benchmark.json` describe this run's execution times; clock times are not fixed scientific targets.

Counts distinguish distinct samples, model-indexed rows and method records. A repeated clean row across the two noise models is not an additional independent sample.
