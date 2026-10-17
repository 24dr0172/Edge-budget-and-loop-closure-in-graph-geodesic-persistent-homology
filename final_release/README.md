# Edge budget and loop closure

This package accompanies `../combined.tex` and `../combined.pdf`. All data, code, configuration and figure assets used by the paper are contained here. Paths below are relative to this folder.

## Reproduce and verify

Use Python 3.12 and `code/requirements.txt`. From this folder:

```sh
python3 code/reproduce_results.py
python3 code/build_figures.py
python3 code/record_current.py
python3 code/verify_release.py
```

`python3 code/update_manuscript.py` refreshes the seven numerical manuscript table bodies from the saved summaries, including the connected bridge-free comparison and the Euclidean, grid-oracle and all-amplitude noise tables. It leaves prose and proofs untouched.

The first command rebuilds numerical summaries from saved measurements, including the noise and ellipse studies and the raw-data check of the large-sample closure rates. The second exports the figure assets. After editing the manuscript, compile its existing source in place before recording its identity. `record_current.py` records the current manuscript identity and file hashes; it does not compile or validate a PDF. The verifier checks those identities, saved diagrams and errors, seeded coordinates and displacement draws, grid selectors, NN-crust recovery, rank events, closure fractions and Wilson intervals. Its report is `results/verification_report.json`. It does not rerun every shortest-path and persistence calculation.

## Evidence

- `results/clean_closure/`, `results/order_checks/`, `results/useful_diagram/`: clean-circle events, error curves, order certificates, controls and post-closure diagrams.
- `results/noise/`: two noise models, six amplitudes, three methods, 20 seeds and the 24-budget grid; per-sample CSV/NPZ checkpoints include every computed diagram and the coordinates and draws needed to reconstruct it.
- `results/matched_noise/`: first-connectivity and before/at-closure radial diagnostics on the same noise seeds.
- `results/ellipse/`, `results/ellipse_checks/`: the full 3,516-graph sweep, 60 coordinate sets, closure diagrams, local/nonlocal certificates, cycle-specific Euclidean death and NN-crust comparisons.
- `results/extra_k/`: the 120 paired n=512 rank measurements and the rank-only size sweep. `results/rank_scaling_audit.json` checks the displayed large-size counts, paired indicators and Wilson intervals against the raw records.
- `results/arc_support_example/exact_calculation.md`, `results/separation_verification.json`: exact six-point calculations, the 1,024 curve-order orientations, and all 120 orders throughout the threshold filtration of Example 1.
- `results/number_audit.json`, `results/helper_validation.json`, `results/cutoff_validation.json`: numerical audits and independent helper checks.
- `results/order_checks/precision_check.json`: exact CkNN ordering and rigorous chord-sum bounds for the failing insertion.
- `results/rank_pairing_replay.json`: paired uniform/trough OR replay at the three largest sizes, including the two cancelling changes at n=32,768.
- `results/runtime_benchmark.json`: fresh serial 24-budget timings at n=512, repeated three times for four representative cases.
- `manuscript_figures/`: the assets used by the manuscript; `config.json` records seeds, conventions and package versions.

## Full calculation

```sh
python3 code/run_studies.py --workers 6
python3 code/check_knn_miss_candidate.py
python3 code/verify_separation.py
python3 code/validate_cutoff.py
```

`run_studies.py` resumes per-sample checkpoints. `--aggregate-only` rebuilds aggregate rows without recomputing graph metrics; `--validate-only` performs the small independent helper comparisons. The rank command also resumes stored samples; `--verify-only` checks stored summaries and selected seed replays. `regenerate_trials.py` supports individual clean/order/closure trials with explicit scratch output. See `code/README.md` for details.

Clean/ellipse seeds are 20310101--20310120; noise seeds are the separate streams 20320101--20320120. Noise amplitudes are 0, 0.0005, 0.001, 0.002, 0.005 and 0.01. Each model scales a common unit displacement vector across amplitudes. The clean observations are shared across models, giving 60 distinct clean method/seed cases, each evaluated by two selectors. At delta=0.01, longest-lifetime selection has median normalized diagram error 0.0085--0.029 across the six model/method groups, versus about 0.32 for the first-nonempty rule. The selected-error summary includes the first and third quartiles, with linear interpolation across 20 seeds. These are finite empirical results.

Edges and shortest paths use float64; Ripser endpoints use float32. Every finite returned interval with death strictly greater than birth is retained, with no positive-lifetime cutoff. The retrospective selector uses raw lifetime tolerance 1e-7 to choose the earliest grid maximizer. Support uses 1e-10, ellipse monotonicity and sufficient-inequality checks use 1e-6, and same-vertex stability checks use raw additive tolerance 5e-8.

The adaptive H1 cutoff is exact for a positive weighted graph metric: above the largest generating edge, every Rips edge is homologous to its shortest graph path by an endpoint fan. Later H1 maps are surjective; once no class survives beyond that edge, the full positive H1 barcode has been found. The fallback is the minimum vertex eccentricity, where a universal vertex makes the complex a cone. The validation report compares 150 random weighted graphs and 126 saved n=512 diagrams with unrestricted Ripser.

Full persistence diagrams are computed only up to n=512. The size sweep up to n=1,000,000 uses local rank counts and does not compute diagrams. `code/benchmark_diagrams.py --repeats 3` measures the graph-order, shortest-path and diagram parts on representative samples; the results are not total study generation times.

## Phase 1 computations

The computations in the current revision plan are complete; manuscript integration belongs to the later phases. The manuscript source and PDF are unchanged by this phase. `results/phase1/README.md` indexes the full connected-bridgeless sweep through n = 32,768, Euclidean Rips baselines, grid oracles, amplitude summaries, exact separation certificate and literature decisions. Existing experimental data and PDF/PNG figure pairs are preserved.

Run `python3 code/verify_phase1.py` for the added saved-evidence checks, or `python3 code/verify_release.py` for the complete package. Recomputing the sparse bridge sweep additionally requires a C++17 compiler; ordinary verification does not. All work remains local.

## Manuscript integration

Table 3 uses the paired connected bridge-free sweep and pointwise Wilson intervals. Tables 5–7 retain selector IQRs and show the Euclidean baseline, the connected-grid oracle, paired wins and both selectors at every noise amplitude. Example 1 includes the two-threshold common-order proof. The general and circle statements retain their separate survival and monotonicity assumptions. The current structural/proof audit is in `results/manuscript_restructure.json`; the integrated verifier also checks all new table values against the released summaries.

The abstract gives the plain, number-free message. The introduction reports the 1,500-seed-per-profile simulation rates, distinguishes them from theorems, and states the circle scope of the error comparison. Its four contributions and novelty table follow the audit, exact circle events, circle error comparison and general criterion. The exact common-order obstruction is stated with its individual-threshold qualification. Formal statements, complete proofs and scientific data are unchanged. Contribution qualifications, clean-data scope, the cost statement and numerical-convention pointers agree with the manuscript.

The proof audit covers every proof environment, the closure proposition, the adaptive cutoff argument, the exact CkNN certificate, the circle small-gap derivation and the complete probability derivation. All 49 original displayed formulas are retained. The current PDF has 39 pages; references begin partway through page 37. Paired seeds, scientific inputs and all PDF/PNG figure pairs are unchanged.

The ten-point consistency review is recorded in `results/manuscript_restructure.json`. Appendix D.2 keeps the seven-point example's coordinates, bar endpoints, direct matching and diagonal matching together, with both normalized costs checked against the saved data. Numerical precision and tolerance conventions have one home in D.2. Procedure 1 reports the representative runtime and full-diagram size limit. The table exporter preserves longtable headers and the shared clean Euclidean baseline row when regenerating the manuscript; its output is unchanged when the current tables are already up to date. Current identity and verification results are in `results/manuscript_revision.json` and `results/verification_report.json`.
