# Code guide

Run commands from `final_release/` with Python 3.12 and the versions in `requirements.txt`. The complete study and numerical conventions are described in `../README.md`.

- `reproduce_results.py`: rebuild tables and audits from saved rows, without graph or persistence recomputation. `--rank-table-only` extracts the four main-text n=512 rank rows.
- `build_figures.py`: export the clean-grid, noise, extra-k and schematic panel assets.
- `verify_release.py`: check manuscript identity, release hashes and scientific results, including `verify_studies.py` and `audit_rank_scaling.py`; write `results/verification_report.json`.
- `record_current.py`: refresh source/PDF identity and the file manifest after successful compilation and numerical checks.
- `run_studies.py`: compute or resume all two-model noise and post-closure ellipse checkpoints. `--aggregate-only` rebuilds aggregate rows; `--validate-only` runs independent helper checks.
- `summarize_studies.py`: rebuild the noise and ellipse summaries and their figures.
- `audit_numbers.py`: reproduce profile closure counts, exact-closure medians, rank replay, post-closure extrema and paired-profile calibration.
- `audit_rank_scaling.py`: verify all raw size-sweep closure fractions and pointwise Wilson intervals, and explain the paired large-size indicators.
- `verify_arc_support.py`: deterministic small-graph checks, independent triangle-boundary calculations over F2/F3/F5, and the 180 ellipse closure certificates. Its saved report is `results/arc_support_verification.json`.
- `verify_separation.py`: verify arc support, Betti numbers and all 1,024 curve-order orientations, and all 120 cyclic orders across the exact threshold filtration of Example 1.
- `validate_cutoff.py`: compare adaptive and unrestricted H1 persistence.
- `update_manuscript.py`: refresh only numerical table bodies from the saved summaries; surrounding prose and proofs stay in the manuscript source.
- `paper1_methods.py`, `finite_diagram.py`: scientific kernel and independent finite-complex reduction/bottleneck matcher.

## Optional individual trial regeneration

```sh
python3 code/regenerate_trials.py --study clean --seed 20310101 --output ../scratch
python3 code/regenerate_trials.py --study useful --seed 20310101 --output ../scratch
python3 code/regenerate_trials.py --study orders --seed 20310101 --output ../scratch
python3 code/regenerate_trials.py --study ellipse --seed 20310101 --output ../scratch
python3 code/regenerate_trials.py --study noise --seed 20320101 --output ../scratch
```

Repeat clean/useful/orders/ellipse for seeds 20310101--20310120 and noise for 20320101--20320120. Scratch output stays outside the delivered package. `export_noise_coordinates.py --output ../noise_coordinates.npz` reconstructs the radial-coordinate arrays from the saved configuration without computing persistence.

`check_knn_miss_candidate.py` computes the rank-only size sweep. `--verify-only` checks its summaries, probability-model values and selected seed replays. To regenerate every sample, move both `results/extra_k/candidate_audit/samples.csv` and `summary.csv` to scratch first. Profiles share uniform draws and phase at each size/seed, so profile closure rates are paired. The two largest sizes have only 100 and 30 replicates per profile; Their Wilson intervals remain in results/extra_k/candidate_audit/summary.csv; the manuscript displays sizes through 32,768 in Table D.2.

## Precision, paired ranks and timing

- `check_cknn_precision.py` checks the failing CkNN insertion using exact rational edge scores and outward integer-square-root cost bounds. Its certificate is `results/order_checks/precision_check.json`; the ideal-circle cost is also reevaluated at 60 and 100 digits.
- `check_rank_pairing.py` replays all 630 uniform/trough OR sample pairs at the three largest sizes and checks the critical inverse-CDF branches; its report is `results/rank_pairing_replay.json`.
- `benchmark_diagrams.py --repeats 3` runs fresh serial 24-budget, three-method blocks at n=512. It writes `results/runtime_benchmark.json`. These timings exclude sampling, references, bottleneck errors, certificates and file I/O; they are not full-study generation times. The full-diagram studies stop at n=512; only the rank study reaches n=1,000,000.

The selector summary includes `q1_e_PH` and `q3_e_PH`, using linear interpolation of the selected errors across the 20 seeds. `update_manuscript.py` writes the median and interquartile range into the selected-rule comparison table.

## Phase 1 computations

All Phase 1 outputs and interpretation are indexed by `results/phase1/README.md`.

- `phase1_bridgeless.py --validate --clean --pilot` independently validates the incremental bridge helper and replays every original clean sample; `--sweep` completes all paired sweep seeds through n = 32,768. Recomputing this study requires g++ with C++17 support. `bridgeless_circle.cpp` is compiled into a temporary directory; no binary is distributed.
- `phase1_noise.py` computes or resumes full unrestricted Euclidean Rips diagrams for the 220 distinct observed clouds and saves paired comparisons to the existing rules.
- `phase1_summaries.py` reproduces sample-wise connected-grid oracles, all-amplitude medians/IQRs, display counts and the recorded literature decisions.
- `phase1_report.py` regenerates the current completion record and results guide from those outputs.
- `verify_phase1.py` checks all new saved evidence, including exact common-order separation, bridge-event replays, Euclidean scores, oracle choices and paired summaries. It is also included in `verify_release.py`. This saved-evidence verification needs no C++ compiler and does not recompute the large persistence study.

The same 24-budget grid is used for both rules and the oracle. Additional exact-event diagnostic budgets are excluded from the oracle. Euclidean baselines share clean samples at zero amplitude and are reused across graph methods.

## Numerical conventions and representative timings

Main trials select SciPy Dijkstra explicitly; the separate arc-support verifier uses SciPy's automatic shortest-path selection. Edge costs, scores and shortest paths use float64; Ripser endpoints use float32. Every finite returned interval with death strictly greater than birth is retained, with no added lifetime threshold. The observed minimum positive lifetime is 1.210719347e-8 raw units; it is not a universal resolution.

Raw comparison tolerances are 5e-8 for post-closure grid/stability checks, 1e-7 for selected-budget verification and 1e-6 for ellipse interval comparisons. They do not change strict endpoint positivity. The maximum clean birth deviation is 3.54e-9 and the largest clean gap is 0.08955721.

The representative benchmark uses an AMD Ryzen 7 7445HS CPU, one serial process, and OpenBLAS/OMP/MKL thread limits of one. Each block constructs all edge orders and evaluates 24 budgets for each of the three methods, including all-pairs Dijkstra and the full finite H1 diagram at connected budgets. Noise uses the validated adaptive cutoff; clean-circle and ellipse cases use unrestricted Ripser. Sampling, intrinsic references, bottleneck errors, certificates, figures and file I/O are excluded. These are representative block times, not total study-generation times.

The following seconds are component-wise medians over three repetitions, followed by the range of total time. Full software versions and per-call timings are in `results/runtime_benchmark.json`.

- Clean circle, uniform: orders 0.04; Dijkstra 4.04; diagrams 110.61; total 114.66, range 114.66–116.12.
- Radial noise, delta=0.01: orders 0.05; Dijkstra 4.31; diagrams 73.91; total 78.27, range 77.43–84.32.
- Disk noise, delta=0.01: orders 0.05; Dijkstra 4.37; diagrams 74.94; total 79.36, range 75.75–80.95.
- Ellipse, uniform: orders 0.04; Dijkstra 4.80; diagrams 111.91; total 116.76, range 115.90–122.42.

Full-diagram studies use n=512. The full released rank-only study reaches n=1,000,000; its primary manuscript display stops at 32,768. All larger-size raw records, summaries and Wilson intervals are retained. The field `condition_5c` in the ellipse data identifies the sufficient ordinary-error test in the arc-support theorem, currently Theorem 3(c); the column name is a machine-readable identifier.

## Manuscript table consistency

`update_manuscript.py` refreshes seven numerical table bodies in the existing `combined.tex`: the connected bridge-free comparison; selected noise errors and IQRs; Euclidean and grid-oracle errors with paired wins; all-amplitude selector medians; noise screens; ellipse results; and the primary rank display. `verify_manuscript_tables.py`, called by `verify_release.py`, reads the saved summaries independently and checks the new displayed values. It neither changes the manuscript nor recomputes persistence. The full release verifier separately validates those summaries against sample-level evidence.

The table exporter preserves the repeated headers and footers of Table 6's longtable and the shared clean Euclidean baseline row in Table 7. The manuscript-table verifier checks that baseline against the saved Euclidean summary as well as the existing row-by-row comparisons. Regenerating the current seven table bodies leaves the source unchanged.
