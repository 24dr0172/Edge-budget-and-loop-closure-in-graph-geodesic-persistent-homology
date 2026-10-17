# Code guide

Run commands from `final_release/` with Python 3.12 and `code/requirements.txt`. Full verification requires `../combined.tex` and its compiled `../combined.pdf`. The [release README](../README.md) gives the complete study workflow.

## Rebuild and verify

- `reproduce_results.py` rebuilds summaries from saved measurements, including Euclidean comparisons and grid oracles. `--rank-table-only` exports the six Table 3 rows, with connection, bridge-free and +3 rates and pointwise 95% Wilson intervals.
- `build_figures.py` exports all six PDF/PNG asset pairs.
- `update_manuscript.py` refreshes seven numerical table bodies in `../combined.tex`. It preserves surrounding text, Table 6's repeated headers and Table 7's shared clean Euclidean baseline.
- `record_current.py` records the source/PDF hashes, current manuscript structure, comparison identities and release manifest. Compile the manuscript before recording it.
- `verify_release.py` checks the manuscript, file hashes and scientific results, and writes `results/verification_report.json`.
- `verify_manuscript_tables.py` compares the displayed values and the six-row Table 3 export with saved summaries. It is included in the release verifier.

## Study calculations

- `run_studies.py` computes or resumes the two-model noise and ellipse studies. `--aggregate-only` combines saved checkpoints; `--validate-only` runs small independent checks.
- `check_knn_miss_candidate.py` computes or resumes the rank-only size sweep. `--verify-only` checks summaries, probability calculations and selected seeds.
- `phase1_bridgeless.py --validate --clean --sweep` validates the bridge calculation, replays the 60 clean samples and runs the 30,000-record comparison through n=32,768. It needs g++ with C++17 support. The helper is compiled in a temporary folder.
- `phase1_noise.py` computes or resumes the 220 distinct Euclidean diagrams. `--summarize-only` rebuilds comparisons from those saved diagrams.
- `phase1_summaries.py` computes the best connected grid budget separately for each sample, then summarizes those errors. It also rebuilds both rules' all-amplitude summaries and the displayed sample counts.
- `phase1_report.py` builds the [comparison guide](../results/phase1/README.md) and its current study record.
- `verify_phase1.py` checks these comparisons and is included in the release verifier.

The grid oracle uses the clean reference to find the smallest error on the prescribed 24-budget grid. It excludes extra diagnostic budgets and serves as a benchmark. Euclidean baselines are shared across graph methods, and clean observations are shared across noise models.

The sample-size sweep uses paired profiles: the same uniform draws and phase are reused at each size and seed. Table D.2 shows four sizes through n=32,768. Results at n=262,144 and n=1,000,000 use 100 and 30 samples per profile. Their pointwise Wilson intervals are in `results/extra_k/candidate_audit/summary.csv`.

## Rerun individual trials or start fresh

These commands write fresh clean-circle, order-check and closure measurements outside the release:

```sh
for seed in $(seq 20310101 20310120)
do
    python3 code/regenerate_trials.py --study clean --seed "$seed" --output ../scratch
    python3 code/regenerate_trials.py --study useful --seed "$seed" --output ../scratch
    python3 code/regenerate_trials.py --study orders --seed "$seed" --output ../scratch
    python3 code/regenerate_trials.py --study ellipse --seed "$seed" --output ../scratch
done
for seed in $(seq 20320101 20320120)
do
    python3 code/regenerate_trials.py --study noise --seed "$seed" --output ../scratch
done
```

The ellipse option replays the connectivity and closure events; `run_studies.py` covers the full ellipse sweep. The noise option replays the matched radial first-connectivity trials; `run_studies.py` covers both models and the full grid.

The main study commands resume saved outputs. To start fresh, work in a copy of the release and move the outputs to a scratch folder first:

- Noise and ellipse: their `checkpoints/` folders under `results/noise/` and `results/ellipse/`. Keep the separate `ellipse_checks/coordinates.npz` input.
- Rank sweep: both `samples.csv` and `summary.csv` in `results/extra_k/candidate_audit/`.
- Bridge comparison: move both `sweep_samples.csv` and `sweep_timings.csv` from `results/phase1/bridgeless/`. The clean comparison is recomputed whenever `--clean` is used.
- Euclidean baseline: `results/phase1/euclidean/checkpoints/`, after regenerating its noise inputs.

Then run the study sequence in the release README. Keep the supplied input coordinates and configuration. `export_noise_coordinates.py --output ../noise_coordinates.npz` reconstructs the matched radial inputs without computing diagrams.

## Supporting checks

- `audit_numbers.py` checks closure counts, error extremes and the small clean study's paired-seed calibration.
- `audit_rank_scaling.py` checks every rank-sweep rate and Wilson interval.
- `check_rank_pairing.py` replays 630 uniform/trough pairs at the three largest sizes and checks how the sampling transformation changes the relevant counts.
- `check_cknn_precision.py` confirms the CkNN exception using exact edge-order comparisons and rigorous bounds on chord sums, with 60- and 100-digit checks.
- `verify_arc_support.py` checks small graphs over F2, selected cases over F3/F5, and all 180 ellipse closure graphs.
- `verify_separation.py` checks Example 1, including 1,024 curve-order orientations and all 120 cyclic orders across 14 threshold graphs.
- `validate_cutoff.py` compares the adaptive H1 calculation with unrestricted persistence.
- `paper1_methods.py` contains sampling, graph and persistence functions. `finite_diagram.py` independently checks the seven-point example using matrix reduction and bottleneck matching.

## Environment and numerical conventions

`requirements.txt` pins the supported environment for summaries and verification, including pandas 2.2.3. The matched-noise configuration and useful-diagram execution record list pandas 3.0.5 with their saved measurements. Those records describe measurement generation; use the pinned requirements for the release commands.

Main trials use SciPy Dijkstra. The arc-support checker lets SciPy choose the shortest-path algorithm. Edge costs, scores and shortest paths use float64; Ripser endpoints use float32.

Every finite interval returned by Ripser with death strictly greater than birth is retained. No positive-lifetime cutoff is applied to those diagrams. The separate seven-point matrix-reduction helper uses a 1e-12 tolerance and reports any smaller positive lifetime separately; no positive interval is discarded in the saved seven-point example. This helper does not supply the main study diagrams.

Raw comparison tolerances are 5e-8 for the clean grid and stability checks, 1e-7 for selected-budget verification, and 1e-6 for ellipse comparisons. Arc support uses 1e-10. Longest-lifetime selection uses raw tolerance 1e-7 to choose the earliest grid maximum. These comparison tolerances do not change the Ripser positivity rule.

The shortest stored positive lifetime is 1.210719347e-8 in raw units. The largest clean birth deviation is 3.54e-9 and the largest clean gap is 0.08955721. These are observed values. Error quartiles use linear interpolation across the 20 seeds. The ellipse field `condition_5c` records the sufficient ordinary-error test in Theorem 3(c).

Noise diagrams use the exact adaptive H1 cutoff justified in Appendix D.2. Clean-circle and ellipse diagrams use unrestricted Ripser. The independent validation compares 150 random weighted graphs and 126 saved n=512 diagrams.

## Runtime and scale

Full diagrams were computed up to n=512; only rank calculations reach n=1,000,000. Run `python3 code/benchmark_diagrams.py --repeats 3` to repeat the representative benchmark.

The saved benchmark uses an AMD Ryzen 7 7445HS CPU, one serial process and one thread for OpenBLAS/OMP/MKL. Each block builds all three edge orders and evaluates 24 budgets per method, including shortest paths and diagrams at connected budgets. Sampling, reference diagrams, error matching, certificates, figures and file I/O are excluded.

Component medians in seconds over three repetitions are:

- Clean circle: orders 0.04; shortest paths 4.04; diagrams 110.61; total 114.66 (range 114.66–116.12).
- Radial noise, delta=0.01: orders 0.05; shortest paths 4.31; diagrams 73.91; total 78.27 (77.43–84.32).
- Disk noise, delta=0.01: orders 0.05; shortest paths 4.37; diagrams 74.94; total 79.36 (75.75–80.95).
- Ellipse: orders 0.04; shortest paths 4.80; diagrams 111.91; total 116.76 (115.90–122.42).

These are representative block times. Full measurements, software versions and timing exclusions are in `results/runtime_benchmark.json`. Bridge and Euclidean timing definitions are in the comparison guide.
