# Edge budget and loop closure

This release contains the data, code and figures for *Edge budget and loop closure in graph-geodesic persistent homology*. The paper asks when a connected graph still misses a sampled curve's loop, and what happens when more edges are added.

Run all commands below from `final_release/`. Full verification also needs `combined.tex` and its compiled `combined.pdf` in the parent folder.

## Set up and check the release

Use Python 3.12 and install the pinned packages:

```sh
python3 -m pip install -r code/requirements.txt
python3 code/verify_release.py
```

The verifier checks file hashes, manuscript tables, saved diagrams, error calculations, sample coordinates, selection rules and graph events. It writes `results/verification_report.json`. It replays selected calculations rather than rerunning every persistence calculation.

The supported environment is recorded in `config.json` and `code/requirements.txt`. Some measurement records list pandas 3.0.5; the supported environment for rebuilding summaries and checking the release uses pandas 2.2.3. See `code/README.md` for the numerical conventions.

## Rebuild summaries and figures

These commands use the saved measurements:

```sh
python3 code/reproduce_results.py
python3 code/build_figures.py
python3 code/record_current.py
python3 code/verify_release.py
```

The first command rebuilds the clean, noise, ellipse, Euclidean-baseline and grid-oracle summaries. It also exports all six rows of Table 3, including the connected bridge-free comparison, and updates the comparison guide.

To refresh the seven numerical table bodies in the manuscript, run `python3 code/update_manuscript.py`. Compile `../combined.tex` in place before running `record_current.py` and the verifier. The recording command updates the manuscript identities, structure record, comparison record and file manifest. Compilation is a separate step.

## Rerun the studies and checks

The following sequence covers the noise and ellipse measurements, rank sweep, connected bridge-free comparison, Euclidean baseline and supporting checks:

```sh
python3 code/run_studies.py --workers 6
python3 code/check_knn_miss_candidate.py
python3 code/phase1_bridgeless.py --validate --clean --sweep
python3 code/phase1_noise.py
python3 code/reproduce_results.py
python3 code/verify_separation.py
python3 code/verify_arc_support.py
python3 code/check_cknn_precision.py
python3 code/check_rank_pairing.py
python3 code/validate_cutoff.py
python3 code/build_figures.py
python3 code/record_current.py
python3 code/verify_release.py
```

The study scripts reuse completed checkpoints. The bridge calculation needs g++ with C++17 support; checking the supplied release does not need a compiler. Clean-circle and matched radial-noise trials can be rerun separately with `regenerate_trials.py`. The seed loops and instructions for starting fresh in a working copy are in `code/README.md`.

## Find the evidence

- `results/clean_closure/`: clean-circle connectivity, closure and error measurements.
- `results/order_checks/`: edge-order checks, the CkNN exception and the shuffled control.
- `results/useful_diagram/`: saved diagrams after closure and the seven-point capped-error example.
- `results/extra_k/`: paired extra-neighbor counts and the sample-size sweep.
- `results/noise/`: both noise models, all amplitudes and both selection rules.
- `results/matched_noise/`: radial-noise checks at connectivity and around latent closure.
- `results/ellipse/` and `results/ellipse_checks/`: ellipse measurements, coordinates and geometric checks.
- `results/phase1/`: bridge-free comparisons, Euclidean baselines, grid oracles and all-amplitude summaries.
- `results/arc_support_example/`: the exact calculation for Example 1.
- `manuscript_figures/`: six PDF assets, each with a matching 300-dpi PNG.

The [results guide](results/README.md) explains the evidence and its limits. The [code guide](code/README.md) describes individual commands, precision and runtimes.

## Read the results

The clean study uses 20 seed streams and three chosen sampling profiles, giving 60 distinct samples; profiles share draws within each seed. The noise study uses 20 separate streams. Within each noise model, the same displacement vectors are scaled across amplitudes; the two models share the zero-noise observations.

The clean-data guarantees require the conditions stated in the paper. The two selection rules coincide under Corollary 2's assumptions. Their noise comparisons are empirical. At the largest tested noise amplitude, longest-lifetime selection has median normalized diagram error 0.0085–0.029 across the six groups, compared with about 0.32 for first-nonempty selection.

Full diagrams were computed up to n=512. The larger study, reaching n=1,000,000, computes ranks only. Table D.2 displays sizes through n=32,768; the release supplies the two larger sizes and their uncertainty intervals.
