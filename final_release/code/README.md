# Code guide

Use Python 3.12 and [requirements.txt](requirements.txt); run commands from the release root. The numerical pipeline is independent of the separately supplied manuscript. The release README maps the saved data and figures to the manuscript. Study settings record the generation-source checksum and, when distinct, the verification-source checksum.

## Verify saved results

```sh
python3 code/verify_release.py --paper-only --manuscript-dir /path/to/paper --report /tmp/edge-budget-paper-check.json
python3 code/verify_release.py --manuscript-dir /path/to/paper --report /tmp/edge-budget-complete-check.json
```

The first checks the paper, file hashes and selected evidence; the second includes additional studies and replays. To skip the manuscript comparison, use `--data-only` instead. The paper-only run also performs additional evidence checks. Neither command regenerates all study cohorts. [tools/README.md](../tools/README.md) explains how to update the records after a change.

## Rebuild from saved data

- `build_retention_table.py` exports circle rank counts and pointwise Wilson intervals.
- `check_tied_ranks.py` checks the spacing tests and evaluates the approximation in Theorem 4 at the Table 1 sample sizes.
- `python3 code/radial_noise_knn.py --output /tmp/radial_noise_knn.csv` regenerates the 400 paired point sets for Table 2.
- `post_closure_rate.py --plot-only` redraws Figure 3 from saved records; `--verify` checks them and replays selected graphs.
- `build_peanut_comparison.py` reconstructs selected edges for the 300 saved closure graphs, checks U_G, refreshes the summary and draws Figure 4. Its JSON also includes winding-based checks as additional evidence. Add `--manuscript-dir /path/to/paper` to copy the rebuilt PNG to that paper's `figures/` directory. Regenerated image bytes can differ across plotting environments; rerecord the manuscript only after compiling the intended source and reviewing its PDF.
- `radius_enlargement_check.py --verify`, `radius_enlargement_extra.py --verify` and `smooth_curve_enlargement.py --verify` check the saved circle and smooth-curve enlargement studies.
- `verify_first_nonempty.py --replay` checks the paired peanut stopping study and selected diagrams. `first_nonempty_rule.first_nonempty_circle` implements the circle bisection rule of Corollary 2; it assumes its hypotheses and does not certify them from unordered input.
- `verify_intrinsic_circle.py` checks the exact winding-triple reference formula by rational triangle-boundary reduction over F2, F3 and F5; the saved summary is `results/intrinsic_circle_verification.json`.
- `rational_certificates.py`, `check_seven_point_example.py`, `verify_natural_instance.py`, `verify_noncircular_locality.py` and `verify_arc_support.py` cover the paper's finite examples and supporting calculations.
- `retention_extra_layers.py` rebuilds the additional extra-layer plot in `results/additional_figures/`.

Some of these commands write generated results. After rebuilding, follow the [recording and verification steps](../tools/README.md); generating a figure does not itself update manuscript hashes or its compiled PDF.

## Full regeneration

```sh
python3 code/reproduce.py --check
python3 code/reproduce.py --output /path/to/new-run --workers 6
```

The full run starts from an empty output directory, generates the paper and additional studies, and can take hours. It retains the supplied Figure 4 PNG and PDF while rebuilding its underlying data, so the external manuscript comparison uses the same frozen figure. Use `--resume` with unchanged code to continue an interrupted run. The bridge-free generator requires a C++17 compiler. The recorded `code/config.json` describes sampling, seeds, grids and tolerances. Numerical diagrams keep every returned positive bar; shortest paths use float64 and Ripser endpoints use float32.

Clean-circle checks use tolerances of `1e-6` for normalized diagram error and `1e-10` for arc support in raw distance units. Peanut closure checks use `5e-8` for saved endpoints in raw distance units and `5e-12` for the edgewise ratio. These tolerances do not impose a lifetime cutoff.
