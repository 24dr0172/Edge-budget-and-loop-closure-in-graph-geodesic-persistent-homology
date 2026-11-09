"""Write the release's generated guides; the root README is edited separately."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
GUIDES={
'code': r'''# Code guide

Use Python 3.12 and [requirements.txt](requirements.txt); run commands from the release root. The numerical pipeline is independent of the separately supplied manuscript. The release README maps the saved data and figures to the manuscript. Study settings record the generation-source checksum and, when distinct, the verification-source checksum.

## Verify saved results

```sh
python3 code/verify_release.py --paper-only --manuscript-dir /path/to/paper --report /tmp/edge-budget-paper-check.json
python3 code/verify_release.py --manuscript-dir /path/to/paper --report /tmp/edge-budget-complete-check.json
```

The first checks the paper, file hashes and selected evidence; the second includes additional studies and replays. To skip the manuscript comparison, use `--data-only` instead. The paper-only run also performs additional evidence checks. Neither command regenerates all study cohorts. [tools/README.md](../tools/README.md) explains how to update the records after a change.

## Rebuild from saved data

- `build_retention_table.py` exports circle rank counts and pointwise Wilson intervals.
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
''',
'results': r'''# Results

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
''',
'manuscript_figures': r'''# External manuscript figures

The manuscript draws Figure 1 directly in its LaTeX source. These three files are its external PNGs and match the manuscript's `figures/` files byte for byte:

- Figure 2: `arc_support_fan.png`, the triangle-fan drawing.
- Figure 3: `post_closure_rate.png`, from [saved rate records](../results/post_closure_rate/records.csv).
- Figure 4: `peanut_comparison.png`, the peanut comparison, including the edgewise certificate bound where applicable, from [saved diagrams and bounds](../results/peanut_neck/closure_error_bound.json).

[PDF exports](../results/figure_exports/README.md) are optional; the saved Figure 4 PDF is a raster export of the manuscript PNG. Further plots are in [additional_figures](../results/additional_figures/README.md). See the [code guide](../code/README.md) for rebuilding.
''',
'results/figure_exports': r'''# PDF figure exports

`arc_support_fan.pdf` and `post_closure_rate.pdf` are optional exports of Figures 2 and 3. `peanut_comparison.pdf` is a raster PDF export of the manuscript's Figure 4 PNG. The manuscript reads the three PNGs in [manuscript_figures](../../manuscript_figures/README.md); Figure 1 is drawn inline. Further plot exports are in [additional_figures](../additional_figures/README.md).
''',
'results/additional_figures': r'''# Additional-study figures

These plots describe additional investigations. `peanut_winding_comparison` shows a winding-bound comparison, `retention_extra_layers` shows retention after added neighbor layers, and `closure_schematic_comparison` illustrates connectivity and closure. See the [results guide](../README.md) and [code guide](../../code/README.md) for saved data and scripts. The manuscript's external PNGs are in [manuscript_figures](../../manuscript_figures/README.md).
''',
'tools': r'''# Release tools

Run these commands from the release root. The manuscript is supplied separately.

- `python3 tools/write_guides.py` regenerates the guides in `code/`, `results/`, `manuscript_figures/` and their figure subfolders. Edit this script to change those guides; the root README is edited directly.
- `python3 tools/record_current.py --manuscript-dir /path/to/paper` records source and PDF hashes, the three external figure hashes, current result numbering, data counts and the release manifest. `--data-only` preserves the manuscript record without checking it.
- `python3 tools/verify_manuscript_tables.py --manuscript-dir /path/to/paper` checks manuscript structure, displayed numerical table values, selected claims, and external figures against saved release data.
- `python3 tools/build_manuscript.py --source /path/to/combined.tex --output /path/to/empty-folder` copies the source and its active external figure files to an empty folder. `--compile` optionally runs an installed LaTeX engine three times.
- `external_manuscript.py` expands active source inputs for comparison.
- `python3 code/radius_enlargement_extra.py --verify` replays the two independent million-cloud circle checks; the script describes its conditional exponential order-statistic method.

After changing the release, update guides and hashes, then verify:

```sh
python3 tools/write_guides.py
python3 code/comparison_report.py --manuscript-dir /path/to/paper
python3 tools/record_current.py --manuscript-dir /path/to/paper
python3 code/verify_release.py --paper-only --manuscript-dir /path/to/paper --report /tmp/edge-budget-paper-check.json
python3 code/verify_release.py --manuscript-dir /path/to/paper --report /tmp/edge-budget-complete-check.json
```

A data-only record does not confirm that a changed paper matches the release. Reports for `--paper-only` must be outside the release; `file_manifest.json` excludes itself, the complete verification report and diagnostic logs. Record again after any later release edit.
'''
}

def main():
    for folder, text in GUIDES.items():
        target=ROOT/folder
        target.mkdir(parents=True, exist_ok=True)
        (target/'README.md').write_text(text)

if __name__=='__main__':
    main()
