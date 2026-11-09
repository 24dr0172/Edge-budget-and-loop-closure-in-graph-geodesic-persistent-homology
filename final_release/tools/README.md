# Release tools

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
