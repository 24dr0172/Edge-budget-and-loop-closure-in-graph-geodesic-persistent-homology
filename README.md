# Edge budget and loop closure in graph-geodesic persistent homology

The code, data, settings and figures for the paper are in [`final_release/`](final_release/). The release contains three folders:

- [`code/`](final_release/code/): experiment generators, checks and settings.
- [`results/`](final_release/results/): saved measurements, diagrams, tables and verification reports.
- [`manuscript_figures/`](final_release/manuscript_figures/): all six figures in PDF and PNG format.

The [code guide](final_release/code/README.md) gives installation instructions and commands to generate all scientific results from seeds in a new, empty output folder.

After installing the dependencies, check the included results from this repository folder with:

```sh
python3 final_release/code/verify_release.py --data-only
```

The manuscript and source notes are kept separately and are not needed to generate the numerical results. The code guide explains how to add an optional comparison with an external manuscript. The saved verification report identifies the manuscript used for the paper comparison.
