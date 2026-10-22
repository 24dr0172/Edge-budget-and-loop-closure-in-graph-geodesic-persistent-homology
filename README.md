# Edge budget and loop closure in graph-geodesic persistent homology

The verified code, saved data, settings and figures for the focused noise-free paper are in [`final_release/`](final_release/). Start with the [release guide](final_release/README.md) for the paper-to-evidence map, installation and reproduction commands.

- [`code/`](final_release/code/): experiment generators, verification scripts and settings.
- [`tools/`](final_release/tools/): manuscript comparison, portable-build and figure-export helpers.
- [`results/`](final_release/results/): saved measurements, diagrams, tables and verification records.
- [`manuscript_figures/`](final_release/manuscript_figures/): exactly the three PNG figures used in the paper: `closure_schematic.png`, `arc_support_fan.png` and `retention_extra_layers.png`.

Plots from additional studies are separate in [`results/additional_figures/`](final_release/results/additional_figures/). Vector figure masters are in [`results/figure_exports/`](final_release/results/figure_exports/). These are image exports, not the paper or LaTeX source.

The manuscript and Overleaf project are kept separately and are not published in this code/data distribution. The release contains no `.tex` or `.bib` files, paper PDF or manuscript ZIP. The manuscript itself compiles from one TeX file plus three PNGs. Its tables, proofs, appendices and bibliography are inline. No separate supplement is required or distributed. The release retains additional-study evidence, clearly distinguished from the focused paper's claims.

After installing the pinned Python dependencies, check the paper's saved evidence from this repository folder with:

```sh
python3 final_release/code/verify_release.py --paper-only --data-only --report /tmp/edge-budget-paper-check.json
```

For every saved-data check, including the additional studies:

```sh
python3 final_release/code/verify_release.py --data-only --report /tmp/edge-budget-complete-check.json
```

The [release guide](final_release/README.md) also explains how to compare an external manuscript. See the release [citation](final_release/CITATION.cff) and [license](final_release/LICENSE).
