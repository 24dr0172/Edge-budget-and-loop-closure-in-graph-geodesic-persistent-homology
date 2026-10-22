# Manuscript figures

This folder contains exactly the three PNGs read by the portable manuscript:

- **Figure 1:** `closure_schematic.png`, the supplied 600-dpi closure schematic.
- **Figure 2:** `arc_support_fan.png`, the supplied 600-dpi supporting-fan drawing; vertex labels match Definition 4.
- **Figure 3:** `retention_extra_layers.png`, generated at 300 dpi from [saved rank records](../results/extra_k/rank_sweep/samples.csv) by `code/retention_extra_layers.py`. The 144 group/layer records and Wilson intervals are in [retention_extra_layers](../results/retention_extra_layers/retention.csv).

The manuscript copies these into `figures/`. Captions, labels and cross-references remain in `combined.tex`. It needs no native drawings, PDF graphics, external tables or appendix files. The figures match the release byte for byte.

Native drawing sources belong to the separate Overleaf project and are excluded here. [Archival PDF exports](../results/figure_exports/README.md) are included as image masters. [Additional study plots](../results/additional_figures/README.md) are kept separately and are not figures of the focused paper.
