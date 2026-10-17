# Manuscript figure assets

All manuscript assets are in this folder. Every PDF asset has a matching 300-dpi PNG with the same base name. All three panels of Figure 1 are drawn in the LaTeX source, using the same serif fonts and circle guides. `unsupported-shortcut.pdf` and `unsupported-shortcut.png` provide a standalone companion illustration of the third panel.

- Figure 2: `Fig1.pdf` (and PNG preview), clean-circle grid by profile: connectivity, maximum-distance error and diagram error.
- Figure 3: `noise-eph.pdf` (and PNG counterpart), median diagram error at delta=0.01 under radial and isotropic noise.
- Release error curves: `noise-all-errors.pdf` (and PNG counterpart), all three errors at representative noise amplitudes.
- Figure D.1: `appendix_extra_k.pdf` (and PNG preview), paired extra-k empirical CDFs.
- Additional diagnostic asset: `Fig2.pdf` (and PNG preview), radius/OR first-connectivity radial noise counts and lifetime ratios, from results/matched_noise/summary.csv.

`figure_index.json` maps the displayed manuscript figures and release assets to their sources. Run `python3 code/build_figures.py` from final_release to rebuild all assets. PDF timestamps are omitted for deterministic exports in the recorded environment.
