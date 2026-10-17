# Figure assets

Each PDF in this folder has a matching 300-dpi PNG with the same base name.

- Figure 1 is drawn directly in `../combined.tex`, relative to `final_release/`. `unsupported-shortcut.pdf` and its PNG provide a separate illustration of its third panel.
- Figure 2: `Fig1.pdf` and `Fig1.png` show the clean-circle budget grid.
- Figure 3: `noise-eph.pdf` and `noise-eph.png` show median diagram error at delta=0.01 for both noise models.
- Figure D.1: `appendix_extra_k.pdf` and its PNG show the paired extra-neighbor distributions.
- `noise-all-errors.pdf` and its PNG show all three error measures at representative noise amplitudes.
- `Fig2.pdf` and `Fig2.png` show the matched radial-noise first-connectivity checks.

`figure_index.json` maps each asset to its data source. Run `python3 code/build_figures.py` from `final_release/` to rebuild the assets. PDF timestamps are omitted so PDF exports are repeatable in the recorded environment.
