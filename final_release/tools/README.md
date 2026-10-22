# Tools

Run these from the release folder. The manuscript is an optional, separate input.

- `write_guides.py` writes the guides in `code/`, `results/`, `manuscript_figures/` and `tools/`. Edit its text, then run it; do not edit those guides directly.
- `record_current.py --data-only` records the hash of every release file and the output counts, without reading a manuscript. With `--manuscript-dir /path/to/paper` it also records the manuscript source, PDF, inputs and graphics. Those paths are relative to the manuscript folder.
- `verify_manuscript_tables.py --manuscript-dir /path/to/paper` compares the manuscript with the release: table and figure contents, numbering, the inequalities of Examples 1 and 2, the computational conventions and the reported values. It also checks the title, authors, URL and table and figure map in the root README, and the numbered results, absence of external text inputs, exactly three matching PNGs and current paper-to-evidence map. The complete verifier runs it when a manuscript is supplied.
- `external_manuscript.py` reads the manuscript with its `\input` files expanded and comments removed, for those comparisons.
- `build_manuscript.py --source /path/to/combined.tex --output /path/to/empty-folder` copies an external manuscript and its PNGs into an empty folder for local compilation (`--compile` runs `pdflatex`). It is not part of the numerical pipeline.

After changing files, rewrite the guides, record once, then run the complete verification, in the order given in the [code guide](../code/README.md#record-and-verify-after-a-change).
