# Loop detection in graph-geodesic persistent homology

Code and numerical evidence for the manuscript by Shivam Singh, Annesha Sen and S. P. Tiwari. The accompanying paper studies loop detection in graph-geodesic persistent homology for samples of simple closed curves, including circle orderings, radius enlargement, radial noise, and numerical locality checks.

## Install and verify

Use Python 3.12 and the versions in [code/requirements.txt](code/requirements.txt). The manuscript check also uses the Poppler commands `pdfinfo` and `pdftotext`. From this folder:

```sh
python3 -m venv /path/to/geodesic-env
. /path/to/geodesic-env/bin/activate
python3 -m pip install -r code/requirements.txt
python3 tools/verify_release.py
```

This command checks release hashes, rank summaries and selected seed replays, noise and winding records, peanut diagrams, and numerical-check provenance. To check the manuscript and Overleaf project as well, run `python3 tools/verify_paper.py` from the complete companion package. That check reads the manuscript in `../main-paper/`, compares its displayed results with the saved evidence, and checks the PDF and Overleaf copies. `--manuscript-dir /path/to/main-paper` selects another manuscript directory. These checks do not prove the theorems.

## Data map

- **Table 1 and Theorem 4 comparisons:** [paired rank sweep](results/detection_rank_sweep/) and its [separate uniform-circle comparison](results/detection_rank_sweep/uniform_100k/). The rank arrays contain the connecting and closure ranks and largest gap for each sample. The first three sizes use 20,000 point sets per profile; the last uses 5,000. The uniform comparison uses 100,000 per size.
- **Clean-circle full diagrams:** [per-sample events and counts](results/clean_closure/) for 60 point sets and the three orderings.
- **Table 2, radial noise:** [paired results](results/radial_noise_knn.csv) for 400 point sets, four amplitudes and both kNN rules; [winding checks](results/radial_noise_winding.csv) for the 258 long-interval graphs at amplitude 0.01.
- **Table 3, radius enlargement:** [circle gaps and summaries](results/radius_enlargement/), [larger circle checks](results/radius_enlargement_extra/), and [smooth-curve batches and diagrams](results/smooth_curve_enlargement/). Ellipse and peanut table entries are consecutive-chord proxies; the spot checks compare these with full diagrams.
- **Appendix A, locality and peanut detection:** [locality results](results/noncircular_locality/), [peanut closure records and diagrams](results/peanut_neck/), and [first-nonempty-prefix records and diagrams](results/first_nonempty/).
- **Appendix E, numerical checks:** [exceptional point sets and interpolation checks](results/smooth_curve_checks/). Geometric constants are grid estimates, and the interpolation comparisons are numerical rather than interval certificates.
- **Finite examples:** [circle reductions](results/intrinsic_circle_verification.json), [polygon certificate](results/polygon_certificate.json), and [small-example checks](results/small_examples_verification.json).
- **Figure 2:** [arc_support_fan.png](manuscript_figures/arc_support_fan.png). Figure 1 is drawn directly in the manuscript source.

Load `.npz` files with `numpy.load(path, allow_pickle=False)`. Seeds, tolerances and sampling rules are documented in [code/config.json](code/config.json), the individual scripts, and the study-specific protocol files. `manifest.json` records SHA-256 hashes of the release files.

## Reproduce

The [code guide](code/README.md) lists commands for every study. Individual generators that write to `results/` replace their corresponding outputs; copy the release first if the supplied records must remain unchanged. Full runs can take hours. Selected seeds can be replayed separately for the clean-circle, peanut and radial-noise studies.

## AI use, citation and license

Generative-AI tools (Anthropic, OpenAI) were used for language editing, code prototyping and checking for gaps and errors in proofs. OpenAI Codex also assisted manuscript revisions, including theorem wording, a proof calculation and experimental code. The authors checked the definitions, theorem statements, proofs and experimental results and take full responsibility for the manuscript.

Cite Singh, Sen and Tiwari, *Loop detection in graph-geodesic persistent homology*; see [CITATION.cff](CITATION.cff). Code is under MIT; data, figures and documentation are under CC BY 4.0. See [LICENSE](LICENSE).
