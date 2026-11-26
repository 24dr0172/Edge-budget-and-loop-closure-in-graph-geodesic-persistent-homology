# Loop detection in graph-geodesic persistent homology

Code and numerical evidence for the manuscript by Shivam Singh, Annesha Sen and S. P. Tiwari. The accompanying paper studies loop detection in graph-geodesic persistent homology for samples of simple closed curves, including circle orderings, radius enlargement, radial noise, and numerical locality checks.

This release accompanies the 32-page manuscript. The manuscript record identifies the LaTeX source, compiled PDF and figure by SHA-256. The manuscript and its Overleaf project are distributed in the companion manuscript package. [Manuscript record](manuscript_record.json) and [v1.0.2 repository snapshot](https://github.com/24dr0172/Edge-budget-and-loop-closure-in-graph-geodesic-persistent-homology/tree/v1.0.2/final_release).

## Install and verify

The commands below create a Python 3.12 virtual environment alongside the release. The manuscript check also uses the Poppler commands `pdfinfo` and `pdftotext`. [Pinned Python dependencies](code/requirements.txt).

```sh
python3.12 -m venv ../.geodesic-env
. ../.geodesic-env/bin/activate
python3 -m pip install -r code/requirements.txt
python3 tools/verify_release.py
```

This command checks release hashes, rank summaries and selected seed replays, noise and winding records, peanut diagrams, and numerical-check provenance. Run `python3 tools/verify_paper.py` from this folder in the complete companion package to compare the manuscript's reported results with the saved evidence and verify the PDF and Overleaf files. The default manuscript directory is `../main-paper/`. Use `--manuscript-dir /path/to/main-paper` to verify another copy of the recorded manuscript alongside the same Overleaf project.

## Data map

- **Table 1.** Table 1 uses the uniform rows of the paired rank sweep. Its arrays retain the connecting rank, closure rank and largest gap for each sample, including the cosine and trough profiles. The first three sizes use 20,000 point sets per profile and the last uses 5,000. [Paired rank sweep](results/detection_rank_sweep/).
- **Appendix E.** Full diagrams for 60 clean-circle point sets at the first connected prefix, just before closure, and at closure in each ordering. The records also include the smallest connecting integer k for the two kNN rules. [Clean-circle events and counts](results/clean_closure/).
- **Table 2.** Paired radial-noise results for 400 clean point sets, four amplitudes and both kNN rules. The release also contains graph winding checks for the 258 graphs at amplitude 0.01 that met the half-lifetime cutoff. [Paired results](results/radial_noise_knn.csv) and [winding checks](results/radial_noise_winding.csv).
- **Table 3 and Appendix D.** Circle misses are calculated from the consecutive chords. The entries at α = 2 for n = 4,096 and 32,768 each come from one million circle point sets. The ellipse and peanut entries use the consecutive-chord proxy with one million samples per curve and size. Full-diagram spot checks and the exceptional-sample diagrams assess agreement with diagram emptiness. [Circle runs](results/radius_enlargement/), [larger circle runs](results/radius_enlargement_extra/), [smooth-curve batches and diagrams](results/smooth_curve_enlargement/) and [exceptional samples and interpolation checks](results/smooth_curve_checks/).
- **Appendix A.** The locality runs include the grid estimates used in Appendix D. [Locality results](results/noncircular_locality/).
- **Limitations.** At a = 0.9 and n = 512, the two rules use the same twenty samples. The first non-empty prefix came before cycle closure seven times for symmetric kNN and twice for mutual kNN. [First non-empty prefixes](results/first_nonempty/) and [paired closure diagrams](results/peanut_neck/).
- **Finite examples.** The eight-point certificate for Remark 2 is under `dense_shortcut`. The same report contains checks for Appendix B and Example 1. [Small examples](results/small_examples_verification.json), [certificate for Remark 1](results/polygon_certificate.json) and [circle reductions](results/intrinsic_circle_verification.json).
- **Figure 2.** Figure 1 is drawn directly in the manuscript source. [Arc-support figure](manuscript_figures/arc_support_fan.png).

Load `.npz` files with `numpy.load(path, allow_pickle=False)`. The configuration and individual studies record seeds, tolerances and sampling rules. The manifest records SHA-256 hashes for the release files. [Configuration](code/config.json) and [manifest](manifest.json).

## Sampling and numerical checks

The two nearest-neighbour rules use the same point sets. For the additional profiles, a uniform phase shifts the cosine density 1 + 0.8 cos(2π(s − φ)) and centres the quarter-circle trough at φ + 1/2. The trough density is 1/7 inside the trough and 9/7 elsewhere. Profiles within a seed share the same uniform draws. The separate full-diagram experiment used 20 new samples per profile at n = 512.

In the radial-noise study, each displacement is the amplitude times an independent uniform draw on [−1, 1] in the outward unit normal direction. Both rules use the same noisy point set at each amplitude, and each clean sample reuses its noise draws across amplitudes. The lifetime cutoff uses the noise-free arc-length interval of that same sample.

The ellipse and peanut use separate random streams for radius enlargement, with a total of eight million point sets across both curves and all sizes. The three enlargement factors reuse each point set. At n = 128, the grid cutoff flagged two ellipse samples and 39 peanut samples. Diagram emptiness agreed with the chord proxy in all 123 exceptional cases and in 13,050 separate checks at n = 128 and 512. For each curve, the separate checks used the first 2,000 point sets at n = 128 and the first 150 plus 25 additional α = 2 proxy misses at n = 512. The additional misses were selected from the first 20,000 point sets, excluding the first 150. Every selected point set was checked at all three factors.

The arc tables use 400,001 points and linear inverse interpolation, with trapezoidal integration for the ellipse and cumulative polyline lengths for the peanut. For the n = 32,768 interpolation check, Simpson integration and cubic inverse interpolation on grids of 1,600,001 and 3,200,001 points gave coordinate differences below 5.4 × 10⁻¹² on 200,001 probe points per curve. The eight largest-gap chords on 1,000 new samples per curve agreed within 4.8 × 10⁻¹². All 6,000 classifications agreed. Near a peanut inflection, interpolation error can exceed the chord–arc difference.

## Reproduce

The code guide gives full-cohort commands and single-seed replays for the clean-circle, peanut and radial-noise studies. Full runs can take hours. Run generators in a copy to preserve the supplied outputs under `results/`. [Code guide](code/README.md).

## AI use, citation and license

Generative-AI tools (Anthropic, OpenAI) were used for language editing, code prototyping and checking for gaps and errors in proofs. The authors checked the definitions, theorem statements, proofs and experimental results and take full responsibility for the manuscript.

Cite Singh, Sen and Tiwari, *Loop detection in graph-geodesic persistent homology*; see [CITATION.cff](CITATION.cff). Code is under MIT; data, figures and documentation are under CC BY 4.0. See [LICENSE](LICENSE).
