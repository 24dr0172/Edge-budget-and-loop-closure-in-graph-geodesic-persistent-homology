# Edge budget and loop closure in graph-geodesic persistent homology

Code, saved data and verification records for the paper by **Shivam Singh, Annesha Sen and S. P. Tiwari**.

The focused paper contains all proofs and the numerical certificates needed to follow its examples. It studies finite, noise-free samples of one closed curve. This release contains its evidence and reproducible scripts, together with additional studies on selection, noise and noncircular geometry. Those additional studies are not claims of the focused paper. No separate supplement is required or distributed. The manuscript and Overleaf project are supplied separately for optional comparison. This release contains no paper, LaTeX document source, bibliography source or manuscript ZIP.

## Terms and scope

- The **ring** joins consecutive samples in the given curve order. **Ring closure** is the first budget containing all its edges.
- An **edge budget M** selects the first M sample pairs in an edge order; shortest paths in that graph give the distances used by persistence.
- **Radius**, **OR-kNN** and **mutual-kNN** orders use chord length, the smaller endpoint neighbor rank and the larger endpoint rank, respectively (Definition 2).
- A **layer** inserts all pairs with the same rank key or chord length. Exact insertions and whole layers are distinct.
- **k_conn**, **k_close** and **k_br** are the first connected layer, the first ring-containing layer and the first connected layer with no bridge.
- **Retention** means the ring is already present. In the circle rank sweep every generated gap is below one third, so it is also positive-loop retention.
- **Short-hops-first** concerns the insertion order before closure (Definition 3); **arc support** concerns triangle fans replacing edges after closure (Definition 4).
- **e_PH** is bottleneck diagram error normalized by the intrinsic sample diameter; **P_star** is the cost of reporting no loop.
- The radius-enlargement proposition concerns independent uniform circle samples and a proportional increase in radius. It does not give a fixed-extra-k guarantee.

The main notation is in Table 2. Proofs are in the paper; numerical checks corroborate them.

## Install

Use Python 3.12 in a clean environment:

```sh
python3 -m venv /path/to/edge-budget-env
. /path/to/edge-budget-env/bin/activate
python3 -m pip install -r code/requirements.txt
```

Checking saved data requires no compiler. Regenerating the bridge-free study requires a C++17 compiler. The [Dockerfile](Dockerfile) supplies the pinned environment and compiler:

```sh
docker build -t edge-budget .
docker run --rm edge-budget
```

## Verify

Run commands from this release folder:

```sh
# Evidence for the focused paper and every payload hash
python3 code/verify_release.py --paper-only --data-only --report /tmp/edge-budget-paper-check.json
# Every saved-data check, including the additional studies
python3 code/verify_release.py --data-only --report /tmp/edge-budget-complete-check.json
# Also compare the external manuscript with this release
python3 code/verify_release.py --manuscript-dir /path/to/paper --report /tmp/edge-budget-manuscript-check.json
# Validate the generation pipeline without generating measurements
python3 code/reproduce.py --check
```

These checks do not regenerate sample cohorts. Complete verification independently replays selected graph diagrams and all 240 immediately preceding cross-neck prefixes of the additional varying-waist study. The focused mode omits the additional noise and neck cohorts; shared certificate checks also verify retained additional exports. Write focused reports outside the release. A complete run without `--report` updates `results/verification_report.json`.

The portable manuscript consists of `combined.tex` and three PNGs in `figures/`: `closure_schematic.png`, `arc_support_fan.png` and `retention_extra_layers.png`. All tables, proofs, appendices and bibliography are inside `combined.tex`; it has no local TeX inputs, external bibliography or PDF graphics. A standard LaTeX installation is sufficient; the paper does not require this release or TikZ. For optional verification, also supply the compiled `combined.pdf`. The three PNGs match [manuscript_figures](manuscript_figures/README.md) byte for byte, and the inlined Table 3 is compared with its numerical CSV export. Paths in [manuscript_record.json](results/manuscript_record.json) are relative to the manuscript folder.

## Paper tables and figures

| Paper item | Source and saved evidence | Rebuild or verify |
|---|---|---|
| Table 1 (established ingredients and verifications) | written directly in the manuscript | Analytical exposition |
| Table 2 (notation, placed before the main proofs) | written directly in the manuscript | Definitions and cross-references |
| Table 3 (retention at first connectivity, four sample sizes) | [retention_compact.csv](results/retention_table/retention_compact.csv); [counts and Wilson intervals](results/retention_table/retention_counts.csv), from [raw ranks](results/extra_k/rank_sweep/samples.csv) and [bridge records](results/comparisons/bridgeless/sweep_samples.csv) | `python3 code/build_retention_table.py` |
| Figure 1 (closure and shortcut schematic) | [closure_schematic.png](manuscript_figures/closure_schematic.png) | Supplied image; editable source is kept separately with the Overleaf project |
| Figure 2 (supporting fan) | [arc_support_fan.png](manuscript_figures/arc_support_fan.png) | Supplied image; editable source is kept separately with the Overleaf project |
| Figure 3 (extra neighbor layers) | [retention_extra_layers.png](manuscript_figures/retention_extra_layers.png), from [144 group/layer records](results/retention_extra_layers/retention.csv); [settings](results/retention_extra_layers/config.json) | `python3 code/retention_extra_layers.py` |

Other evidence used by the paper:

- **Examples 1 and 2:** [exact separation certificate](results/separation_verification.json), [rational ellipse certificate](results/natural_instance_verification.json), and [small-example checks](results/small_examples_verification.json). The proofs and bounds remain in the paper.
- **Proposition 1:** analytical radius-enlargement proof. `verify_paper_evidence.py` checks the beta integral numerically against its gamma expression.
- **Section 7.1:** [configuration](code/config.json), [pinned dependencies](code/requirements.txt), and [cutoff validation](results/cutoff_validation.json). The clean normalized-error tolerance is 5e-8; arc-support tolerance is 1e-10 raw units. Positive returned bars are not removed by a lifetime cutoff.
- **Section 7.2:** [24 primary retention cells](results/rank_table.csv), [raw ranks and seeds](results/extra_k/rank_sweep/config.json), and [full-diagram event checks](results/clean_closure/event_table.csv) on 60 separate n=512 clouds.
- **Section 7.3:** [fixed-k counts](results/retention_table/retention_fixed_k.csv), [bridge-free summaries](results/comparisons/bridgeless/sweep_summary.csv), and the extra-layer records behind Figure 3.
- **Section 7.4 and Corollaries 1-2:** [checked error extremes](results/practical_claims_verification.json), [clean-circle cells](results/clean_closure/cells.csv), and [saved post-closure prefixes](results/useful_diagram/saved_prefixes.csv). The 0.044 observation concerns the stored budget grid. The 0.115 dense-sample limit is not a finite-sample universal bound.
- **Section 8.2:** [independent-key predictions](results/extra_k/rank_sweep/iid_predictions.csv). The displayed equation is a heuristic: actual keys share spacings.
- **Appendix A:** [two-curve locality records](results/noncircular_locality/locality_check_results.csv) and [threshold summaries](results/noncircular_locality/locality_check_summary.json). The continuum thresholds are grid estimates, not certified global minima.
- **Appendix B:** [small examples](results/small_examples_verification.json), including polygon, rectangle and rank-layer ties.
- **Appendix C:** [seven-point witness](results/useful_diagram/rational_counterexample.json), [recomputed values](results/seven_point_check/seven_point_numbers.json), and [independent calculations](results/independent_checks/independent_checks_output.txt).

## Rebuild saved-data outputs

```sh
python3 code/reproduce_results.py --paper-table-only
python3 code/verify_separation.py
python3 code/verify_natural_instance.py
python3 code/verify_small_examples.py
python3 code/check_practical_claims.py
python3 code/build_retention_table.py
python3 code/retention_extra_layers.py
python3 code/check_seven_point_example.py
python3 code/verify_retention_cost_witness.py --write
python3 tools/write_guides.py
python3 tools/record_current.py --manuscript-dir /path/to/paper
python3 code/verify_release.py --manuscript-dir /path/to/paper
```

Without a manuscript, use `--data-only` for the final two commands. The recorder retains the saved manuscript identities in that mode. The compact table and three-layer display are numerical CSV exports; no document source is generated in the release. [The code guide](code/README.md#record-and-verify-after-a-change) gives the complete additional-study rebuild sequence.

The Overleaf project, paper source and editable drawing source are kept separate from this release. Image masters remain in `results/figure_exports/`. `tools/build_manuscript.py` can copy an external manuscript and its three PNGs into an empty folder for optional local compilation; the release itself contains no LaTeX documents or paper PDF.

## Additional studies

The release preserves the existing measurements, figures and complete checks for noise, diagram-based selectors, CkNN and shuffled scales, the ellipse sweep, the peanut/trough cohort, the varying-waist neck study, broader size sweeps and comparisons. Their selectors and figures are additional-study labels, not numbered results of the focused paper. Their plots are kept separately in [results/additional_figures](results/additional_figures/README.md), outside the three-file manuscript PNG inventory. No mathematical claim in the paper depends on an unprovided document.

- [Additional study comparisons](results/comparisons/README.md), including [noise selectors](results/noise/selectors.csv).
- [Varying-waist study settings](results/peanut_neck/config.json) and [all event records](results/peanut_neck/trials.csv).
- [Peanut/trough settings](results/peanut_trough/config.json) and [closure records](results/peanut_trough/closure.csv).
- [Code guide](code/README.md), [results inventory](results/README.md), [figure inventory](manuscript_figures/README.md), and [tools](tools/README.md).

## Full reproduction from seeds

```sh
python3 code/reproduce.py --check
python3 code/reproduce.py --output /path/to/new-run --workers 6
```

The second command regenerates every study and takes hours. Use `--resume` for an interrupted run. `python3 code/smoke_test.py` is a one-sample check that writes nothing.

## Use of generative-AI tools

Generative-AI tools (Anthropic, OpenAI) were used for language editing, code prototyping and code review, and to produce early drafts of some proofs. These drafts were substantially rewritten by the authors. All definitions, theorem statements, final proofs and experimental verifications were written, checked and finalized by the authors, who take full responsibility for the manuscript.

## Citation and license

Singh, S., Sen, A., and Tiwari, S. P. *Edge budget and loop closure in graph-geodesic persistent homology*. Code and data release, 2026. See [CITATION.cff](CITATION.cff).

The release is hosted on [GitHub](https://github.com/24dr0172/Edge-budget-and-loop-closure-in-graph-geodesic-persistent-homology). Code is under MIT; data, figures and documentation are under CC BY 4.0. See [LICENSE](LICENSE). Third-party dependencies retain their licenses.
