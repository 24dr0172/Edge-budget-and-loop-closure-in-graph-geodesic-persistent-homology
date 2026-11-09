# Edge budget and loop closure in graph-geodesic persistent homology

Code and saved data for the manuscript by **Shivam Singh, Annesha Sen and S. P. Tiwari**. The manuscript and its Overleaf project are supplied separately; this folder contains code, data, checks and figure assets.

The manuscript studies noise-free samples of a simple closed curve. It separates first connectivity from the addition of the last consecutive-sample edge, proves circle and smooth-curve results for loop retention, and compares persistence-diagram error with the intrinsic metric on the same sample. Its scope and limitations are stated in the manuscript.

The circle is an exact model for connectivity and closure; its Euclidean and intrinsic Rips filtrations differ by a monotone change of threshold. The peanut comparison targets intrinsic diagrams and reports positive Euclidean bar counts as well as errors.

## Verify the release

Use Python 3.12 and the exact versions in [code/requirements.txt](code/requirements.txt). From this folder:

```sh
python3 -m venv /path/to/edge-budget-env
. /path/to/edge-budget-env/bin/activate
python3 -m pip install -r code/requirements.txt
python3 code/verify_release.py --paper-only --manuscript-dir /path/to/paper --report /tmp/edge-budget-paper-check.json
```

The manuscript directory must contain `combined.tex`, `combined.pdf` and its `figures/` folder. The verifier compares their recorded hashes, theorem and figure structure, and displayed table values with the saved data. It checks release-file hashes and selected scientific evidence; it does not regenerate the full experiments or formally verify the proofs. The paper-only check also runs further evidence checks; its report identifies what was checked.

For all saved-data checks, including the additional studies, use `python3 code/verify_release.py --manuscript-dir /path/to/paper --report /tmp/edge-budget-complete-check.json`. To check the data without a manuscript, use `--data-only` instead of `--manuscript-dir`. The [Dockerfile](Dockerfile) provides a pinned environment. See the [code guide](code/README.md) for reproduction and the [tools guide](tools/README.md) for refreshing records after changes.

## Manuscript results and evidence

- **Intrinsic bar:** the winding-triple formula gives the exact death and positivity test for the finite intrinsic sample. [Exact boundary reductions](results/intrinsic_circle_verification.json), replayed by [verify_intrinsic_circle.py](code/verify_intrinsic_circle.py), check every distance threshold on 190 rational samples over F2, F3 and F5.
- **Perturbation remark:** strict chord and route-versus-edge comparisons with minimum margin eta remain valid for point motions smaller than eta/(2n), with order, minor arcs and the geometric hypotheses retained. Each preserved connected-prefix metric and its Rips diagram move by at most 2(n-1) times the point-motion bound.
- **Theorems 1–3 and Corollary 2:** under the stated insertion condition, cutting a connected prefix before ring closure at its missing ring edge makes the shortest-path metric Robinsonian, so positive-degree persistence is empty. For every chord-weighted graph that contains the ring, the ring deaths satisfy B <= tau_E <= tau_G, and adding edges cannot raise the graph ring death. Arc support makes the H1 persistence module the interval [B, tau_G), possibly empty. On the circle, the first nonempty connected exact prefix identifies closure under the stated positive-bar and distinct-chord assumptions. The [rank sweep](results/extra_k/rank_sweep/samples.csv) and [retention counts](results/retention_table/retention_counts.csv) support the circle measurements in Table 3.
- **Theorem 4:** for independent uniform samples on a fixed smooth closed curve, the first connected radius graph misses the loop with high probability. On that event, a connected radius graph has a nonempty loop diagram if and only if it contains the ring, so the first such graph is the closure graph. The event statements fail with probability at most n(1-u_0)^(n-1), where u_0 is the minimum of the locality radius and one third of the global chord–arc ratio. The theorem also gives finite-sample two-sided bounds for radius enlargement; at a fixed factor alpha > 1, the miss probability is asymptotic to Gamma(alpha+1) n^(1-alpha). [Circle results](results/radius_enlargement/summary.csv) and [ellipse and peanut proxy results](results/smooth_curve_enlargement/summary.csv) support Table 5; the latter are chord-event proxies, with full-diagram spot checks.
- **Lemma 9 and Corollaries 3–4:** for the circle graph orders, the optimum over connected exact prefixes is the smaller of the closure error and the empty-diagram cost; every post-closure error is below the empty-diagram cost exactly when g < psi(b_star), with psi(b) = 2 sin(pi b)/pi - b. Under that condition error is nondecreasing after closure and closure is a minimizer; a sufficient gap threshold is g_0 approximately 0.1513. The error bound after closure is quadratic in the longest shorter curve arc spanned by a retained edge; the ceiling is the Euclidean direct-matching cost. [Saved closure errors](results/exact_optimum/records.csv), [rate records](results/post_closure_rate/records.csv) and [exact certificates](results/rational_certificates/certificates.json) support the numerical statements. The exact equally spaced example and its perturbation argument are in the paper.
- **Corollary 5 and the peanut study:** the paper's known-curve certificate uses the edgewise distortion bound U_G. [Closure measurements](results/peanut_neck/closure_error_bound.csv) provide the four aggregated noncircular rows of Table 4; [first-nonempty records](results/first_nonempty/trials.csv) provide the stopping and radius-doubling counts. The seeds are paired across widths and methods; the accuracy target is the intrinsic diagram. Euclidean positive-bar counts are 1 at waist parameter 0, range from 2 to 3 at 0.4 and equal 2 at the narrower waists. At waist parameter a >= 0.6, Euclidean error equals the empty-diagram cost in every recorded case. The saved files also contain winding-based bounds as additional data.
- **Retention heuristic:** [the practical-claims record](results/practical_claims_verification.json) contains the independent-key model for ties among consecutive-pair ranks. Fixed-neighbor and bridge-free measurements are separate study data.
- **Examples and locality:** [the exact polygon certificate](results/rational_certificates/certificates.json), [seven-point witness](results/useful_diagram/rational_counterexample.json), [ellipse certificate](results/natural_instance_verification.json), and [two-curve locality records](results/noncircular_locality/locality_check_results.csv) accompany Remark 1 and Appendices A–E. Locality thresholds are grid-search estimates with interpolation and are not certified bounds. Proposition 1 also gives an explicit sufficient radius when a lower chord–arc ratio and an upper curvature bound are certified.

## Tables and figures

| Manuscript item | Source or evidence |
|---|---|
| Table 1: established ingredients | Written in the manuscript |
| Table 2: notation | Written in the manuscript |
| Table 3: circle retention | [retention_compact.csv](results/retention_table/retention_compact.csv), [raw ranks](results/extra_k/rank_sweep/samples.csv) |
| Table 4: peanut closure, one row per noncircular waist | [closure_error_bound.csv](results/peanut_neck/closure_error_bound.csv) and [trials](results/peanut_neck/trials.csv) |
| Table 5: finite-size enlargement | [circle rates](results/radius_enlargement/summary.csv), [additional million-cloud circle checks](results/radius_enlargement_extra/summary.csv), [smooth-curve rates](results/smooth_curve_enlargement/summary.csv) |
| Figure 1: connectivity and closure | Drawn inline in the manuscript source |
| Figure 2: supporting fan | [arc_support_fan.png](manuscript_figures/arc_support_fan.png) |
| Figure 3: error rate | [post_closure_rate.png](manuscript_figures/post_closure_rate.png), [records](results/post_closure_rate/records.csv) |
| Figure 4: peanut comparison | [peanut_comparison.png](manuscript_figures/peanut_comparison.png), [saved diagrams and bounds](results/peanut_neck/closure_error_bound.json) |

The three PNGs used by the manuscript match its external figure files byte for byte. Figure 1 is inline. [PDF exports](results/figure_exports/README.md) are optional. Additional plots, including a winding-bound comparison and extra neighbor layers, are in [additional_figures](results/additional_figures/README.md).

## Reproduction and additional data

`python3 code/reproduce.py --check` checks the environment and inputs. `python3 code/reproduce.py --output /path/to/new-run --workers 6` generates the studies in a new folder and takes hours. It preserves the supplied Figure 4 PNG and PDF while rebuilding its data, so the manuscript comparison can use the frozen asset. The [code guide](code/README.md) distinguishes rebuilding plots from regenerating samples.

The release also includes [noise and graph comparisons](results/comparisons/README.md), [other plots](results/additional_figures/README.md), and winding, bridge-free and ellipse-window studies. These studies are separate from the theorems of the manuscript. The [results guide](results/README.md) describes the provenance and file records.

## Use of generative-AI tools

Generative-AI tools (Anthropic, OpenAI) were used for language editing, code prototyping and to review the proofs. All definitions, theorem statements, final proofs and experimental verifications were written, checked and finalized by the authors, who take full responsibility for the manuscript.

## Citation and license

Singh, S., Sen, A., and Tiwari, S. P. *Edge budget and loop closure in graph-geodesic persistent homology*. Code and data release, 2026. See [CITATION.cff](CITATION.cff). Code is under MIT; data, figures and documentation are under CC BY 4.0. See [LICENSE](LICENSE).
