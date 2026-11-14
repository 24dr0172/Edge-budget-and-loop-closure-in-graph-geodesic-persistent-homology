# Edge count and loop detection in graph-geodesic persistent homology

Code and saved data for the manuscript by **Shivam Singh, Annesha Sen and S. P. Tiwari**. The manuscript and its Overleaf project are supplied separately; this folder contains code, data, checks and figure assets.

The manuscript studies noise-free samples of a simple closed curve. It separates first connectivity from the addition of the last consecutive-sample edge, proves circle and smooth-curve results for loop retention, and compares persistence-diagram error with the intrinsic metric on the same sample. The radial-noise experiment measures both empty diagrams and intervals above two lifetime cutoffs. The repository URL uses the manuscript's earlier title.

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

- **Arc-length metric:** Lemma 3 gives the exact death and positivity test for the loop interval of the arc-length metric on a finite circle sample. [Exact boundary reductions](results/intrinsic_circle_verification.json), replayed by [verify_intrinsic_circle.py](code/verify_intrinsic_circle.py), check every distance threshold on 190 rational samples over F2, F3 and F5.
- **Theorems 1–3 and Corollary 2:** under the stated insertion condition, cutting a connected graph before cycle closure at its missing consecutive edge makes the shortest-path metric Robinsonian, so positive-degree persistence is empty. For a graph with Euclidean edge weights containing every consecutive edge, the deaths of the cycle class in the Euclidean and graph metrics satisfy B <= tau_E <= tau_G. Adding edges cannot raise tau_G. Arc support gives the H1 interval [B, tau_G), possibly empty. Under the circle hypotheses, the first connected graph with a nonempty loop diagram is the first graph containing every consecutive edge. The [rank sweep](results/extra_k/rank_sweep/samples.csv) and [retention counts](results/retention_table/retention_counts.csv) support the circle measurements in Table 1.
- **Theorem 4:** for uniform circle point sets, a Poisson approximation gives the probability of a tie at the largest consecutive rank up to a term tending to zero. For either nearest-neighbor rule, the probability of an empty diagram at the smallest connecting k stays bounded away from zero. The numerical check reports approximations for both rules at the four sample sizes in Table 1. [Numerical tie check](code/check_tied_ranks.py)
- **Theorem 5:** for independent uniform samples on a fixed smooth closed curve, the first connected ε-neighborhood graph misses the loop with high probability. On the same event, a connected ε-neighborhood graph detects the loop if and only if it contains every consecutive edge. The event statements fail with probability at most n(1-u_0)^(n-1), where u_0 is the minimum of the locality radius and one third of the global chord–arc ratio. The theorem also gives finite-sample two-sided bounds for radius enlargement; at a fixed factor alpha > 1, the miss probability is asymptotic to Gamma(alpha+1) n^(1-alpha). [Circle results](results/radius_enlargement/summary.csv) and [ellipse and peanut proxy results](results/smooth_curve_enlargement/summary.csv) support Table 4; the latter are chord-event proxies, with full-diagram spot checks.
- **Lemma 9 and Corollaries 3–4:** for the circle edge orderings, the optimum over connected graphs on the first M edges is the smaller of the closure error and the empty-diagram cost; every post-closure error is below the empty-diagram cost exactly when g < psi(b_star), with psi(b) = 2 sin(pi b)/pi - b. Under that condition error is nondecreasing after closure and closure is a minimizer; a sufficient gap threshold is g_0 approximately 0.1513. The error bound after closure is quadratic in the longest shorter curve arc spanned by a retained edge; the ceiling is the Euclidean direct-matching cost. [Saved closure errors](results/exact_optimum/records.csv), [rate records](results/post_closure_rate/records.csv) and [exact certificates](results/rational_certificates/certificates.json) support the numerical statements.
- **Corollary 5 and the peanut study:** the paper's error bound for a known curve uses the edgewise distortion bound U_G. [Closure measurements](results/peanut_neck/closure_error_bound.csv) provide the four aggregated noncircular rows of Table 3; [first-nonempty records](results/first_nonempty/trials.csv) provide the stopping and radius-doubling counts. The seeds are paired across widths and methods; the accuracy target is the intrinsic diagram. Euclidean positive-bar counts are 1 at waist parameter 0, range from 2 to 3 at 0.4 and equal 2 at the narrower waists. At waist parameter a >= 0.6, Euclidean error equals the empty-diagram cost in every recorded case. The saved files also contain winding-based bounds as additional data.
- **Radial noise:** the 400 paired point sets measure empty diagrams, quarter- and half-lifetime counts, and coverage of the underlying consecutive edges at the smallest connecting k. [Generator](code/radial_noise_knn.py), [records](results/radial_noise_knn.csv)
- **Examples and locality:** [the exact polygon certificate](results/rational_certificates/certificates.json), [seven-point witness](results/useful_diagram/rational_counterexample.json), [ellipse certificate](results/natural_instance_verification.json), and [two-curve locality records](results/noncircular_locality/locality_check_results.csv) accompany Remark 1 and Appendices A–E. Locality thresholds are grid-search estimates with interpolation and are not certified bounds. Proposition 1 also gives an explicit sufficient radius when a lower chord–arc ratio and an upper curvature bound are certified. The winding-check configuration preserves the checksum of the script used to generate its saved measurements and records the checksum of the current verification script separately.

## Tables and figures

| Manuscript item | Source or evidence |
|---|---|
| Table 1: circle loop detection | [retention_compact.csv](results/retention_table/retention_compact.csv), [raw ranks](results/extra_k/rank_sweep/samples.csv) |
| Table 2: radial noise | [radial_noise_knn.csv](results/radial_noise_knn.csv) |
| Table 3: peanut closure, one row per noncircular waist | [closure_error_bound.csv](results/peanut_neck/closure_error_bound.csv) and [trials](results/peanut_neck/trials.csv) |
| Table 4: finite-size enlargement | [circle rates](results/radius_enlargement/summary.csv), [additional million-cloud circle checks](results/radius_enlargement_extra/summary.csv), [smooth-curve rates](results/smooth_curve_enlargement/summary.csv) |
| Figure 1: connectivity and closure | Drawn inline in the manuscript source |
| Figure 2: supporting fan | [arc_support_fan.png](manuscript_figures/arc_support_fan.png) |
| Figure 3: error rate | [post_closure_rate.png](manuscript_figures/post_closure_rate.png), [records](results/post_closure_rate/records.csv) |
| Figure 4: peanut comparison | [peanut_comparison.png](manuscript_figures/peanut_comparison.png), [saved diagrams and bounds](results/peanut_neck/closure_error_bound.json) |

The three PNGs used by the manuscript match its external figure files byte for byte. Figure 1 is inline. [PDF exports](results/figure_exports/README.md) are optional. Additional plots, including a winding-bound comparison and larger nearest-neighbour graphs, are in [additional_figures](results/additional_figures/README.md).

## Reproduction and additional data

`python3 code/reproduce.py --check` checks the environment and inputs. `python3 code/reproduce.py --output /path/to/new-run --workers 6` generates the studies in a new folder and takes hours. It preserves the supplied Figure 4 PNG and PDF while rebuilding its data, so the manuscript comparison can use the frozen asset. The [code guide](code/README.md) distinguishes rebuilding plots from regenerating samples.

The release also includes [noise and graph comparisons](results/comparisons/README.md), [other plots](results/additional_figures/README.md), and winding, bridge-free and ellipse-window studies. These studies are separate from the theorems of the manuscript. The [results guide](results/README.md) describes the provenance and file records.

## Use of generative-AI tools

Generative-AI tools (Anthropic, OpenAI) were used for language editing, code prototyping and to review the proofs. OpenAI Codex was used to improve the manuscript's readability and presentation. The authors used OpenAI's ChatGPT (6.1 Sol) to generate some of the LaTeX code for the schematic diagrams and then checked the code. All definitions, theorem statements, final proofs and experimental verifications were written, checked and finalized by the authors, who take full responsibility for the manuscript.

## Citation and license

Singh, S., Sen, A., and Tiwari, S. P. *Edge count and loop detection in graph-geodesic persistent homology*. Code and data release, 2026. See [CITATION.cff](CITATION.cff). Code is under MIT; data, figures and documentation are under CC BY 4.0. See [LICENSE](LICENSE).
