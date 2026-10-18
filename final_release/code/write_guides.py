"""Write the plain-language guides shipped with the three-folder release."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
CODE='''# Code and reproduction

This folder contains the code for **Edge budget and loop closure in graph-geodesic persistent homology**. It can generate the paper's scientific results from the stated seeds and formulas, starting with empty results and figure folders. It does not need another release or experiment folder.

The release has three folders:

- `code`: experiment generators, checks and settings.
- `results`: coordinates, graph events, persistence diagrams, tables and check reports.
- `manuscript_figures`: figures in PDF and PNG format.

## Install

Use Python 3.12 and a C++17 compiler (`g++`). From the release folder, run:

```sh
python3 -m venv /path/to/edge-budget-env
. /path/to/edge-budget-env/bin/activate
python3 -m pip install -r code/requirements.txt
```

The compiler is needed for the bridge comparison. Its C++ source is included. The Python dependencies must be installed once; the experiments then run locally without network access.

## Generate everything from seeds

Choose a new, empty output folder outside this release:

```sh
python3 code/reproduce.py --output /path/to/new-run --workers 6
```

This command copies only `code`, creates empty `results` and `manuscript_figures` folders, and runs all experiments, summaries, figures and checks. It does not copy any measured CSV file, coordinate archive or persistence diagram from this release. The commands and their logs are recorded under `results/run_logs` in the new output.

Full persistence calculations are expensive. Allow several hours and several GB of memory; the time depends on the computer. Lower `--workers` if memory is limited. To continue an interrupted run with the same code and settings:

```sh
python3 code/reproduce.py --output /path/to/new-run --workers 6 --resume
```

Completed sample checkpoints are reused only within that output folder. Changing the code or settings requires a new output folder. The run also measures fresh execution times. Clock times depend on hardware and load, so the recorded seconds in the paper are not exact reproduction targets.

## Check the included results

```sh
python3 code/verify_release.py
```

This checks file hashes, seed reconstruction, saved diagrams and scores, graph events, selector rules, summary statistics and the small exact examples. If the external manuscript is present beside the release, its displayed tables and identity are checked too. It does **not** regenerate every large persistence diagram. Use `reproduce.py` for that.

To compare a manuscript kept elsewhere, supply its directory:

```sh
python3 code/record_current.py --manuscript-dir /path/to/paper
python3 code/verify_release.py --manuscript-dir /path/to/paper
```

That directory must contain `combined.tex` and `combined.pdf`. Neither file is copied into the release. Use `--data-only` with the verifier to check only the numerical release. The report says whether manuscript comparison was run.

To rebuild tables and figures from the included raw measurements:

```sh
python3 code/reproduce_results.py
python3 code/build_figures.py
python3 code/write_guides.py
python3 code/record_current.py
python3 code/verify_release.py
```

## Experiments in the paper

- **Clean circle and order conditions:** `generate_clean.py` samples the three density profiles, builds the radius, OR-kNN and mutual-kNN graphs, and computes connectivity, closure, errors and full diagrams. It also checks ordinary CkNN and a control obtained by shuffling its local scales. Closed-prefix diagrams and the seven-point example support the diagram comparisons.
- **Ranks and additional neighbours:** `check_knn_miss_candidate.py` generates the size sweep from seeds, computes the local rank formulas, and produces closure rates, Wilson intervals and the probability predictions. `phase1_bridgeless.py` generates the connected bridge-free comparison using the included `bridgeless_circle.cpp` helper and an independent graph check.
- **Noise:** `run_studies.py` generates radial and uniform-disk perturbations, evaluates the prescribed budget grid and exact diagnostic events, and computes full diagrams and both selection rules. It also computes the Euclidean baseline. `generate_matched_noise.py` independently generates the radial connectivity and closure diagnostic from the same seeds and noise draws.
- **Ellipse:** `generate_ellipse.py` generates the coordinates by numerical arclength inversion and checks connectivity and closure. `run_studies.py` then evaluates the post-closure graphs, local and nonlocal edges, and NN-crust reconstruction.
- **Exact and numerical checks:** `generate_witness.py` computes the seven-point example; `verify_separation.py` checks Example 1 in all cyclic orders; `verify_arc_support.py` checks the arc condition on small graphs and the ellipse samples; `check_cknn_precision.py` resolves the CkNN failure using exact rational bounds; `check_rank_pairing.py` checks the paired large-sample ranks. `validate_cutoff.py` compares the adaptive calculation with unrestricted Ripser on random graphs and selected sample graphs.
- **Timing:** `benchmark_diagrams.py` runs the four representative full-diagram blocks three times each.

Stored method labels are `R2_global` (radius), `OR_kNN`, `mutual_kNN`, `R1_local` (CkNN) and `R3_permuted` (the shuffled-scale control).

`experiment_core.py`, `paper1_methods.py` and `finite_diagram.py` provide the shared calculations. `generate_measurements.py` coordinates the independent sample jobs. `reproduce_results.py`, `summarize_studies.py` and the summary helpers form tables from those measurements. The `verify_*.py` scripts and audits check the outputs.

## Inputs, conventions and scope

`config.json` records seeds, density profiles, noise amplitudes, budget settings, permutations, numerical tolerances and package versions. `result_schemas.json` contains CSV column names only. This release reproduces the paper’s fixed protocol; its checks expect those settings and sample counts. The source notes are kept outside the release and are not inputs to the numerical experiments. The exact Example 1 certificate is generated by `verify_separation.py`.

There are 60 clean samples (20 seeds and three chosen density profiles). The noise study uses separate seeds. Its two noise models share the clean observations, giving 220 distinct clouds and 240 model/amplitude/seed rows. The same displacement vector is reused across amplitudes within each seed and model. Profiles in the rank sweep also share their uniform draws, so these comparisons are paired.

Distances and shortest paths use float64. Ripser reports float32 endpoints. A positive interval has finite death strictly greater than birth; no positive-lifetime filter is applied. The longest-lifetime selector uses the earliest grid budget within a raw lifetime tolerance of 1e-7 of the maximum. Diagnostic event budgets are excluded from both selectors and the grid oracle.

The adaptive H1 cutoff still computes the full finite H1 diagram: above the largest generating edge, no new H1 class can be born. A cutoff above that edge with no surviving H1 therefore completes the diagram. If needed, the routine extends to a cone radius. The validation script also runs unrestricted Ripser independently.

Full diagrams are computed through n=512. The size sweep through one million points uses ranks only. Bridge calculations extend through n=32,768. Numerical agreement is assessed at the paper's stated tolerances and displayed precision; elapsed times and compressed archive hashes can change between runs.

## Manuscript and notes

The manuscript and source notes are kept outside this release. They are not needed to generate the numerical results. The optional `--manuscript-dir` argument adds a comparison with the external paper. A complete fresh run also accepts this argument; without an external paper, its final report covers the numerical release only.

`build_manuscript.py` is an optional helper that takes an explicit external source and creates a temporary compilation folder outside the release. It does not change or bundle the manuscript. `generate_notes.py` likewise requires explicit external source and output paths and is not part of the numerical reproduction command.
'''
RESULTS='''# Results

These files contain the scientific measurements, diagrams and summaries for the accompanying paper. The complete generation command is in [the code guide](../code/README.md). No other experiment folder is required.

- `clean_closure`: the three clean-circle profiles; coordinates, graph events, budget cells, errors, closure selections and sample checkpoints.
- `order_checks`: ordinary CkNN certificates, the shuffled-scale control, its diagrams and the exact precision check.
- `useful_diagram`: full diagrams of the tested closed graph prefixes, their transitions and the seven-point example.
- `extra_k`: additional integer neighbours needed for closure. `candidate_audit` contains the seed-generated size sweep and probability tables.
- `matched_noise`: the radial first-connectivity diagnostic and the before/at-closure comparison, including coordinates and diagrams.
- `noise`: the two noise models, the prescribed budget grid, exact event checks, both selectors and error summaries.
- `ellipse_checks`: seed-generated ellipse coordinates, analytic bounds and the exact connectivity/closure checks.
- `ellipse`: post-closure diagrams, arc-condition checks, nonlocal edge events and NN-crust results.
- `phase1`: the paper's bridge-free, Euclidean-baseline and grid-oracle comparisons. [The guide](phase1/README.md) explains the measurements.

CSV files contain measured rows or summaries. NPZ files contain arrays and persistence diagrams; load them with NumPy using `allow_pickle=False`. Checkpoint files let an interrupted calculation continue. The code guide defines pairing, numerical precision and selector conventions.

`generation_summary.json` lists the generated studies and their sizes. `verification_report.json` records the full release check. `file_manifest.json` lists SHA-256 hashes. When manuscript comparison is requested, the manuscript records identify the external source and PDF by their hashes. `generation_timing.json` and `runtime_benchmark.json` describe this run's execution times; clock times are not fixed scientific targets.

Counts distinguish distinct samples, model-indexed rows and method records. A repeated clean row across the two noise models is not an additional independent sample.
'''
FIGURES='''# Manuscript figures

Every figure is supplied as a PDF and a 300-dpi PNG. Regenerate all six pairs with `python3 code/build_figures.py` from the release folder, after rebuilding the summaries.

- `Fig1`: clean-study budget curves (Figure 2 in the paper).
- `Fig2`: radial first-connectivity diagnostic; supplied as a supporting figure.
- `noise-eph`: noisy diagram-error curves (Figure 3 in the paper).
- `noise-all-errors`: supporting curves for all three errors.
- `appendix_extra_k`: distributions of additional neighbours needed for closure.
- `unsupported-shortcut`: the example panel used to explain unsupported shortcuts.

The paper's opening schematic is drawn in its LaTeX source. Its other panels are not separate figure files in this folder.
'''
def main():
    for folder,text in [('code',CODE),('results',RESULTS),('manuscript_figures',FIGURES)]:
        (ROOT/folder).mkdir(exist_ok=True);(ROOT/folder/'README.md').write_text(text)
if __name__=='__main__':main()
