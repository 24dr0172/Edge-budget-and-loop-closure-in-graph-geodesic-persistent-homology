# Loop detection in graph-geodesic persistent homology

Code and numerical evidence for the paper by Shivam Singh, Annesha Sen and S. P. Tiwari. The current release is in [final_release/](final_release/).

- [Code and numerical protocols](final_release/code/) contain the generators, verification scripts, settings and pinned Python dependencies.
- [Results](final_release/results/) contain the saved rank sweeps, diagrams, radial-noise and winding checks, and smooth-curve experiments.
- [Manuscript figure](final_release/manuscript_figures/) contains the figure used by the paper.
- [Release guide](final_release/README.md) maps the saved evidence to the paper and gives reproduction commands.

The manuscript and Overleaf project are maintained separately and are not included in this repository.

## Verify the release

With Python 3.12, from the repository root:

```sh
python3 -m venv /path/to/geodesic-env
. /path/to/geodesic-env/bin/activate
python3 -m pip install -r final_release/code/requirements.txt
python3 final_release/tools/verify_release.py
```

The verifier checks the release checksums and saved numerical evidence; it does not check the manuscript or prove the theorems. For full reproduction instructions, see the [code guide](final_release/code/README.md). See also the [citation](CITATION.cff) and [license](LICENSE).
