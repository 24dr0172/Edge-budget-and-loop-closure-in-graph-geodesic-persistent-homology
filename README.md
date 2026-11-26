# Loop detection in graph-geodesic persistent homology

Code and numerical evidence for the 32-page manuscript by Shivam Singh, Annesha Sen and S. P. Tiwari. The [v1.0.2 release](https://github.com/24dr0172/Edge-budget-and-loop-closure-in-graph-geodesic-persistent-homology/releases/tag/v1.0.2) accompanies this version.

- [Release guide](final_release/README.md) maps the manuscript tables and numerical checks to the saved evidence.
- [Code and protocols](final_release/code/) contain the generators, verification scripts, experiment settings and pinned Python dependencies.
- [Results](final_release/results/) contain the saved rank sweeps, diagrams, radial-noise and winding checks, and smooth-curve experiments.
- [Manuscript record](final_release/manuscript_record.json) identifies the LaTeX source, compiled PDF and figure by SHA-256. The manuscript and Overleaf files are distributed in the companion manuscript package.

## Verify the archive

Run these commands from the repository root with Python 3.12.

```sh
python3.12 -m venv .geodesic-env
. .geodesic-env/bin/activate
python3 -m pip install -r final_release/code/requirements.txt
python3 final_release/tools/verify_release.py
```

The verifier checks the release manifest and saved numerical evidence, including selected seed replays. [Reproduction commands](final_release/code/README.md) · [Citation](CITATION.cff) · [License](LICENSE)
