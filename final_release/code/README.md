# Code and numerical protocols

Run commands from the release root after installing [requirements.txt](requirements.txt). [config.json](config.json) summarizes the seeds, sizes and tolerances; the study-specific files under `results/` record the detailed generation protocols. Generators write into `results/` and overwrite their outputs. Copy the release before rerunning a complete cohort if the supplied records should be preserved.

## Verification and sample replays

```sh
python3 tools/verify_release.py
python3 code/check_tied_ranks.py
python3 code/noise_winding_check.py
python3 code/verify_intrinsic_circle.py
python3 code/polygon_certificate.py
python3 code/verify_small_examples.py --report /tmp/geodesic-small-examples.json
python3 code/clean_circle.py --seed 20310101 --profile uniform
python3 code/peanut_detection.py --seed 20361021 --alpha 0.4 --method R2_global --closure-only
python3 code/radius_enlargement_check.py --verify
python3 code/radius_enlargement_extra.py --verify
python3 code/smooth_curve_enlargement.py --verify
```

The first command checks the release hashes and saved evidence. In the complete manuscript package, `python3 tools/verify_paper.py` also compares the manuscript and Overleaf files with the release. `clean_circle.py --all` replays all 60 point sets. Remove `--closure-only` to replay the first-nonempty-prefix search for a peanut case; `--all` replays all 300 method/sample cases. These calculations require substantial time. `radial_noise_knn.evaluate(seed)` returns the eight paired graph records for one of the radial-noise seeds and can be compared with `results/radial_noise_knn.csv`.

The winding check reconstructs the 258 graphs with long intervals at noise amplitude 0.01, computes a spanning tree and its fundamental cycles, and replaces `results/radial_noise_winding.csv`. It checks graph cycles, not which cycles represent the persistence intervals.

## Full study generation

```sh
python3 code/rank_table_sweep.py
g++ -O3 -std=c++17 code/uniform_rank_sweep.cpp -o /tmp/targeted-rank-sweep
python3 code/uniform_rank_check.py
python3 code/radial_noise_knn.py --output results/radial_noise_knn.csv --workers 4
python3 code/check_noncircular_locality.py
python3 code/radius_enlargement_check.py
python3 code/radius_enlargement_extra.py
python3 code/smooth_curve_enlargement.py --workers 4
python3 code/check_exceptional.py
python3 code/check_interpolation.py
```

The two `--all` replay commands above recompute the clean-circle and peanut detection outcomes from their sample definitions. The peanut replay checks full diagrams at the selected prefixes and, for the radius ordering, after doubling the connecting radius. `check_exceptional.py` recomputes the 123 full diagrams for the 41 samples outside the grid-estimated good event. `check_interpolation.py` recomputes the arc-table comparisons. Both scripts write their own results under `results/smooth_curve_checks/`.

The radius-enlargement curve entries are chord-event proxies; the saved full-diagram checks and exceptional-point checks assess their agreement with diagram emptiness. A grid-estimated locality threshold is not a certified global bound. Numerical verification does not replace the manuscript's mathematical proofs.
