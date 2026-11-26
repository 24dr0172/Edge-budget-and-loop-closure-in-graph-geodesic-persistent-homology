# Code and numerical protocols

Run these commands from the release root after installing the pinned dependencies. The configuration summarises the seeds, sizes and tolerances. Study-specific protocols are under `results/`. Generators replace their outputs there. Run them in a copy to preserve the supplied records. [Dependencies](requirements.txt) and [configuration](config.json).

## Verification and sample replays

```sh
python3 tools/verify_release.py
python3 code/check_tied_ranks.py
python3 code/verify_intrinsic_circle.py
python3 code/polygon_certificate.py
python3 code/verify_small_examples.py --report /tmp/geodesic-small-examples.json
python3 code/clean_circle.py --seed 20310101 --profile uniform
python3 code/peanut_detection.py --seed 20361021 --alpha 0.4 --method R2_global --closure-only
python3 code/radius_enlargement_check.py --verify
python3 code/radius_enlargement_extra.py --verify
python3 code/smooth_curve_enlargement.py --verify
PYTHONPATH=code python3 -c 'from radial_noise_knn import evaluate; print(*evaluate(20320301), sep="\n")'
```

The first command checks the release hashes and saved evidence. The tied-rank check also evaluates the symmetric miss approximation at 512 and one billion points, as reported after Theorem 4. The small-example check includes the rational interval certificate for the eight-point shortcut in Remark 2. In the complete manuscript package, `python3 tools/verify_paper.py` also compares the manuscript and Overleaf files with the release.

`clean_circle.py --all` replays all 60 point sets. Remove `--closure-only` to replay the first non-empty prefix search for a peanut case. Use `--all` to replay all 300 method and sample combinations. The radial-noise command prints all eight graph records for one seed in the column order of the saved CSV. Both graph rules use the same noisy point set for each amplitude, and the seed reproduces the uniform draws on [−1, 1] used at every amplitude. [Saved radial-noise records](../results/radial_noise_knn.csv).

The million-sample circle runs at n = 4,096 and 32,768 replay the batch counts exactly. The replay allows an absolute tolerance of 10⁻¹² for the margins used to check the largest-gap condition. The supplied configuration records the generation and verification script hashes. [Larger circle runs](../results/radius_enlargement_extra/).

## Full study generation

```sh
python3 code/rank_table_sweep.py
g++ -O3 -std=c++17 code/uniform_rank_sweep.cpp -o /tmp/targeted-rank-sweep
python3 code/uniform_rank_check.py
python3 code/radial_noise_knn.py --output results/radial_noise_knn.csv --workers 4
python3 code/noise_winding_check.py
python3 code/check_noncircular_locality.py
python3 code/radius_enlargement_check.py
python3 code/radius_enlargement_extra.py
python3 code/smooth_curve_enlargement.py --workers 4
python3 code/check_exceptional.py
python3 code/check_interpolation.py
```

The two `--all` replay commands above recompute the clean-circle and peanut detection outcomes from their sample definitions. The peanut replay checks full diagrams at the selected prefixes and, for the radius ordering, after doubling the connecting radius. `check_exceptional.py` recomputes the 123 full diagrams for the 41 samples outside the grid-estimated good event. `check_interpolation.py` recomputes the arc-table comparisons. Both scripts write their own results under `results/smooth_curve_checks/`.

The winding check reconstructs the 258 graphs that meet the half-lifetime cutoff at noise amplitude 0.01. It records the fundamental graph cycles obtained from a spanning tree and replaces `results/radial_noise_winding.csv`.

For Table 3, the circle generator classifies diagram emptiness using the circle theorem and the ellipse and peanut generators evaluate the consecutive-chord event. The smooth-curve samples use independent draws for each curve and size. The full-diagram spot checks and the exceptional-sample checks compare the chord event with diagram emptiness.
