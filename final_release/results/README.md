# Results

The [release README](../README.md) maps these files to the manuscript. The folder also retains the additional sampling profiles, full peanut sweep and earlier finite checks. In the small-example report, `dense_shortcut` refers to the eight-point example in Remark 2. The six-point shortcut is recorded under `unsupported_shortcut` and `shortcut_decimal_certificate`.

The paired rank sweep stores arrays as `ranks[repetition, profile, method, (k_conn, k_close)]` and `largest_gap[repetition, profile]`. The labels `OR_kNN`, `mutual_kNN` and `R2_global` denote symmetric kNN, mutual kNN and chord-length ordering, respectively. The peanut checkpoint files retain sample coordinates and the diagrams at closure, first detection and twice the connecting radius. Smooth-curve exceptional checkpoints retain coordinates, gaps and the diagrams at the three enlargement factors. Load all NPZ files with `allow_pickle=False`.

The smooth-curve summary uses one million point sets per curve and size, totalling four million per curve and eight million across both curves. The three enlargement factors reuse these point sets.

`radial_noise_winding.csv` records spanning-tree fundamental-cycle checks for the 129 long-interval graphs per rule at amplitude 0.01. The cycle counts refer to graph cycles; they do not identify persistence-interval representatives.

The detailed data-generation protocols contain generation-script hashes. Where a script has separate generation and verification checksums, the former identifies the code used for the saved measurements and the latter identifies the supplied script used to replay them.
