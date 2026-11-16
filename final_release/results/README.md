# Results

The release README maps these files to the manuscript. All point sets, summary counts, saved intervals and numerical checks in this folder support the current paper.

The paired rank sweep stores arrays as `ranks[repetition, profile, method, (k_conn, k_close)]` and `largest_gap[repetition, profile]`. The peanut checkpoint files retain sample coordinates and the diagrams used for the closure, first-detection and doubled-radius claims. Smooth-curve exceptional checkpoints retain coordinates, gaps and the three returned diagrams. Load all NPZ files with `allow_pickle=False`.

`radial_noise_winding.csv` records spanning-tree fundamental-cycle checks for the 129 long-interval graphs per rule at amplitude 0.01. The cycle counts refer to graph cycles; they do not identify persistence-interval representatives.

The detailed data-generation protocols contain generation-script hashes. Where a script has separate generation and verification checksums, the former identifies the code used for the saved measurements and the latter identifies the supplied script used to replay them.
