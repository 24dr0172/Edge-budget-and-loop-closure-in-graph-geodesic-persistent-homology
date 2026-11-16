"""Paired circle rank sweep for the revised manuscript, with replayable seeds.

Run from any directory with NumPy and the released code available. This script
overwrites its saved rank-sweep results in this release. Each sample uses one
uniform sequence and phase across the three profiles; the two rules use the
same sampled coordinates.
"""
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

RELEASE = Path(__file__).resolve().parents[1]
HERE = RELEASE / 'results' / 'detection_rank_sweep'
sys.path.insert(0, str(RELEASE / 'code'))
from rank_sweep import fast_inverse, neighboring_keys, audit_local, wilson

PROFILES = ('uniform', 'strong_cosine', 'localized_trough')
METHODS = ('OR_kNN', 'mutual_kNN')
SETTINGS = ((128, 20000), (512, 20000), (4096, 20000), (32768, 5000))
SEED_BASE = 202611160000000


def sample(n, repetition):
    seed = SEED_BASE + n * 100000 + repetition
    rng = np.random.default_rng(seed)
    u = np.sort(rng.uniform(size=n))
    phase = float(rng.uniform())
    ranks = np.empty((3, 2, 2), dtype=np.int32)
    gaps = np.empty(3, dtype=float)
    for p, profile in enumerate(PROFILES):
        coords = fast_inverse(u, profile, phase)
        gaps[p] = np.max(np.diff(np.r_[coords, coords[0] + 1]))
        if gaps[p] >= .5:
            raise ValueError(f'sample outside local-rank assumptions: {seed}, {profile}')
        keys, _, _ = neighboring_keys(coords)
        for m, method in enumerate(METHODS):
            ranks[p, m] = np.sort(np.partition(keys[method], -2)[-2:])
    return ranks, gaps


def worker(setting):
    n, reps = setting
    started = time.monotonic()
    ranks = np.empty((reps, 3, 2, 2), dtype=np.int32)
    gaps = np.empty((reps, 3), dtype=float)
    for j in range(reps):
        ranks[j], gaps[j] = sample(n, j)
    np.savez_compressed(HERE / f'samples_{n}.npz', ranks=ranks, largest_gap=gaps)
    print(f'n={n}: {reps} paired samples, {time.monotonic()-started:.1f}s', flush=True)
    return n


def summarize():
    rows = []
    for n, reps in SETTINGS:
        with np.load(HERE / f'samples_{n}.npz') as f:
            ranks, gaps = f['ranks'], f['largest_gap']
        assert ranks.shape == (reps, 3, 2, 2) and gaps.shape == (reps, 3)
        assert np.all(ranks[..., 0] <= ranks[..., 1])
        # Replay selected records, across all profiles and rules.
        for j in (0, reps // 2, reps - 1):
            expected_ranks, expected_gaps = sample(n, j)
            np.testing.assert_array_equal(ranks[j], expected_ranks)
            np.testing.assert_array_equal(gaps[j], expected_gaps)
        for p, profile in enumerate(PROFILES):
            for m, method in enumerate(METHODS):
                extra = ranks[:, p, m, 1] - ranks[:, p, m, 0]
                count = int(np.count_nonzero(extra == 0))
                lo, hi = wilson(count, reps)
                rows.append(dict(n=n, samples=reps, profile=profile, method=method,
                    tied=count, tie_percent=100*count/reps,
                    wilson_lo_percent=100*lo, wilson_hi_percent=100*hi,
                    plus1=int(np.count_nonzero(extra <= 1)),
                    plus2=int(np.count_nonzero(extra <= 2)),
                    plus3=int(np.count_nonzero(extra <= 3)),
                    gap_at_least_third=int(np.count_nonzero(gaps[:, p] >= 1/3)),
                    max_gap=float(np.max(gaps[:, p]))))
    with (HERE / 'summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main():
    started = time.monotonic()
    audit = audit_local(29100000, 120)
    print('Local ranks versus full chord ranks:', json.dumps(audit), flush=True)
    with ProcessPoolExecutor(max_workers=4) as pool:
        assert set(pool.map(worker, SETTINGS)) == {n for n, _ in SETTINGS}
    rows = summarize()
    (HERE / 'protocol.json').write_text(json.dumps(dict(
        description='Circle local-rank sweep with paired density profiles and graph rules',
        settings=SETTINGS, profiles=PROFILES, methods=METHODS,
        seed_formula='202611160000000 + 100000*n + repetition (zero-based)',
        generator='NumPy default_rng; sorted uniform coordinates and independent phase within each stream',
        sample_file_layout='ranks[repetition,profile,method,(k_conn,k_close)]; largest_gap[repetition,profile]',
        local_vs_full_audit=audit, sample_replays_per_size=3,
        seconds=time.monotonic()-started, rows=rows), indent=2)+'\n')
    print('All summaries and seed replays complete.', flush=True)

if __name__ == '__main__':
    main()
