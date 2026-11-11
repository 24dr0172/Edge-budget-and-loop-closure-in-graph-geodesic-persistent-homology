"""Compute the two smallest-connected kNN diagrams in the paired radial study."""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import graph_methods as graph
from run_studies import finite_graph_diagram

N = 512
SEEDS = range(20320301, 20320701)
AMPLITUDES = (0., 0.0005, 0.002, 0.01)
METHODS = ('OR_kNN', 'mutual_kNN')
FIELDS = ('seed', 'delta', 'method', 'k', 'empty', 'longest_ratio',
          'at_least_quarter', 'at_least_half', 'all_consecutive')


def evaluate(seed):
    arc, clean = graph.sample_cloud(seed, 'uniform', n=N)
    reference = graph.diagram(graph.intrinsic(arc))
    assert len(reference) == 1
    reference_lifetime = float(reference[0, 1] - reference[0, 0])
    normal = clean * (2 * np.pi)
    draws = np.random.default_rng(seed + graph.NOISE_SEED_OFFSET).uniform(-1, 1, N)
    rows = []
    for delta in AMPLITUDES:
        points = clean + delta * draws[:, None] * normal
        _, ii, jj, lengths, orders, ranks = graph.build_orders(points)
        for method in METHODS:
            order = orders[method]
            first = graph.first_connected(order, ii, jj, N)
            ordered_ranks = ranks[method][order]
            k = int(ordered_ranks[first - 1])
            size = int(np.searchsorted(ordered_ranks, k, side='right'))
            covers = graph.closure(order, ii, jj, N) <= size
            metric = graph.graph_metric(order, ii, jj, lengths, size, N)
            diagram = finite_graph_diagram(metric, float(lengths[order[:size]].max()))
            if delta == 0:
                assert bool(len(diagram)) == covers
            longest = float(np.max(diagram[:, 1] - diagram[:, 0])) if len(diagram) else 0.
            ratio = longest / reference_lifetime
            rows.append((seed, delta, method, k, int(len(diagram) == 0), ratio,
                         int(ratio >= .25), int(ratio >= .5), int(covers)))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    temporary = args.output.with_suffix('.tmp')
    with temporary.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(FIELDS)
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            for completed, rows in enumerate(executor.map(evaluate, SEEDS), 1):
                writer.writerows(rows)
                if completed % 50 == 0:
                    print(f'{completed}/{len(SEEDS)} paired point sets', flush=True)
    temporary.replace(args.output)


if __name__ == '__main__':
    main()
