"""Verify the arc-support extension using finite F2 algebra and saved ellipses.

The small graphs are deterministic validation samples generated here. Ellipse coordinates
and diagrams are reused; only their graph distances and support are recomputed.
No large-sample persistence calculation is performed. Geometry is float64.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True

import numpy as np
from scipy.sparse.csgraph import shortest_path

from finite_diagram import matching_distance, reduction_h1
from graph_methods import build_orders, closure, graph_metric, intrinsic

SUPPORT_TOL = 1e-10
CHECK_TOL = 1e-9


def metric(chords, edges):
    weights = np.full(chords.shape, np.inf)
    np.fill_diagonal(weights, 0.0)
    for a, b in edges:
        weights[a, b] = weights[b, a] = chords[a, b]
    return shortest_path(weights, directed=False)


def supporting_arc(distances, chords, a, b):
    """Return an arc and endpoint witness, or None, for the stated definition."""
    n = len(chords)
    cost = chords[a, b]
    for direction in (1, -1):
        steps = ((b - a) * direction) % n
        arc = (a + direction * np.arange(steps + 1)) % n
        if np.any(chords[arc[:-1], arc[1:]] > cost + SUPPORT_TOL):
            continue
        for endpoint in (a, b):
            if np.all(distances[endpoint, arc] <= cost + SUPPORT_TOL):
                return arc, endpoint
    return None


def supported(distances, chords, edges):
    return all(supporting_arc(distances, chords, a, b) is not None
               for a, b in edges)


def ring_death(distances, birth):
    """First threshold >= birth where the ring lies in the triangle span.

    Unlike the supplied ZIP's PH routine, this accumulates triangle boundaries
    directly and tests the distinguished ring, without pairing persistence bars.
    """
    n = len(distances)
    edges = list(itertools.combinations(range(n), 2))
    index = {e: k for k, e in enumerate(edges)}
    cycle = 0
    for a in range(n):
        cycle ^= 1 << index[tuple(sorted((a, (a + 1) % n)))]
    triangles = []
    for vertices in itertools.combinations(range(n), 3):
        pairs = list(itertools.combinations(vertices, 2))
        column = sum(1 << index[e] for e in pairs)
        triangles.append((max(distances[e] for e in pairs), column))
    triangles.sort()
    basis = {}

    def reduce(column):
        while column:
            pivot = column.bit_length() - 1
            if pivot not in basis:
                return column
            column ^= basis[pivot]
        return 0

    # Include B even if no triangle is born exactly at that threshold.
    thresholds = sorted({float(birth)} | {float(t) for t, _ in triangles if t >= birth})
    position = 0
    for threshold in thresholds:
        while position < len(triangles) and triangles[position][0] <= threshold:
            column = reduce(triangles[position][1])
            if column:
                basis[column.bit_length() - 1] = column
            position += 1
        if reduce(cycle) == 0:
            return threshold
    raise AssertionError('The complete Rips complex must fill the ring')


def prime_rank(matrix, prime):
    matrix = np.asarray(matrix, dtype=np.int64).copy() % prime
    row = 0
    for col in range(matrix.shape[1]):
        candidates = np.flatnonzero(matrix[row:, col])
        if not len(candidates):
            continue
        pivot = row + candidates[0]
        matrix[[row, pivot]] = matrix[[pivot, row]]
        matrix[row] = matrix[row] * pow(int(matrix[row, col]), -1, prime) % prime
        for other in range(row + 1, matrix.shape[0]):
            if matrix[other, col]:
                matrix[other] = (matrix[other] - matrix[other, col] * matrix[row]) % prime
        row += 1
        if row == matrix.shape[0]:
            break
    return row


def betti_at(distances, threshold, prime):
    """Independent signed triangle-boundary rank over an odd prime field."""
    n = len(distances)
    edges = [e for e in itertools.combinations(range(n), 2) if distances[e] <= threshold]
    index = {e: k for k, e in enumerate(edges)}
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            a = parent[a]
        return a

    for a, b in edges:
        parent[find(a)] = find(b)
    components = len({find(a) for a in range(n)})
    columns = []
    for a, b, c in itertools.combinations(range(n), 3):
        if all(e in index for e in ((a, b), (a, c), (b, c))):
            column = np.zeros(len(edges), dtype=np.int64)
            column[index[(b, c)]] = 1
            column[index[(a, c)]] = -1
            column[index[(a, b)]] = 1
            columns.append(column)
    boundary = np.column_stack(columns) if columns else np.empty((len(edges), 0), dtype=np.int64)
    return len(edges) - n + components - prime_rank(boundary, prime)


def curve_points(kind, s, rng):
    if kind == 'circle':
        return np.column_stack((np.cos(2*np.pi*s), np.sin(2*np.pi*s))) / (2*np.pi)
    theta = np.linspace(0., 2*np.pi, 40001)
    if kind == 'ellipse':
        points = np.column_stack((2*np.cos(theta), np.sin(theta)))
    else:
        if kind == 'peanut':
            radius = 1 + .6*np.cos(2*theta)
        else:
            a = rng.uniform(-.25, .25, 4)
            radius = 1 + a[0]*np.cos(2*theta) + a[1]*np.sin(3*theta) + a[2]*np.cos(4*theta) + a[3]*np.sin(5*theta)
        points = np.column_stack((radius*np.cos(theta), radius*np.sin(theta)))
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    perimeter = lengths.sum()
    cumulative = np.r_[0., np.cumsum(lengths)] / perimeter
    return np.column_stack([np.interp(s, cumulative, points[:,i]) for i in range(2)]) / perimeter


def check_small_graphs(samples):
    counts = dict(samples=0, budgets=0, supported_budgets=0, empty_supported=0,
                  euclidean_loop_supported=0, nested_supported_pairs=0,
                  small_gap_budgets=0, unsupported_multibar=0, odd_field_checks=0)
    for kind_index, kind in enumerate(('circle', 'ellipse', 'blob', 'peanut')):
        rng = np.random.default_rng(2026101201 + kind_index)
        accepted = 0
        while accepted < samples:
            n = int(rng.integers(6, 11))
            # Guarantee that part (c) is exercised, rather than merely counted.
            dense = kind == 'circle' and accepted < 4
            if dense:
                n = 12
                s = (np.arange(n) + rng.uniform(-.15, .15, n) + .5) / n
            else:
                s = np.sort(rng.random(n))
            reference = intrinsic(s)
            ref_bars, _ = reduction_h1(reference)
            if len(ref_bars) != 1 or np.diff(np.r_[s, s[0]+1]).max() >= .45:
                continue
            points = curve_points(kind, s, rng)
            chords, iu, ju, lengths, orders, _ = build_orders(points)
            ring = [tuple(sorted((a, (a+1)%n))) for a in range(n)]
            B = max(chords[e] for e in ring)
            tau = ring_death(chords, B)
            g, death_ref = ref_bars[0]
            diameter = reference.max()
            baseline = (death_ref-g) / (2*diameter)
            ring_bars, _ = reduction_h1(metric(chords, ring))
            ring_cap = min(matching_distance(ring_bars, ref_bars)/diameter, baseline)
            families = {name: [(int(iu[e]), int(ju[e])) for e in order]
                        for name, order in orders.items()}
            others = [e for e in itertools.combinations(range(n), 2) if e not in ring]
            families['random'] = ring + [others[i] for i in rng.permutation(len(others))]
            growing, remaining = list(ring), list(others)
            while remaining:
                D = metric(chords, growing)
                choices = [e for e in remaining if supporting_arc(D, chords, *e) is not None]
                if not choices:
                    break
                edge = choices[int(rng.integers(len(choices)))]
                growing.append(edge)
                remaining.remove(edge)
            families['supported_growth'] = growing
            for name, order in families.items():
                start = max(order.index(e) for e in ring) + 1
                previous = None
                for M in range(start, len(order)+1):
                    edges = order[:M]
                    D = metric(chords, edges)
                    bars, _ = reduction_h1(D)
                    counts['budgets'] += 1
                    if not supported(D, chords, edges):
                        counts['unsupported_multibar'] += len(bars) > 1
                        continue
                    counts['supported_budgets'] += 1
                    d = ring_death(D, B)
                    assert tau-CHECK_TOL <= d <= death_ref+CHECK_TOL
                    if d <= B + CHECK_TOL:
                        assert len(bars) == 0
                        counts['empty_supported'] += 1
                    else:
                        assert bars.shape == (1, 2) and np.max(np.abs(bars[0]-[B,d])) < CHECK_TOL
                    counts['euclidean_loop_supported'] += tau > B + CHECK_TOL
                    error = matching_distance(bars, ref_bars) / diameter
                    cap = min(error, baseline)
                    F = max(g-B, death_ref-d) / diameter
                    assert abs(cap-min(F, baseline)) < CHECK_TOL
                    assert cap >= ring_cap-CHECK_TOL
                    if previous is not None:
                        assert d <= previous[0]+CHECK_TOL and cap >= previous[1]-CHECK_TOL
                        counts['nested_supported_pairs'] += 1
                    previous = d, cap
                    if tau > (g+death_ref)/2+CHECK_TOL and B > (3*g-death_ref)/2+CHECK_TOL:
                        assert error < baseline
                        counts['small_gap_budgets'] += 1
                    if accepted == 0 and M in (start, len(order)):
                        # Interval interiors avoid numerical threshold ties.
                        for t in (B/2, (B+d)/2, (d+float(D.max()))/2):
                            expected = int(d > B+CHECK_TOL and B <= t < d)
                            for prime in (3, 5):
                                assert betti_at(D, t, prime) == expected
                                counts['odd_field_checks'] += 1
            accepted += 1
            counts['samples'] += 1
    assert counts['empty_supported'] and counts['small_gap_budgets'] and counts['unsupported_multibar']
    return {name: int(value) for name, value in counts.items()}


def check_witness():
    s = np.array([0., .051, .442, .491, .543, .953])
    ds = intrinsic(s)
    C = np.sin(np.pi*ds)/np.pi
    ring = [tuple(sorted((a, (a+1)%6))) for a in range(6)]
    D = metric(C, ring+[(0,3)])
    bars, _ = reduction_h1(D)
    assert bars.shape == (2,2)
    assert supporting_arc(D, C, 0, 3) is None
    assert all(betti_at(D, .33, p) == 2 for p in (2,3,5))
    return dict(s=s.tolist(), reference=reduction_h1(ds)[0].tolist(),
                ring_bars=reduction_h1(metric(C, ring))[0].tolist(),
                extra_edge_bars=bars.tolist(), extra_edge_cost=float(C[0,3]),
                euclidean_ring_death=float(ring_death(C,max(C[e] for e in ring))),
                beta1_at_point33_over_F2_F3_F5=2)


def check_ellipses(release):
    folder = release/'results/ellipse_checks'
    rows = list(csv.DictReader((folder/'checks.csv').open()))
    with np.load(folder/'coordinates.npz') as coords, np.load(folder/'diagrams.npz') as diagrams:
        cache = {}
        for row in rows:
            key = f"{row['seed']}_{row['profile']}"
            if key not in cache:
                cache[key] = build_orders(coords[key+'_xy'])
            C, iu, ju, lengths, orders, _ = cache[key]
            order = orders[row['method']]
            M = int(row['M_close'])
            assert closure(order, iu, ju, len(C)) == M
            D = graph_metric(order, iu, ju, lengths, M, len(C))
            edges = [(int(iu[e]),int(ju[e])) for e in order[:M]]
            assert supported(D, C, edges)
            assert max(C[e] for e in edges) < .18
            bars = diagrams[f"{key}_{row['method']}_M{M}"]
            assert bars.shape == (1,2)
            assert abs(bars[0,0]-float(row['B'])) <= float(row['endpoint_tolerance'])
    return dict(coordinate_sets=len(cache), closure_graphs=len(rows),
                arc_supported=len(rows), saved_single_bar_diagrams=len(rows),
                scope='Distances recomputed at closure; persistence diagrams reused, no ellipse grid or new ellipse persistence.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--release', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--samples-per-curve', type=int, default=20)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    assert args.samples_per_curve >= 4
    report = dict(date='2026-10-12', status='passed',
                  seed_base=2026101201, samples_per_curve=args.samples_per_curve,
                  support_tolerance=SUPPORT_TOL, comparison_tolerance=CHECK_TOL,
                  small_graphs=check_small_graphs(args.samples_per_curve),
                  witness=check_witness(), ellipses=check_ellipses(args.release.resolve()))
    report['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.report:
        args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
