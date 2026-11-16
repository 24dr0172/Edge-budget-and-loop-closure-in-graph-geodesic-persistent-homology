"""Finite arc-support and homology checks used by the examples."""
import itertools
import numpy as np
from scipy.sparse.csgraph import shortest_path
from finite_diagram import reduction_h1
from graph_methods import intrinsic
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

    This accumulates triangle boundaries directly and tests the distinguished
    ring, without pairing persistence bars.
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


