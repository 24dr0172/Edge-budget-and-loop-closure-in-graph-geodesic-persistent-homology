"""Test whether sampled post-closure Rips graphs can be cyclic in the curve order.

Background (Adams--Coldren--Willmot). A directed graph on cyclically ordered
vertices is cyclic when every arc a->b implies a->w and w->b for every vertex w
strictly inside the forward cyclic arc from a to b. An undirected threshold graph
is cyclic in a given order when each of its edges can be oriented this way.

For a graph-geodesic metric D on n samples listed in curve (latent) order, and
for each unordered pair {a,b} with b=a+k (mod n), let

    innerF(a,b) = max over w strictly inside the forward arc a..b of
                  max(D(a,w), D(w,b)),
    innerR(a,b) = innerF(b,a)  (the complementary arc).

Convention: the Rips threshold graph at t has edge {x,y} iff D(x,y) < t.
Orienting {a,b} as a->b at t needs every w of the forward arc adjacent to both
endpoints, i.e. innerF(a,b) < t; similarly b->a needs innerR(a,b) < t. With the
tolerance TOL=1e-10, the pair is a *certified failure* at every t with

    D(a,b) + TOL < t <= min(innerF, innerR) - TOL,

where the edge is present and both orientations are blocked by a sample at
distance >= t from an endpoint. A certified failure proves that the threshold
graph is not cyclic in the curve order. Passing this local test is only a
necessary condition; such graphs are recorded as 'not certified failing', never
as cyclic. A weaker test that uses only the shorter arc is also recorded.

H1 bars are read from the saved release diagrams (no persistence is recomputed):
circle rows from results/useful_diagram/saved_prefixes.csv, ellipse rows from
results/ellipse/checkpoints/*.npz. Under the strict convention the longest bar is
alive for t in (birth, death].

Optionally (--two-sat) an exact 2-SAT orientation check (each edge receives
exactly one orientation) is run at a threshold grid for the arc-supported
single-bar ellipse graphs whose shorter-arc failure precedes their bar's death.

Outputs: results/curve_order/curve_order_checks.csv and summary.json.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.dont_write_bytecode = True

import numpy as np

from graph_methods import METHODS, build_orders, closure, graph_metric

ROOT = Path(__file__).resolve().parents[1]
TOL = 1e-10
COLUMNS = ['curve', 'seed', 'profile', 'method', 'kind', 'M', 'M_close', 'n',
           'arc_supported', 'bars', 'birth', 'death',
           'pairs', 'certified_pairs', 'certified_failure',
           'first_certified_threshold', 'first_certified_window_end',
           'last_certified_threshold', 'certified_before_death',
           'first_certified_before_death', 'certified_while_alive',
           'first_certified_while_alive', 'last_certified_while_alive',
           'certified_alive_windows', 'certified_alive_measure',
           'shorter_arc_failing_pairs', 'shorter_arc_failure',
           'shorter_arc_first_threshold', 'shorter_arc_before_death']


def inner_maxima(D):
    """Return F with F[a,b] = innerF(a,b) for the vertex order 0..n-1 (O(n^2)).

    F[a,(a+k)%n] = max(max_{0<j<k} D(a,a+j), max_{0<j<k} D(a+k-j, a+k));
    -inf when the arc has no interior vertex (k=1)."""
    n = len(D)
    idx = np.arange(n)
    fwd = D[idx[:, None], (idx[:, None] + idx[None, :]) % n]   # fwd[a,j]=D(a,a+j)
    bwd = D[idx[:, None], (idx[:, None] - idx[None, :]) % n]   # bwd[b,j]=D(b,b-j)
    A = np.full((n, n), -np.inf)
    B = np.full((n, n), -np.inf)
    A[:, 2:] = np.maximum.accumulate(fwd[:, 1:n - 1], axis=1)   # A[a,k], j<k
    B[:, 2:] = np.maximum.accumulate(bwd[:, 1:n - 1], axis=1)
    F = np.full((n, n), -np.inf)
    for k in range(2, n):
        b = (idx + k) % n
        F[idx, b] = np.maximum(A[:, k], B[b, k])
    return F


def merge(lo, hi):
    """Union of half-open intervals (lo, hi] as sorted disjoint [a, b] pairs."""
    out = []
    for a, b in sorted(zip(lo.tolist(), hi.tolist())):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def certified_at(D, t):
    """True when some pair is a certified two-direction failure at t."""
    F = inner_maxima(D)
    m = np.minimum(F, F.T)
    return bool(((D + TOL < t) & (t <= m - TOL)).any())


def analyse(D, s, birth, death):
    """Certified two-direction failures and the shorter-arc comparison."""
    n = len(D)
    F = inner_maxima(D)
    m = np.minimum(F, F.T)                     # min(innerF, innerR), symmetric
    iu, ju = np.triu_indices(n, 1)
    d, mm = D[iu, ju], m[iu, ju]
    lo, hi = d + TOL, mm - TOL                 # certified t in (lo, hi]
    cert = lo < hi
    out = dict(pairs=len(d), certified_pairs=int(cert.sum()), certified_failure=bool(cert.any()))
    out['first_certified_threshold'] = float(d[cert].min()) if cert.any() else None
    out['first_certified_window_end'] = float(mm[cert][np.argmin(d[cert])]) if cert.any() else None
    out['last_certified_threshold'] = float(hi[cert].max()) if cert.any() else None
    before = cert & (lo < np.minimum(hi, death))
    out['certified_before_death'] = bool(before.any())
    out['first_certified_before_death'] = float(d[before].min()) if before.any() else None
    alive = cert & (np.maximum(lo, birth) < np.minimum(hi, death))
    out['certified_while_alive'] = bool(alive.any())
    out['first_certified_while_alive'] = float(np.maximum(lo, birth)[alive].min()) if alive.any() else None
    out['last_certified_while_alive'] = float(np.minimum(hi, death)[alive].max()) if alive.any() else None
    windows = merge(np.maximum(lo, birth)[alive], np.minimum(hi, death)[alive])
    out['certified_alive_windows'] = len(windows)
    out['certified_alive_measure'] = float(sum(b - a for a, b in windows))
    # Weaker test: shorter latent arc only, from a to b=a+k with arc < 1/2.
    arc = (s[ju] - s[iu]) % 1.0
    inner_short = np.where(arc < .5, F[iu, ju], F[ju, iu])
    short_fail = inner_short - d > TOL
    out['shorter_arc_failing_pairs'] = int(short_fail.sum())
    out['shorter_arc_failure'] = bool(short_fail.any())
    out['shorter_arc_first_threshold'] = float(d[short_fail].min()) if short_fail.any() else None
    out['shorter_arc_before_death'] = bool(short_fail.any() and d[short_fail].min() < death)
    return out


def two_sat(D, t):
    """Exact check: can the graph {D<t} be oriented cyclically in index order?

    Literal 2e is 'edge e oriented i->j' (i<j, forward arc i+1..j-1), 2e+1 is
    'j->i' (arc j+1..n-1,0..i-1). Returns True when satisfiable."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    n = len(D)
    adj = D < t
    np.fill_diagonal(adj, False)
    iu, ju = np.nonzero(np.triu(adj, 1))
    eid = -np.ones((n, n), np.int64)
    eid[iu, ju] = np.arange(len(iu))
    lit = np.where(np.arange(n)[:, None] < np.arange(n)[None, :], 2 * eid, 2 * eid.T + 1)
    src, dst = [], []                  # lit[x,y]: literal for orientation x->y of {x,y}
    for e, (i, j) in enumerate(zip(iu, ju)):
        for here, arc in ((2 * e, np.arange(i + 1, j)),
                          (2 * e + 1, np.r_[np.arange(j + 1, n), np.arange(0, i)])):
            a, b = (i, j) if here % 2 == 0 else (j, i)
            if len(arc) == 0:
                continue
            if not (adj[a, arc].all() and adj[arc, b].all()):
                src.append(np.array([here])); dst.append(np.array([here ^ 1]))   # unit clause
                continue
            need = np.r_[lit[a, arc], lit[arc, b]]
            src += [np.full(len(need), here), need ^ 1]
            dst += [need, np.full(len(need), here ^ 1)]
    N = 2 * len(iu)
    src = np.concatenate(src) if src else np.zeros(0, np.int64)
    dst = np.concatenate(dst) if dst else np.zeros(0, np.int64)
    g = coo_matrix((np.ones(len(src), np.int8), (src, dst)), shape=(N, N)).tocsr()
    _, label = connected_components(g, directed=True, connection='strong')
    return bool((label[0::2] != label[1::2]).all())


def circle_tasks(root):
    data = root / 'results'
    cells = list(csv.DictReader(open(data / 'clean_closure/cells.csv')))
    prefixes = {(r['seed'], r['profile'], r['method'], int(r['M'])): r
                for r in csv.DictReader(open(data / 'useful_diagram/saved_prefixes.csv'))}
    grid = defaultdict(dict)
    for r in cells:
        if r['noise'] == 'clean':
            grid[(r['seed'], r['profile'], r['method'])][int(r['M'])] = r
    with np.load(data / 'clean_closure/latent_coordinates.npz') as z:
        for seed, profile in sorted({k[:2] for k in grid}):
            yield dict(curve='circle', seed=seed, profile=profile, s=z[f'{seed}_{profile}'],
                       grid={m: sorted(grid[(seed, profile, m)]) for m in METHODS},
                       prefixes={k: v for k, v in prefixes.items() if k[:2] == (seed, profile)})


def ellipse_tasks(root):
    data = root / 'results'
    rows = defaultdict(list)
    for r in csv.DictReader(open(data / 'ellipse/cells.csv')):
        rows[(r['seed'], r['profile'])].append(r)
    with np.load(data / 'ellipse_checks/coordinates.npz') as z:
        for (seed, profile), g in sorted(rows.items()):
            yield dict(curve='ellipse', seed=seed, profile=profile, s=z[f'{seed}_{profile}_latent'],
                       xy=z[f'{seed}_{profile}_xy'], rows=g,
                       diagrams=str(data / f'ellipse/checkpoints/{profile}_{seed}.npz'))


def longest_bar(bars):
    bars = np.asarray(bars, float)
    life = bars[:, 1] - bars[:, 0]
    i = int(np.argmax(life))
    return len(bars), float(bars[i, 0]), float(bars[i, 1])


def run_task(task):
    s = task['s']
    assert np.all(np.diff(s) > 0), 'samples must be stored in curve order'
    n = len(s)
    out = []
    if task['curve'] == 'circle':
        xy = np.column_stack((np.cos(2 * np.pi * s), np.sin(2 * np.pi * s))) / (2 * np.pi)
        _, iu, ju, L, orders, _ = build_orders(xy)
        for method in METHODS:
            o = orders[method]
            Mc = closure(o, iu, ju, n)
            graphs = [(Mc, 'closure')] + [(M, 'grid') for M in task['grid'][method] if M >= Mc]
            for M, kind in graphs:
                p = task['prefixes'][(task['seed'], task['profile'], method, M)]
                assert int(p['M_close']) == Mc
                bars, birth, death = int(p['bars']), float(p['birth']), float(p['death'])
                D = graph_metric(o, iu, ju, L, M, n)
                rec = dict(curve='circle', seed=task['seed'], profile=task['profile'], method=method,
                           kind=kind, M=M, M_close=Mc, n=n, arc_supported='', bars=bars,
                           birth=birth, death=death)
                rec.update(analyse(D, s, birth, death))
                out.append(rec)
    else:
        _, iu, ju, L, orders, _ = build_orders(task['xy'])
        with np.load(task['diagrams']) as z:
            assert np.array_equal(z['latent'], s) and np.array_equal(z['observed'], task['xy'])
            diagrams = {k: z[k] for k in z.files}
        for r in task['rows']:
            method, M = r['method'], int(r['M'])
            o = orders[method]
            assert closure(o, iu, ju, n) == int(r['M_close'])
            bars, birth, death = longest_bar(diagrams[f'{method}_M{M}'])
            assert bars == int(r['bars']) and abs((death - birth) - float(r['longest'])) < 1e-6
            kind = ('closure' if r['closure'] == 'True' else 'grid' if r['grid'] == 'True'
                    else 'first_nonlocal' if r['first_nonlocal'] == 'True' else 'other')
            D = graph_metric(o, iu, ju, L, M, n)
            rec = dict(curve='ellipse', seed=task['seed'], profile=task['profile'], method=method,
                       kind=kind, M=M, M_close=int(r['M_close']), n=n, arc_supported=r['arc_supported'],
                       bars=bars, birth=birth, death=death)
            rec.update(analyse(D, s, birth, death))
            out.append(rec)
    return out


def _two_sat_graph(job):
    root, r = job
    with np.load(root / 'results/ellipse_checks/coordinates.npz') as coords:
        xy = coords[f"{r['seed']}_{r['profile']}_xy"]
    _, iu, ju, L, orders, _ = build_orders(xy)
    D = graph_metric(orders[r['method']], iu, ju, L, r['M'], len(xy))
    first = r['first_certified_threshold'] if r['first_certified_threshold'] is not None else r['death']
    ts = [r['birth'] + i * (first - r['birth']) / 8 for i in range(1, 8)]
    ts += [first - 1e-6, first + 1e-9, (first + r['death']) / 2, r['death']]
    start = time.time()
    sat = [(float(t), two_sat(D, t)) for t in ts]
    local = [certified_at(D, t) for t in ts]
    assert all(not ok for (_, ok), c in zip(sat, local) if c), 'certified threshold must be unsatisfiable'
    below = [t for t, ok in sat if t < first and not ok]
    return dict(seed=r['seed'], profile=r['profile'], method=r['method'], M=r['M'],
                thresholds=[t for t, _ in sat], satisfiable=[ok for _, ok in sat], locally_certified=local,
                unsatisfiable_without_local_certificate=sum((not ok) and not c for (_, ok), c in zip(sat, local)),
                unsatisfiable_below_first_certified=below,
                unsatisfiable_from_first_certified_to_death=all(not ok for t, ok in sat if t > first),
                seconds=round(time.time() - start, 1))


def two_sat_scan(root, records, workers):
    """Exact 2-SAT at a threshold grid for ellipse graphs whose shorter-arc
    failure precedes their bar's death: seven points in (birth, first certified),
    just below and just above the first certified threshold, the midpoint to the
    death and the death."""
    targets = [r for r in records if r['curve'] == 'ellipse' and r['shorter_arc_before_death']
               and r['arc_supported'] == 'True' and r['bars'] == 1]
    with Pool(workers) as pool:
        return pool.map(_two_sat_graph, [(root, r) for r in targets])


def rng(values):
    values = [v for v in values if v is not None]
    return [min(values), max(values)] if values else None


def summarise(records):
    out = {}
    for curve in ('circle', 'ellipse'):
        g = [r for r in records if r['curve'] == curve]
        cert = [r for r in g if r['certified_failure']]
        before = [r for r in g if r['certified_before_death']]
        short = [r for r in g if r['shorter_arc_failure']]
        out[curve] = dict(
            graphs=len(g), certified_failure_graphs=len(cert),
            certified_failure_graphs_first_before_death=sum(r['first_certified_threshold'] < r['death'] for r in cert),
            certified_before_death_graphs=len(before),
            certified_while_alive_graphs=sum(r['certified_while_alive'] for r in g),
            first_certified_threshold_range=rng([r['first_certified_threshold'] for r in cert]),
            first_certified_before_death_range=rng([r['first_certified_before_death'] for r in before]),
            death_range=rng([r['death'] for r in g]),
            certified_alive_measure_fraction_range=rng([r['certified_alive_measure'] / (r['death'] - r['birth'])
                                                        for r in g if r['certified_while_alive']]),
            shorter_arc_failure_graphs=len(short),
            shorter_arc_first_threshold_range=rng([r['shorter_arc_first_threshold'] for r in short]),
            shorter_arc_before_death_graphs=sum(r['shorter_arc_before_death'] for r in g),
            shorter_arc_failing_but_not_certified=sum(r['shorter_arc_failure'] and not r['certified_failure'] for r in g),
            certified_but_not_shorter_arc_failing=sum(r['certified_failure'] and not r['shorter_arc_failure'] for r in g))
    keys = ('seed', 'profile', 'method', 'M', 'kind', 'M_close', 'arc_supported', 'bars', 'birth', 'death',
            'shorter_arc_first_threshold', 'certified_failure', 'first_certified_threshold',
            'certified_before_death', 'certified_while_alive', 'first_certified_while_alive',
            'last_certified_while_alive', 'certified_alive_windows', 'certified_alive_measure', 'certified_pairs')
    out['ellipse_shorter_arc_before_death'] = [
        {k: r[k] for k in keys} for r in records if r['curve'] == 'ellipse' and r['shorter_arc_before_death']]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--release', type=Path, default=ROOT)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--two-sat', action='store_true', help='Also run the exact 2-SAT scan (minutes).')
    args = ap.parse_args()
    root = args.release.resolve()
    start = time.time()
    tasks = list(circle_tasks(root)) + list(ellipse_tasks(root))
    with Pool(args.workers) as pool:
        records = [rec for part in pool.imap(run_task, tasks) for rec in part]
    target = root / 'results/curve_order'
    target.mkdir(exist_ok=True)
    with (target / 'curve_order_checks.csv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=COLUMNS, lineterminator='\n')
        w.writeheader()
        for r in records:
            w.writerow({k: ('' if r[k] is None else repr(r[k]) if isinstance(r[k], float) else r[k]) for k in COLUMNS})
    summary = dict(
        status='passed', script='code/verify_curve_order.py',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        tolerance=TOL,
        convention=('Rips threshold graph at t has edge {x,y} iff D(x,y) < t. Pair {a,b} is a certified '
                    'failure at t when D(a,b)+tol < t <= min(innerF,innerR)-tol. The longest saved H1 bar '
                    'is alive for t in (birth, death]. "before death": some certified t <= death; '
                    '"while alive": some certified t in (birth, death]. Graphs without a certified pair '
                    'are not certified failing; they are not claimed cyclic.'),
        diagram_sources={'circle': 'results/useful_diagram/saved_prefixes.csv',
                         'ellipse': 'results/ellipse/checkpoints/*.npz'},
        **summarise(records))
    if args.two_sat:
        summary['two_sat'] = dict(
            scope='Exact orientation-consistency check (one orientation per edge) for the arc-supported '
                  'single-bar ellipse graphs whose shorter-arc failure precedes the bar death: seven '
                  'thresholds in (birth, first certified), first certified -1e-6 and +1e-9, the midpoint '
                  'to the death, and the death.',
            graphs=two_sat_scan(root, records, max(1, min(args.workers, 3))))
    summary['seconds'] = round(time.time() - start, 1)
    (target / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: summary[k] for k in ('circle', 'ellipse')}, indent=2))


if __name__ == '__main__':
    main()
