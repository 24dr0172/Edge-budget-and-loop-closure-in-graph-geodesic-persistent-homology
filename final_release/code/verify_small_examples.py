"""Corroborate the small examples of the paper numerically.

The ellipse companion uses six fixed coordinates; no saved ellipse sweep,
noise, bridge, probability or large-rank output is read.
Numerical checks do not replace the analytical, all-field manuscript proofs.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import sys
sys.dont_write_bytecode = True
import numpy as np
from finite_diagram import reduction_h1, POSITIVE_TOL
from verify_arc_support import metric, supported, supporting_arc, betti_at, check_witness
from verify_separation import build_report as separation_report
from verify_natural_instance import build_report as natural_report

ROOT = Path(__file__).resolve().parents[1]


def circle(s):
    s = np.asarray(s, float)
    d = np.abs(s[:, None]-s[None, :])
    intrinsic = np.minimum(d, 1-d)
    return intrinsic, np.sin(np.pi*intrinsic)/np.pi


def ring(n):
    return sorted({tuple(sorted((i, (i+1) % n))) for i in range(n)})


def check_hexagon():
    _, chords = circle(np.arange(6)/6)
    path = [(i, i+1) for i in range(5)]
    edges = path+[(1, 5)]
    d = metric(chords, edges)
    bars, tiny = reduction_h1(d)
    a, b = 1/(2*np.pi), np.sqrt(3)/(2*np.pi)
    assert np.allclose(bars, [[b, 2*a]], rtol=0, atol=1e-12) and tiny == 0
    assert (0, 5) not in edges
    # The unique shorter arc 1,0,5 crosses the missing ring edge 0,5.
    assert 2/6 < 1/2 and (0, 1) in edges and (0, 5) not in edges
    threshold = (b+2*a)/2
    present = {e for e in itertools.combinations(range(6), 2) if d[e] <= threshold}
    assert present == set(edges)
    assert all(betti_at(d, threshold, p) == 1 and
               betti_at(d, 2*a+1e-12, p) == 0 for p in (2, 3, 5))
    old = metric(chords, path)
    assert abs(old[1, 5]-4*a) < 1e-12
    closed = metric(chords, ring(6))
    assert abs(closed[1, 5]-2*a) < 1e-12
    return dict(status='passed',missing_ring_edge=[0, 5],added_edge=[1, 5],
                short_hops_first=False,bar=bars[0].tolist(),
                betti_fields=[2, 3, 5],old_open_path_distance=float(old[1, 5]),
                old_closed_ring_distance=float(closed[1, 5]),
                scope='Specified arbitrary prefix; not a standard radius/kNN prefix.')


def check_ties():
    # Four-point reference and all endpoint ranks, independently of local counts.
    d, chords = circle([0, .34, .561, .778])
    ranks = np.argsort(np.argsort(np.where(np.eye(4), np.inf, chords), axis=1), axis=1)+1
    keys = {}
    for method, combine in [('OR', min), ('mutual', max)]:
        values = sorted(combine(int(ranks[i, j]), int(ranks[j, i])) for i, j in ring(4))
        keys[method] = values
    assert keys == dict(OR=[1, 1, 1, 2], mutual=[1, 2, 2, 2])
    bars, tiny = reduction_h1(d)
    assert np.allclose(bars, [[.34, .438]], atol=1e-12, rtol=0) and tiny == 0
    # One regular polygon and one non-square rectangle corroborate the tie mechanisms.
    _, regular = circle(np.arange(6)/6)
    regular_d = metric(regular, ring(6))
    assert betti_at(regular_d, 1.1/(2*np.pi), 3) == 1
    xy = np.array([(2., 1.), (-2., 1.), (-2., -1.), (2., -1.)])
    rectangle = np.linalg.norm(xy[:, None]-xy[None, :], axis=2)
    assert all(len(set(rectangle[i, j] for j in range(4) if i != j)) == 3 for i in range(4))
    assert {e for e in itertools.combinations(range(4), 2) if rectangle[e] <= 4} == set(ring(4))
    rectangle_bars, _ = reduction_h1(metric(rectangle, ring(4)))
    assert np.allclose(rectangle_bars, [[4., 6.]], atol=1e-12, rtol=0)
    return dict(status='passed',four_point_sorted_keys=keys,
                four_point_reference=bars.tolist(),rectangle_bar=rectangle_bars.tolist(),
                scope='Finite corroborating cases; general polygon/rectangle arguments are in Appendix B.')


def shortcut_decimal_certificate():
    """Rational sine and min-plus bounds for Remark 2, not float comparisons."""
    from fractions import Fraction as Q
    from verify_natural_instance import Interval, PI, trig
    arc=list(map(Q,['0','0.051','0.442','0.491','0.543','0.953']))
    invpi=Interval(1/PI.high,1/PI.low)
    chord=[[Interval(0) for _ in arc] for _ in arc]
    for i,j in itertools.combinations(range(6),2):
        gap=abs(arc[i]-arc[j]);gap=min(gap,1-gap)
        chord[i][j]=chord[j][i]=trig(180*gap)*invpi
    distances=[[Interval(0) if i==j else None for j in arc] for i in arc]
    for i,j in ring(6)+[(0,3)]:distances[i][j]=distances[j][i]=chord[i][j]
    for k in range(6):
        for i in range(6):
            for j in range(6):
                if distances[i][k] is None or distances[k][j] is None:continue
                candidate=distances[i][k]+distances[k][j];old=distances[i][j]
                distances[i][j]=candidate if old is None else Interval(min(old.low,candidate.low),min(old.high,candidate.high))
    bounds={(0,2):Q('0.35060'),(0,4):Q('0.35249'),(3,1):Q('0.34862'),(3,5):Q('0.35742')}
    margins=[distances[i][j].low-bound for (i,j),bound in bounds.items()]
    assert min(margins)>Q('0.00001')
    assert Q('0.31818')<chord[0][3].low<chord[0][3].high<Q('0.31819')
    return dict(status='passed',method='Machin pi, rational sine Taylor bounds, min-plus path bounds',
                shortcut_cost=chord[0][3].bounded(),
                distances={f'{i+1}{j+1}':distances[i][j].bounded() for i,j in bounds},
                lower_bound_margin=float(min(margins)))


def build_report():
    separation = separation_report()
    assert not separation['common_orders_at_live_thresholds']
    witness = check_witness()
    return dict(status='passed',scope='Small examples of the paper only.',
                fields_checked=[2, 3, 5],float64_comparison_tolerance=1e-12,
                independent_reduction_positive_tolerance=POSITIVE_TOL,
                separation=dict(proof_thresholds=separation['proof_thresholds'],
                    common_orders=separation['common_orders_at_live_thresholds'],
                    positive_interval_radius_one=separation['positive_interval_radius_one']),
                hexagon=check_hexagon(),unsupported_shortcut=witness,ties=check_ties(),
                shortcut_decimal_certificate=shortcut_decimal_certificate(),
                natural_ellipse=natural_report(),
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', type=Path, default=ROOT)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    report = build_report()
    target = args.report or args.release/'results/small_examples_verification.json'
    target.write_text(json.dumps(report, indent=2)+'\n')
    print('Small examples passed: live separation, natural ellipse, hexagon, unsupported shortcut and tie cases.')
