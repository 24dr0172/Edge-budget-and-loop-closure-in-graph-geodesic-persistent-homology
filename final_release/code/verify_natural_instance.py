"""Certify the six-point ellipse companion using rational enclosing bounds.

The manuscript proof uses arc support, triangle counts and forced consecutive
vertices. Enumeration and floating-point persistence are corroboration only.
The ring closes at insertion 7 for the radius order and at insertion 6 for
the OR and mutual orders. No saved ellipse sweep, manuscript helper or extra
dependency is used.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import itertools
import json
import math
import sys
sys.dont_write_bytecode = True
import numpy as np
from finite_diagram import reduction_h1
from verify_arc_support import metric, betti_at
from verify_separation import cyclic_orientation

ROOT = Path(__file__).resolve().parents[1]
ANGLES = (30, 85, 150, 255, 290, 335)
PAIRS = list(itertools.combinations(range(6), 2))
RING = {(i, i+1) for i in range(5)} | {(0, 5)}
EDGES = sorted(RING | {(0, 4), (1, 4)})


class Interval:
    def __init__(self, low, high=None):
        self.low = F(low)
        self.high = self.low if high is None else F(high)
        assert self.low <= self.high

    def __add__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        return Interval(self.low+other.low, self.high+other.high)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.high, -self.low)

    def __sub__(self, other):
        return self + (-other if isinstance(other, Interval) else -F(other))

    def __mul__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        products = [a*b for a in (self.low, self.high) for b in (other.low, other.high)]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def square(self):
        low = 0 if self.low <= 0 <= self.high else min(self.low**2, self.high**2)
        return Interval(low, max(self.low**2, self.high**2))

    def sqrt(self):
        assert self.low >= 0
        scale = 10**18
        lower = math.isqrt(self.low.numerator*scale*scale//self.low.denominator)
        upper = math.isqrt(self.high.numerator*scale*scale//self.high.denominator)+1
        return Interval(F(lower, scale), F(upper, scale))

    def bounded(self):
        # Decimal endpoints are rounded outwards, not to nearest.
        scale = 10**12
        low = self.low.numerator*scale//self.low.denominator
        high = -((-self.high.numerator*scale)//self.high.denominator)
        return dict(lower=f'{low/scale:.12f}', upper=f'{high/scale:.12f}')


def atan_bounds(denominator, terms=25):
    x = F(1, denominator)
    partial = sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(terms))
    next_value = partial+(-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return Interval(min(partial, next_value), max(partial, next_value))


PI = 16*atan_bounds(5)-4*atan_bounds(239)


def trig(degrees, cosine=False):
    """Taylor polynomial plus a uniform remainder and input-error bound.

    Machin's identity gives pi. For a rational midpoint m and an exact
    angle in its interval, sine and cosine change by at most the input
    error, since their derivatives have absolute value at most one.
    Taylor's theorem bounds the omitted tail by |m|^(degree+1)/(degree+1)!.
    """
    x = PI*F(degrees, 180)
    m = (x.low+x.high)/2
    start = 0 if cosine else 1
    powers = [start+2*k for k in range(30)]
    partial = sum((-1)**k*m**power/F(math.factorial(power)) for k, power in enumerate(powers))
    degree = powers[-1]
    error = abs(m)**(degree+1)/F(math.factorial(degree+1))+(x.high-x.low)/2
    return Interval(partial-error, partial+error)


def enclosure_data():
    points = [(2*trig(a, True), trig(a)) for a in ANGLES]
    chords = [[Interval(0) for _ in ANGLES] for _ in ANGLES]
    for i, j in PAIRS:
        dx = points[i][0]-points[j][0]
        dy = points[i][1]-points[j][1]
        chords[i][j] = chords[j][i] = (dx.square()+dy.square()).sqrt()
    numeric = np.array([[float((v.low+v.high)/2) for v in row] for row in chords])
    order = sorted(PAIRS, key=lambda e: numeric[e])
    assert all(chords[a][b].high < chords[c][d].low
               for (a, b), (c, d) in zip(order, order[1:]))
    ranks = np.zeros((6, 6), int)
    for i in range(6):
        for k, j in enumerate(sorted((j for j in range(6) if j != i), key=lambda j: numeric[i, j]), 1):
            ranks[i, j] = k
    orders = dict(radius=order,
                  OR=sorted(PAIRS, key=lambda e: (min(ranks[e], ranks[e[::-1]]), numeric[e])),
                  mutual=sorted(PAIRS, key=lambda e: (max(ranks[e], ranks[e[::-1]]), numeric[e])))
    assert all(sorted(order[:8]) == EDGES for order in orders.values())
    closures = {name:1+max(order.index(e) for e in RING) for name, order in orders.items()}
    assert closures == dict(radius=7, OR=6, mutual=6)
    # Riemann integration is unnecessary: Cauchy--Schwarz bounds each
    # arc by sqrt(delta * integral(1+3 sin(theta)^2)); L >= 4 sqrt(5)
    # follows from the four inscribed quadrant chords.
    gap_upper = []
    for a, b in zip(ANGLES, ANGLES[1:]+(ANGLES[0]+360,)):
        delta = PI*F(b-a, 180)
        speed_squared_integral = F(5, 2)*delta-F(3, 4)*(trig(2*b)-trig(2*a))
        gap_upper.append((delta*speed_squared_integral).sqrt().high)
    perimeter_lower = 4*Interval(5).sqrt().low
    gap_bound = max(gap_upper)/perimeter_lower
    assert gap_bound < F(1, 3)
    # Enclose all graph distances by min-plus Floyd--Warshall. Each bound
    # encloses all coordinate correlations; present edges equal chords.
    distances = [[Interval(0) if i == j else None for j in range(6)] for i in range(6)]
    for i, j in EDGES:
        distances[i][j] = distances[j][i] = chords[i][j]
    for k in range(6):
        for i in range(6):
            for j in range(6):
                if distances[i][k] is None or distances[k][j] is None:
                    continue
                candidate = distances[i][k]+distances[k][j]
                current = distances[i][j]
                distances[i][j] = candidate if current is None else Interval(
                    min(current.low, candidate.low), min(current.high, candidate.high))
    for i, j in EDGES:
        distances[i][j] = distances[j][i] = chords[i][j]
    return chords, numeric, distances, orders, closures, gap_upper, gap_bound


def margin_certificate(c, d, gap_bound):
    """Rigorous lower bounds on the margins of the strict comparisons used.

    Structural comparisons are those the argument needs: chord orderings, support
    inequalities, the threshold memberships at 21/10 and the gap bound. Enclosure
    margins measure how far each quoted decimal bound is from the true value; they
    are small because the decimals are rounded.
    """
    ordered = sorted(PAIRS, key=lambda e: float(c[e[0]][e[1]].low))
    structural = {
        'closest_two_chords': min(c[b[0]][b[1]].low-c[a[0]][a[1]].high
                                  for a, b in zip(ordered, ordered[1:])),
        'c15_above_c14': c[1][5].low-c[1][4].high,
        'x0x4_step_below_chord': c[0][4].low-max(c[u][v].high for u, v in [(4, 5), (0, 5)]),
        'x1x4_step_below_chord': c[1][4].low-max(c[u][v].high for u, v in [(4, 5), (0, 5), (0, 1)]),
        'x1x4_fan_below_chord': c[1][4].low-max(d[4][0].high, d[4][5].high),
        'threshold_present_edges': F(21, 10)-max(d[i][j].high for i, j in EDGES),
        'threshold_absent_pairs': min(d[i][j].low for i, j in PAIRS if (i, j) not in EDGES)-F(21, 10),
        'gap_bound_below_one_third': F(1, 3)-gap_bound,
    }
    quoted = [((0, 4), '1.7807', '1.7808'), ((1, 4), '2.0018', '2.0019'),
              ((1, 2), '1.9698', '1.9699')]
    enclosure = {}
    for (i, j), low, high in quoted:
        enclosure[f'c{i}{j}_above_{low}'] = c[i][j].low-F(low)
        enclosure[f'c{i}{j}_below_{high}'] = F(high)-c[i][j].high
    enclosure['c15_above_2.1672'] = c[1][5].low-F('2.1672')
    enclosure['steps_below_1.242'] = F('1.242')-max(c[u][v].high for u, v in [(4, 5), (0, 5)])
    enclosure['steps_below_1.635'] = F('1.635')-max(c[u][v].high for u, v in [(4, 5), (0, 5), (0, 1)])
    enclosure['fan_below_1.781'] = F('1.781')-max(d[4][0].high, d[4][5].high)
    enclosure['absent_graph_distance_above_2.443'] = structural['threshold_absent_pairs']+F(21, 10)-F('2.443')
    assert min(structural.values()) > F(3, 100)
    assert min(enclosure.values()) > 0
    return dict(structural_lower_bounds={k: float(v) for k, v in structural.items()},
                enclosure_lower_bounds={k: float(v) for k, v in enclosure.items()},
                smallest_structural_margin=float(min(structural.values())),
                smallest_enclosure_margin=float(min(enclosure.values())),
                note='Structural margins exceed 0.03 raw units; rounded decimal enclosures can be about 2e-5 from the true value.')


def build_report():
    c, numeric, d, orders, closures, gap_upper, gap_bound = enclosure_data()
    assert max(gap_upper) < F('2.554')
    for edge, low, high in [((0,4), '1.7807', '1.7808'),
                            ((1,4), '2.0018', '2.0019'),
                            ((1,2), '1.9698', '1.9699')]:
        i,j=edge
        assert F(low) < c[i][j].low <= c[i][j].high < F(high)
    assert c[1][5].low > F('2.1672')
    assert max(c[u][v].high for u,v in [(4,5),(0,5)]) < F('1.242')
    assert max(c[u][v].high for u,v in [(4,5),(0,5),(0,1)]) < F('1.635')
    assert max(d[4][0].high,d[4][5].high) < F('1.781')
    assert min(d[i][j].low for i,j in PAIRS if (i,j) not in EDGES) > F('2.443')
    supports = []
    for edge, arc in [((0, 4), [4, 5, 0]), ((1, 4), [4, 5, 0, 1])]:
        a, b = edge
        outer = c[a][b]
        for u, v in zip(arc, arc[1:]):
            assert c[u][v].high < outer.low
        for v in arc:
            if v not in edge:
                assert d[4][v].high < outer.low
        supports.append(dict(edge=edge, arc=arc, chosen_endpoint=4,
                             max_step_upper=float(max(c[u][v].high for u, v in zip(arc, arc[1:]))),
                             interior_distance_upper=float(max(d[4][v].high for v in arc if v not in edge))))
    # The shorter ellipse arc from 4 to 1 passes through 5 and 0, since
    # the ellipse's speed is pi-periodic and the reverse angular span is 205 degrees.
    assert c[1][5].low > c[1][4].high
    numeric_d = metric(numeric, EDGES)
    assert all(float(d[i][j].low)-1e-13 <= numeric_d[i,j] <= float(d[i][j].high)+1e-13
               for i,j in PAIRS)
    expected = [sorted(RING | {(0,4)}), EDGES]
    threshold_records = []
    for t, want in zip([F(199,100), F(21,10)], expected):
        selected = []
        for i, j in PAIRS:
            value = d[i][j]
            assert value.high < t or value.low > t, 'Threshold needs more precision'
            if value.high < t:
                selected.append((i,j))
        assert selected == want
        triangles = [tri for tri in itertools.combinations(range(6),3)
                     if all(tuple(sorted(e)) in selected for e in itertools.combinations(tri,2))]
        permitted = [order for order in [(0,)+q for q in itertools.permutations(range(1,6))]
                     if cyclic_orientation(order, selected) is not None]
        fields = {str(p):int(betti_at(numeric_d,float(t),p)) for p in (2,3,5)}
        assert set(fields.values()) == {1}
        if t == F(199,100):
            assert triangles == [(0,4,5)]
            assert permitted == [(0,1,2,3,4,5),(0,5,4,3,2,1)]
        else:
            assert triangles == [(0,1,4),(0,4,5)] and not permitted
            # No triangle contains 12,23 or 34: these force one consecutive
            # block 1,2,3,4. Each direction of 14 then requires an absent edge.
            in_triangles = {e for tri in triangles for e in itertools.combinations(tri,2)}
            assert all(e not in in_triangles for e in [(1,2),(2,3),(3,4)])
            assert all(e not in selected for e in [(1,3),(2,4),(1,5)])
        threshold_records.append(dict(threshold=str(t),threshold_float=float(t),edges=selected,
                                      triangles=triangles,admissible_orders=permitted,
                                      betti_1_by_field=fields))
    margins = margin_certificate(c, d, gap_bound)
    bars, tiny = reduction_h1(numeric_d)
    assert bars.shape == (1,2) and tiny == 0
    assert abs(bars[0,0]-max(numeric[e] for e in RING)) < 1e-12
    assert bars[0,0] < 1.99 < 2.1 < bars[0,1]
    return dict(status='passed',angles_degrees=ANGLES,units='Raw ellipse (2 cos(theta), sin(theta)); divide by perimeter for length one.',
                graph_edges=EDGES,ring_closure_insertions=closures,common_prefix_edges=8,
                edge_orders=orders,arc_support=supports,arc_gap_raw_upper_bounds=[float(q) for q in gap_upper],
                normalized_gap_upper_bound=float(gap_bound),positive_reference_guaranteed=True,
                nonlocal_witness=dict(outer_chord=c[1][4].bounded(),endpoint_to_interior_chord=c[1][5].bounded(),
                                      minor_arc=[4,5,0,1]),
                chords={f'{i}{j}':c[i][j].bounded() for i,j in PAIRS},
                distances={f'{i}{j}':d[i][j].bounded() for i,j in PAIRS},
                margins=margins,
                proof_thresholds=threshold_records,positive_interval_raw=bars[0].tolist(),
                comparison_arithmetic='Rational Machin arctangent bounds for pi; Taylor remainders; integer outward square-root bounds; min-plus interval distances.',
                proof_scope='The paper gives all-field independence and an all-order impossibility argument; F2/F3/F5 calculations and enumeration corroborate them.',
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    report=build_report()
    target=args.report or ROOT/'results/natural_instance_verification.json'
    target.write_text(json.dumps(report,indent=2)+'\n')
    print('Natural ellipse companion passed: common eighth prefix; closures 7/6/6; arc support; live noncyclic threshold; enclosing arithmetic.')
