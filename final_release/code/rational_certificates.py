"""Exact rational certificates for two constructed facts used with the smooth-curve results.

1. The polygon example: a connected radius graph on a simple closed
   six-vertex polygon, formed before ring closure, whose graph-distance diagram is
   nonempty. All squared chords are rational and compared exactly. Graph
   distances are sums of square roots of rationals; they are enclosed by
   integer square roots at a fixed scale, and each one is identified with an
   exact symbolic value. Equal symbolic values are exact ties; all other
   comparisons are certified by disjoint intervals. Betti numbers are exact
   ranks over Q, F2, F3 and F5.
2. The uniform gap threshold g0: the unique root in [0, 1/6] of
   h(g) = (2/pi) sin(pi (1/3 + g)) - 1/3 - 2g. Pi is enclosed by Machin's
   formula, and sine by a Taylor sum with an alternating tail. A rational
   bracket of width 1e-13 is certified, together with the margins used by the
   strict-optimality condition.
No floating-point value decides any certified comparison.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import argparse, hashlib, itertools, json

ROOT = Path(__file__).resolve().parents[1]
FOLDER = 'results/rational_certificates'
SCALE = 10**40
POINTS = [(F(0), F(0)), (F(1), F(0)), (F(1), F(1)), (F(0), F(9, 10)), (F(-1, 10), F(101, 100)), (F(-1, 5), F(9, 5))]
G0_BRACKET = (F(15127097021660, 10**14), F(15127097021670, 10**14))


# ---------- interval helpers ----------
def sqrt_interval(q):
    """Integer-scaled enclosure [lo, hi]/SCALE of sqrt(q) for a nonnegative rational q."""
    num = q.numerator*SCALE*SCALE
    lo = isqrt(num//q.denominator)
    hi = lo+1
    assert F(lo, SCALE)**2 <= q <= F(hi, SCALE)**2
    return F(lo, SCALE), F(hi, SCALE)


def squarefree(q):
    """Write sqrt(q) = c*sqrt(m) with rational c and square-free integer m >= 1."""
    a = q.numerator*q.denominator
    coefficient, m, p = 1, 1, 2
    while p*p <= a:
        while a % (p*p) == 0:
            a //= p*p; coefficient *= p
        if a % p == 0:
            a //= p; m *= p
        p += 1
    m *= a
    return F(coefficient, q.denominator), m


def symbolic(squares):
    """Exact normal form of a sum of square roots: {square-free radicand: rational coefficient}."""
    form = {}
    for q in squares:
        c, m = squarefree(q)
        form[m] = form.get(m, F(0))+c
    return tuple(sorted((m, c) for m, c in form.items() if c))


def form_text(form):
    return ' + '.join(str(c) if m == 1 else f'{c}*sqrt({m})' for m, c in form)


# ---------- exact linear algebra ----------
def rank(rows, prime=None):
    rows = [list(r) for r in rows]
    if prime is not None:
        rows = [[x % prime for x in r] for r in rows]
    r, cols = 0, len(rows[0]) if rows else 0
    for c in range(cols):
        pivot = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                if prime is None:
                    f = F(rows[i][c])/rows[r][c]
                    rows[i] = [a-f*b for a, b in zip(rows[i], rows[r])]
                else:
                    f = rows[i][c]*pow(rows[r][c], -1, prime) % prime
                    rows[i] = [(a-f*b) % prime for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


def betti1(n, edges, triangles, prime=None):
    """First Betti number of the flag complex with these edges and triangles."""
    index = {e: k for k, e in enumerate(edges)}
    d1 = [[0]*len(edges) for _ in range(n)]
    for k, (a, b) in enumerate(edges):
        d1[a][k] -= 1; d1[b][k] += 1
    d2 = [[0]*len(triangles) for _ in edges]
    for k, (a, b, c) in enumerate(triangles):
        d2[index[(b, c)]][k] += 1; d2[index[(a, c)]][k] -= 1; d2[index[(a, b)]][k] += 1
    r1 = rank(d1, prime) if edges else 0
    r2 = rank(d2, prime) if triangles else 0
    return len(edges)-r1-r2


def orientation(p, q, r):
    return (q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0])


def segments_cross(p1, p2, q1, q2):
    """Closed segments intersect (exact rational test, including touching)."""
    def on(p, q, r):
        return min(p[0], r[0]) <= q[0] <= max(p[0], r[0]) and min(p[1], r[1]) <= q[1] <= max(p[1], r[1])
    o = [orientation(p1, p2, q1), orientation(p1, p2, q2), orientation(q1, q2, p1), orientation(q1, q2, p2)]
    if (o[0] > 0) != (o[1] > 0) and (o[2] > 0) != (o[3] > 0) and all(o):
        return True
    return any(v == 0 and on(*t) for v, t in zip(o, [(p1, q1, p2), (p1, q2, p2), (q1, p1, q2), (q1, p2, q2)]))


def six_point_certificate():
    n = len(POINTS)
    pairs = list(itertools.combinations(range(n), 2))
    sq = {(i, j): (POINTS[i][0]-POINTS[j][0])**2+(POINTS[i][1]-POINTS[j][1])**2 for i, j in pairs}
    ring = [(i, i+1) for i in range(n-1)]+[(0, n-1)]
    r2 = sq[(2, 3)]                       # r = c(3,4) in the 1-based labels of the text
    assert r2 == F(101, 100)
    graph = sorted(e for e in pairs if sq[e] <= r2)
    expected = sorted([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 3), (3, 5)])
    assert graph == expected and (0, n-1) not in graph and sq[(0, n-1)] == F(328, 100) > r2
    smallest_excluded = min(sq[e] for e in pairs if e not in graph)
    # Connectivity of the radius graph.
    seen, stack = {0}, [0]
    while stack:
        a = stack.pop()
        for b in range(n):
            if b not in seen and tuple(sorted((a, b))) in graph:
                seen.add(b); stack.append(b)
    assert len(seen) == n
    # Simple polygon, longest segment below half the perimeter, all arc gaps below one half.
    nonadjacent = [(e, f) for e, f in itertools.combinations(ring, 2) if not set(e) & set(f)]
    assert not any(segments_cross(POINTS[e[0]], POINTS[e[1]], POINTS[f[0]], POINTS[f[1]]) for e, f in nonadjacent)
    seg = {e: sqrt_interval(sq[e]) for e in ring}
    perimeter = (sum(v[0] for v in seg.values()), sum(v[1] for v in seg.values()))
    longest = max(ring, key=lambda e: sq[e])
    assert longest == (0, n-1) and 2*seg[longest][1] < perimeter[0]
    # Short-hops-first fails: the shortcut {1,4} has shorter arc 4-5-6-1, which uses the missing edge.
    arc_456_1 = [(3, 4), (4, 5), (0, 5)]
    arc_1234 = [(0, 1), (1, 2), (2, 3)]
    lo_a, hi_a = sum(seg[e][0] for e in arc_456_1), sum(seg[e][1] for e in arc_456_1)
    lo_b = sum(seg[e][0] for e in arc_1234)
    assert hi_a < lo_b and (0, n-1) not in graph
    # Graph distances: enumerate every simple path, enclose its length and keep its exact form.
    weight = {e: sqrt_interval(sq[e]) for e in graph}
    adjacency = {v: [b for b in range(n) if tuple(sorted((v, b))) in graph] for v in range(n)}
    distance = {}
    for a, b in pairs:
        paths, stack = [], [(a, [a])]
        while stack:
            v, path = stack.pop()
            if v == b:
                paths.append(path); continue
            for w in adjacency[v]:
                if w not in path:
                    stack.append((w, path+[w]))
        intervals = []
        for path in paths:
            edges = [tuple(sorted(p)) for p in zip(path, path[1:])]
            intervals.append((sum(weight[e][0] for e in edges), sum(weight[e][1] for e in edges),
                              symbolic([sq[e] for e in edges]), path))
        intervals.sort(key=lambda v: v[0])
        best = intervals[0]
        ties = [v for v in intervals if v[2] == best[2]]
        others = [v for v in intervals if v[2] != best[2]]
        assert all(best[1] < v[0] for v in others), (a, b)
        distance[(a, b)] = dict(lo=best[0], hi=best[1], form=best[2], path=[p+1 for p in best[3]], paths=len(intervals),
                                tied_paths=len(ties))
    # Certified order of distinct values; identical exact forms are exact ties.
    values = sorted({d['form'] for d in distance.values()}, key=lambda f: next(d['lo'] for d in distance.values() if d['form'] == f))
    enclosure = {f: next((d['lo'], d['hi']) for d in distance.values() if d['form'] == f) for f in values}
    assert all(enclosure[x][1] < enclosure[y][0] for x, y in zip(values, values[1:]))
    level = {e: values.index(distance[e]['form']) for e in pairs}
    r_form = symbolic([r2])
    birth_level = values.index(r_form)
    assert distance[(2, 3)]['form'] == r_form
    death_form = symbolic([F(1), F(81, 100)])                 # 1 + 9/10 = 19/10
    death_level = values.index(death_form)
    assert death_form == ((1, F(19, 10)),) and distance[(1, 3)]['form'] == death_form
    betti = {}
    for k in range(len(values)):
        edges = [e for e in pairs if level[e] <= k]
        triangles = [t for t in itertools.combinations(range(n), 3) if all(tuple(sorted(p)) in edges for p in itertools.combinations(t, 2))]
        betti[k] = {name: betti1(n, edges, triangles, prime) for name, prime in [('Q', None), ('F2', 2), ('F3', 3), ('F5', 5)]}
        assert len(set(betti[k].values())) == 1
    one = [k for k in betti if betti[k]['Q'] == 1]
    assert one == list(range(birth_level, death_level)) and all(betti[k]['Q'] == 0 for k in betti if k not in one)
    at_r = [e for e in pairs if level[e] <= birth_level]
    triangles_at_r = [t for t in itertools.combinations(range(n), 3) if all(tuple(sorted(p)) in at_r for p in itertools.combinations(t, 2))]
    assert sorted(at_r) == expected and triangles_at_r == [(3, 4, 5)]
    label = lambda e: f'{e[0]+1}{e[1]+1}'
    return dict(status='passed', points=[[str(x), str(y)] for x, y in POINTS], labels='1-based, as in the text',
        radius_squared=str(r2), radius_enclosure=[str(v) for v in sqrt_interval(r2)],
        radius_graph_edges=[label(e) for e in graph], missing_ring_edge='16', missing_ring_edge_squared=str(sq[(0, n-1)]),
        smallest_excluded_squared_chord=str(smallest_excluded), connected=True, simple_polygon=True,
        longest_segment='16', longest_segment_below_half_perimeter=True,
        shortcut='14', shortcut_shorter_arc=[4, 5, 6, 1], shortcut_arc_uses_missing_edge=True, short_hops_first=False,
        rips_edges_at_r=[label(e) for e in at_r], triangles_at_r=['456'],
        betti1_at_r={k: v for k, v in betti[birth_level].items()},
        bar=dict(birth='sqrt(101/100)', death='19/10', birth_level=birth_level, death_level=death_level,
                 fields=['Q', 'F2', 'F3', 'F5']),
        distinct_distance_values=len(values),
        smallest_adjacent_separation=str(min(enclosure[y][0]-enclosure[x][1] for x, y in zip(values, values[1:]))),
        distances={label(e): dict(value=form_text(d['form']), path=d['path'], lower=f"{float(d['lo']):.12f}",
                                  upper=f"{float(d['hi']):.12f}") for e, d in distance.items()},
        scope='The polygon only. Smoothing it to a C2 curve of length one is not certified here.')


# ---------- the threshold g0 ----------
def arctan_inverse(k, terms=60):
    """Enclosure of arctan(1/k) by the alternating Taylor series."""
    x = F(1, k)
    total = F(0)
    for j in range(terms):
        total += (-1)**j*x**(2*j+1)/(2*j+1)
    tail = x**(2*terms+1)/(2*terms+1)
    return (total, total+tail) if terms % 2 == 0 else (total-tail, total)


def pi_interval():
    a_lo, a_hi = arctan_inverse(5)
    b_lo, b_hi = arctan_inverse(239)
    return 16*a_lo-4*b_hi, 16*a_hi-4*b_lo


def sin_interval(x, terms=40):
    """Enclosure of sin(x) for rational 0 <= x <= 2."""
    assert 0 <= x <= 2
    term, total = x, x
    for k in range(1, terms):
        term = -term*x*x/((2*k)*(2*k+1)); total += term
    tail = abs(term*x*x/((2*terms)*(2*terms+1)))
    return total-tail, total+tail


def decimal(v, digits=30, upward=False):
    """Round an exact rational outward (down, or up if upward) to a decimal string for the record."""
    scaled = v*10**digits
    k = -((-scaled.numerator)//scaled.denominator) if upward else scaled.numerator//scaled.denominator
    sign, k = ('-' if k < 0 else ''), abs(k)
    return f"{sign}{k//10**digits}.{k % 10**digits:0{digits}d}"


def h_bounds(g, pi):
    """Enclosure of h(g) = (2/pi) sin(pi (1/3+g)) - 1/3 - 2g, for 0 <= g <= 1/6.

    The angle pi(1/3+g) lies in [pi/3, pi/2], where sine increases, so the lower
    and upper ends of the pi interval give the two sine bounds."""
    lo_pi, hi_pi = pi
    s_lo = sin_interval(lo_pi*(F(1, 3)+g))[0]
    s_hi = sin_interval(hi_pi*(F(1, 3)+g))[1]
    assert hi_pi*(F(1, 3)+g) < 2 and lo_pi*(F(1, 3)+g) >= 0
    return 2*s_lo/hi_pi-F(1, 3)-2*g, 2*s_hi/lo_pi-F(1, 3)-2*g


def threshold_certificate():
    pi = pi_interval()
    assert pi[1]-pi[0] < F(1, 10**80) and F(314159, 10**5) < pi[0] < pi[1] < F(314160, 10**5)
    lo, hi = G0_BRACKET
    at_lo, at_hi = h_bounds(lo, pi), h_bounds(hi, pi)
    at_zero, at_sixth = h_bounds(F(0), pi), h_bounds(F(1, 6), pi)
    assert at_lo[0] > 0 > at_hi[1] and at_zero[0] > 0 > at_sixth[1]
    # h'(g) = 2cos(pi(1/3+g)) - 2 <= -1 on [0, 1/6], since the angle is in [pi/3, pi/2]: h is strictly decreasing.
    simple = (2/pi[1]-F(1, 2), 2/pi[0]-F(1, 2))     # 2/pi - 1/2, a simpler closed-form threshold below g0
    psi_third = (2*sin_interval(pi[0]/3)[0]/pi[1]-F(1, 3), 2*sin_interval(pi[1]/3)[1]/pi[0]-F(1, 3))
    clean_gap = F(90, 1000)
    step_two = F(1, 4)+pi[1]**2/3*F(1, 64)
    assert simple[1] < lo and clean_gap < lo and psi_third[1] < F(1, 4) and step_two < F(1, 3)
    return dict(status='passed', function='h(g) = (2/pi) sin(pi(1/3+g)) - 1/3 - 2g on [0, 1/6]',
        root_bracket=[str(lo), str(hi)], root_bracket_decimal=['0.15127097021660', '0.15127097021670'],
        h_at_lower_end_lower_bound=decimal(at_lo[0], 40), h_at_upper_end_upper_bound=decimal(at_hi[1], 40, upward=True),
        h_at_zero_lower_bound=decimal(at_zero[0]), h_at_one_sixth_upper_bound=decimal(at_sixth[1], upward=True),
        strictly_decreasing='h\'(g) = 2cos(pi(1/3+g)) - 2 <= -1 on [0,1/6]', unique_root=True,
        pi_enclosure_width_below='1e-80', pi_method='Machin: 16 arctan(1/5) - 4 arctan(1/239), alternating tails',
        sine_method='Taylor sum with alternating tail bound',
        simpler_threshold_2_over_pi_minus_half_upper=decimal(simple[1], upward=True), margin_over_simpler_threshold_lower=decimal(lo-simple[1]),
        clean_circle_largest_gap_bound='0.090', margin_over_clean_circle_gap_lower=decimal(lo-clean_gap),
        psi_one_third_upper=decimal(psi_third[1], upward=True), psi_one_third_below_one_quarter=True,
        step_two_bound='1/4 + (pi^2/3)(1/4)^3 < 1/3', step_two_value_upper=decimal(step_two, upward=True))


def build():
    return dict(six_point_example=six_point_certificate(), threshold_root=threshold_certificate())


def verify(root, require):
    root = Path(root)
    folder = root/FOLDER
    saved = json.loads((folder/'certificates.json').read_text())
    require('rational_certificates_script_hash', saved['script_sha256'] == hashlib.sha256((root/'code/rational_certificates.py').read_bytes()).hexdigest())
    again = build()
    require('rational_certificates_six_point', again['six_point_example'] == saved['six_point_example'] and saved['six_point_example']['status'] == 'passed')
    require('rational_certificates_threshold', again['threshold_root'] == saved['threshold_root'] and saved['threshold_root']['status'] == 'passed')
    return dict(status='passed', six_point_bar=[saved['six_point_example']['bar']['birth'], saved['six_point_example']['bar']['death']],
                g0_bracket=saved['threshold_root']['root_bracket_decimal'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.verify:
        def require(name, condition, detail=None):
            if not bool(condition): raise AssertionError(name)
        print(json.dumps(verify(ROOT, require), indent=2))
    else:
        report = build()
        report['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        (ROOT/FOLDER).mkdir(parents=True, exist_ok=True)
        (ROOT/FOLDER/'certificates.json').write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(report, indent=2))
