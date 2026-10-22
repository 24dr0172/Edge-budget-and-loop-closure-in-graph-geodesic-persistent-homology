"""Recompute the seven-point example of Appendix C independently (no release helpers imported).

Convention (as in the manuscript's chord formula and in code/finite_diagram.py and
code/generate_witness.py): points x_i = (cos 2 pi s_i, sin 2 pi s_i)/(2 pi) on the
circle of circumference 1, so c(i,j) = chi(d_star) with chi(r) = sin(pi r)/pi.
Vertices are numbered x_0..x_6 here; the manuscript numbers them x_1..x_7.
Checks: the intrinsic bar by our own F2 reduction and by Ripser; radius order and
closure; graph metrics at M=13 and 14; bars by our own F2 reduction, by Ripser and
by explicit threshold-graph and triangle facts; both bottleneck matchings; and arc
support of every edge (both arcs, both endpoints).
Output: results/seven_point_check/seven_point_numbers.json
"""
from pathlib import Path
import itertools, json, math
from fractions import Fraction
import numpy as np
from scipy.sparse.csgraph import shortest_path
from ripser import ripser

S = ['0', '0.08', '0.15', '0.161', '0.29', '0.60', '0.66']
s = np.array(S, float)
n = len(s)
chi = lambda r: math.sin(math.pi * r) / math.pi
Q = [Fraction(x) for x in S]
dq = {(i, j): min(abs(Q[i] - Q[j]), 1 - abs(Q[i] - Q[j])) for i in range(n) for j in range(n)}
dstar = np.array([[float(dq[i, j]) for j in range(n)] for i in range(n)])
xy = np.c_[np.cos(2 * np.pi * s), np.sin(2 * np.pi * s)] / (2 * np.pi)
euc = np.linalg.norm(xy[:, None] - xy[None], axis=2)
chord = np.vectorize(chi)(dstar)
out = dict(coordinates=S, convention='x_i=(cos 2pi s_i, sin 2pi s_i)/(2pi); circumference 1; c=chi(d_star), chi(r)=sin(pi r)/pi; vertices 0..6 in curve order')
checks = {}
def check(name, ok, detail=None):
    checks[name] = dict(passed=bool(ok), detail=detail)
check('embedding_chord_equals_chi_of_intrinsic', np.abs(euc - chord).max() < 1e-15, float(np.abs(euc - chord).max()))
pairs = list(itertools.combinations(range(n), 2))
check('intrinsic_pair_distances_distinct_exact', len({dq[p] for p in pairs}) == len(pairs))

def h1_bars(D, tol=1e-12):
    """Plain F2 boundary-matrix reduction of the Rips filtration (dims <= 2)."""
    simp = []
    for k in (1, 2, 3):
        for v in itertools.combinations(range(len(D)), k):
            simp.append((max([D[a, b] for a, b in itertools.combinations(v, 2)], default=0.0), k, v))
    simp.sort()
    idx = {v: t for t, (_, _, v) in enumerate(simp)}
    low = {}
    bars = []
    for t, (val, k, v) in enumerate(simp):
        col = set(idx[v[:j] + v[j + 1:]] for j in range(k)) if k > 1 else set()
        while col and max(col) in low:
            col ^= low[max(col)]
        if col:
            p = max(col); low[p] = col
            if simp[p][1] == 2 and val - simp[p][0] > tol:
                bars.append((simp[p][0], val))
    return sorted(bars)

def ripser_h1(D):
    dg = ripser(D, distance_matrix=True, maxdim=1)['dgms'][1]
    return sorted((float(a), float(b)) for a, b in dg if b > a)

# Reference
gaps = [float((Q[(i + 1) % n] - Q[i]) % 1) for i in range(n)]
g = max(gaps)
ref = h1_bars(dstar)
refR = ripser_h1(dstar)
diam = float(dstar.max())
check('g_is_0.34', abs(g - 0.34) < 1e-15 and gaps.index(g) == 6, dict(gaps=gaps))
check('intrinsic_bar_[0.34,0.37)', len(ref) == 1 and abs(ref[0][0] - 0.34) < 1e-12 and abs(ref[0][1] - 0.37) < 1e-12, ref)
check('ripser_intrinsic_bar', len(refR) == 1 and np.allclose(refR[0], ref[0], atol=1e-6), refR)
check('intrinsic_diameter_0.499', abs(diam - 0.499) < 1e-15 and dq[3, 6] == Fraction('0.499'))
b_star = ref[0][1]
# Radius order: chord order (all distinct)
order = sorted(pairs, key=lambda p: (chord[p], p))
ring = {(i, (i + 1) % n) if i < (i + 1) % n else ((i + 1) % n, i) for i in range(n)}
pos = {p: k + 1 for k, p in enumerate(order)}
Mclose = max(pos[e] for e in ring)
check('radius_closure_M13', Mclose == 13 and order[12] == (0, 6) and order[13] == (4, 6),
      dict(order=[list(p) for p in order[:14]]))

def graph_metric(M):
    W = np.zeros((n, n))
    for (a, b) in order[:M]:
        W[a, b] = W[b, a] = chord[a, b]
    D, pred = shortest_path(W, directed=False, return_predecessors=True)
    return D, pred

def path(pred, i, j):
    p = [j]
    while p[-1] != i:
        p.append(int(pred[i, p[-1]]))
    return p[::-1]

res = {}
for M in (12, 13, 14):
    D, pred = graph_metric(M)
    if not np.isfinite(D).all():
        res[M] = dict(connected=False); continue
    bars = h1_bars(D); barsR = ripser_h1(D)
    paths = {f'{i}-{j}': path(pred, i, j) for i, j in pairs}
    res[M] = dict(connected=True, bars_F2=bars, bars_ripser=barsR, D=D, paths=paths)
# M=12: connected (all 7 vertices), only the ring edge 0-6 (manuscript x7x1) still missing, so no loop: empty H1, error P_star.
missing12 = sorted(e for e in ring if pos[e] > 12)
check('M12_connected_no_loop_only_ring_edge_0-6_missing', res[12]['connected'] and res[12]['bars_F2'] == [] and len(res[12]['bars_ripser']) == 0 and missing12 == [(0, 6)],
      dict(connected=res[12]['connected'], bars_F2=res[12].get('bars_F2'), missing_ring_edges=[list(e) for e in missing12]))
c34, c31, c06, c37 = chi(0.34), chi(0.31), chi(0.06), chi(0.37)
d13 = res[13]['bars_F2']; d14 = res[14]['bars_F2']
check('M13_bar_[chi(.34),chi(.31)+chi(.06))', len(d13) == 1 and abs(d13[0][0] - c34) < 1e-15 and abs(d13[0][1] - (c31 + c06)) < 1e-15, d13)
check('M14_bar_[chi(.34),chi(.37))', len(d14) == 1 and abs(d14[0][0] - c34) < 1e-15 and abs(d14[0][1] - c37) < 1e-15, d14)
check('ripser_agrees_M13_M14', all(len(res[M]['bars_ripser']) == 1 and np.allclose(res[M]['bars_ripser'][0], res[M]['bars_F2'][0], atol=1e-6) for M in (13, 14)),
      {M: res[M]['bars_ripser'] for M in (13, 14)})
check('D13(4,6)_path_4-5-6', res[13]['paths']['4-6'] == [4, 5, 6] and abs(res[13]['D'][4, 6] - (c31 + c06)) < 1e-15)
check('D14(4,6)_direct', res[14]['paths']['4-6'] == [4, 6] and abs(res[14]['D'][4, 6] - c37) < 1e-15)

# Hand certificate facts (threshold graphs of the graph metrics)
def thr_edges(D, t):
    return sorted([i, j] for i, j in pairs if D[i, j] <= t + 1e-15)
def tri(D, t):
    return [list(v) for v in itertools.combinations(range(n), 3) if all(D[a, b] <= t + 1e-15 for a, b in itertools.combinations(v, 2))]
cert = {}
for M in (13, 14):
    D = res[M]['D']; b, d = res[M]['bars_F2'][0]
    vals = sorted(set(np.round(D[np.triu_indices(n, 1)], 15)))
    below_birth = max(v for v in vals if v < b - 1e-15)
    below_death = max(v for v in vals if v < d - 1e-15)
    absent_at_birth = [[i, j] for i, j in pairs if D[i, j] > b + 1e-15]
    absent_below_death = [[i, j] for i, j in pairs if D[i, j] > below_death + 1e-15]
    cert[M] = dict(birth=b, death=d, largest_value_below_birth=below_birth, largest_value_below_death=below_death,
                   pairs_absent_at_birth=absent_at_birth, pairs_absent_just_below_death=absent_below_death,
                   pair_entering_at_death=[[i, j] for i, j in pairs if abs(D[i, j] - d) < 1e-15],
                   triangles_containing_0_6_at_death=[v for v in tri(D, d) if 0 in v and 6 in v],
                   distances={f'{i}-{j}': dict(value=float(D[i, j]), path=res[M]['paths'][f'{i}-{j}']) for i, j in pairs})
    # H1 rank checks just below birth/death and at death via own reduction on sub-threshold complexes
# Hand certificate: on [birth, death) the threshold graph is the clique on {0,..,4} plus the path 4-5-6-0;
# below birth the edge 0-6 is missing (contractible: clique plus pendant path 4-5-6);
# at death the pair 4-6 enters, filling triangles 4-5-6 and 0-4-6 (all of {0,..,6} then cones off).
hand = {}
for M in (13, 14):
    D = res[M]['D']; b, d = res[M]['bars_F2'][0]
    t = (b + d) / 2
    nb = {v: sorted(j for j in range(n) if j != v and D[v, j] <= t) for v in range(n)}
    clique = all(D[i, j] <= t for i, j in itertools.combinations(range(5), 2))
    only_values = sorted({float(D[i, j]) for i, j in pairs if b - 1e-15 < D[i, j] < d - 1e-15})
    hand[M] = dict(neighbours_of_5=nb[5], neighbours_of_6=nb[6], clique_0_to_4=clique,
                   birth_is_chord_0_6=abs(b - chord[0, 6]) < 1e-15 and D[0, 6] == chord[0, 6],
                   no_pair_value_strictly_inside_bar=(only_values == [] or only_values == [b]),
                   death_is_D_4_6=abs(d - D[4, 6]) < 1e-15,
                   triangles_456_and_046_present_at_death=all(D[i, j] <= d + 1e-15 for i, j in [(4, 5), (5, 6), (4, 6), (0, 4), (0, 6)]),
                   vertex_6_dominated_by_4_at_death=all(D[4, j] <= d + 1e-15 for j in range(n) if D[6, j] <= d + 1e-15))
    ok = nb[5] == [4, 6] and nb[6] == [0, 5] and all(v for k, v in hand[M].items() if not k.startswith('neigh'))
    check(f'hand_certificate_M{M}', ok, hand[M])
out['hand_certificate'] = hand
out['certificate_facts'] = cert

# Matchings
def costs(bar, refbar, diam):
    direct = max(abs(bar[0] - refbar[0]), abs(bar[1] - refbar[1]))
    diag = max((bar[1] - bar[0]) / 2, (refbar[1] - refbar[0]) / 2)
    return dict(direct=direct, direct_normalized=direct / diam, both_to_diagonal=diag, both_to_diagonal_normalized=diag / diam,
                bottleneck_normalized=min(direct, diag) / diam)
def bottleneck_persim(a, b):
    import persim
    return float(persim.bottleneck(np.array(a), np.array(b)))
P_star = (b_star - g) / (2 * diam)
m13 = costs(d13[0], ref[0], diam); m14 = costs(d14[0], ref[0], diam)
m13['persim_bottleneck_normalized'] = bottleneck_persim(d13, ref) / diam
m14['persim_bottleneck_normalized'] = bottleneck_persim(d14, ref) / diam
check('error_M13_0.0440643', abs(m13['bottleneck_normalized'] - 0.0440643) < 5e-8 and abs(m13['persim_bottleneck_normalized'] - m13['bottleneck_normalized']) < 1e-12, m13['bottleneck_normalized'])
check('error_M14_0.0300601_eq_P_star', abs(m14['bottleneck_normalized'] - P_star) < 1e-15 and abs(P_star - 0.0300601) < 5e-8 and abs(m14['persim_bottleneck_normalized'] - P_star) < 1e-12, m14['bottleneck_normalized'])
check('direct_cost_0.122371_to_0.156051', abs(m13['direct_normalized'] - 0.122371) < 5e-7 and abs(m14['direct_normalized'] - 0.156051) < 5e-7,
      [m13['direct_normalized'], m14['direct_normalized']])
check('M13_diagonal_matching_optimal_and_M14_reference_half_lifetime_dominates',
      m13['both_to_diagonal'] < m13['direct'] and m14['both_to_diagonal'] == (b_star - g) / 2 and (d14[0][1] - d14[0][0]) / 2 < (b_star - g) / 2)
release = json.load(open(Path(__import__('os').environ.get('RELEASE_ROOT', Path(__file__).resolve().parent.parent)) / 'results/useful_diagram/rational_counterexample.json'))
rr = release['methods']['R2_global']
check('matches_release_json', abs(rr['before']['e_PH'] - m13['bottleneck_normalized']) < 1e-12 and abs(rr['after']['e_PH'] - m14['bottleneck_normalized']) < 1e-12
      and abs(rr['before']['death'] - d13[0][1]) < 1e-12 and abs(rr['after']['death'] - d14[0][1]) < 1e-12)

# Arc support (Definition arc-support): for each edge, some arc and endpoint p with
# consecutive chords on arc <= c(a,b) and D_G(p,w) <= c(a,b) for all w on the arc.
def arcs(a, b):
    fwd = [a]; k = a
    while k != b:
        k = (k + 1) % n; fwd.append(k)
    bwd = [a]; k = a
    while k != b:
        k = (k - 1) % n; bwd.append(k)
    return [fwd, bwd]
def arc_support(M, D):
    rep = {}
    for (a, b) in order[:M]:
        c = chord[a, b]; ok = []
        for arc in arcs(a, b):
            steps = all(chord[arc[i - 1], arc[i]] <= c + 1e-15 for i in range(1, len(arc)))
            for p in (a, b):
                if steps and all(D[p, w] <= c + 1e-15 for w in arc):
                    ok.append(dict(arc=arc, p=p))
        rep[f'{a}-{b}'] = ok[0] if ok else None
    return rep
for M in (13, 14):
    rep = arc_support(M, res[M]['D'])
    check(f'arc_support_M{M}', all(v is not None for v in rep.values()), rep)
# also every prefix M >= 13 up to complete
allok = True
for M in range(13, len(pairs) + 1):
    D, _ = graph_metric(M)
    allok &= all(v is not None for v in arc_support(M, D).values())
check('arc_support_all_prefixes_M13_to_21', allok)

# Small-gap corollary hypothesis and cost ceiling (Cor. cost-ceiling)
check('small_gap_hypothesis_fails', g >= 2 / math.pi - 0.5 and g >= 1 / 3, dict(g=g, threshold=2 / math.pi - 0.5))
ceiling = max(g - chi(g), b_star - chi(b_star)) / diam
check('death_bracket_chi(b*)<=d_G<=b*', all(chi(b_star) - 1e-15 <= res[M]['bars_F2'][0][1] <= b_star for M in (13, 14)))
check('errors_below_finite_ceiling', m13['bottleneck_normalized'] <= ceiling and m14['bottleneck_normalized'] <= ceiling, ceiling)

# kNN orders (secondary): rank keys with label tie-break, as in Definition rank-key
rank = np.zeros((n, n), int)
for i in range(n):
    others = sorted([j for j in range(n) if j != i], key=lambda j: (chord[i, j], j))
    for r, j in enumerate(others):
        rank[i, j] = r + 1
knn = {}
for name, f in [('OR_kNN', min), ('mutual_kNN', max)]:
    o = sorted(pairs, key=lambda p: (f(rank[p[0], p[1]], rank[p[1], p[0]]), chord[p], p))
    mc = max(o.index(e) + 1 for e in ring)
    j = o.index((4, 6)) + 1
    def gm(M):
        W = np.zeros((n, n))
        for (a, b) in o[:M]:
            W[a, b] = W[b, a] = chord[a, b]
        return shortest_path(W, directed=False)
    knn[name] = dict(M_close=mc, insertion_of_4_6=j, bar_before=h1_bars(gm(j - 1)), bar_after=h1_bars(gm(j)))
check('kNN_release_transitions', knn['OR_kNN']['M_close'] == 10 and knn['OR_kNN']['insertion_of_4_6'] == 13 and
      knn['mutual_kNN']['M_close'] == 13 and knn['mutual_kNN']['insertion_of_4_6'] == 17 and
      all(np.allclose(knn[k]['bar_before'], d13) and np.allclose(knn[k]['bar_after'], d14) for k in knn), knn)

out.update(g=g, gaps=gaps, reference_bar=ref[0], intrinsic_diameter=diam, P_star=P_star,
           radius_order=[list(p) for p in order], M_close_radius=Mclose,
           chi=dict(chi_034=c34, chi_031=c31, chi_006=c06, chi_037=c37, chi_031_plus_chi_006=c31 + c06),
           bars={'M13': d13[0], 'M14': d14[0]}, ripser_bars={'M13': res[13]['bars_ripser'], 'M14': res[14]['bars_ripser']},
           matchings={'M13': m13, 'M14': m14}, finite_ceiling_normalized=ceiling,
           small_gap_threshold=2 / math.pi - 0.5, knn=knn, checks=checks)
OUT = Path(__file__).resolve().parent.parent / 'results' / 'seven_point_check'
OUT.mkdir(parents=True, exist_ok=True)
with open(OUT / 'seven_point_numbers.json', 'w') as stream:
    json.dump(out, stream, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else float(x))
for k, v in checks.items():
    print(('PASS ' if v['passed'] else 'FAIL ') + k)
print(json.dumps(dict(bars=out['bars'], matchings=out['matchings'], P_star=P_star, ceiling=ceiling), indent=1))
print(json.dumps({M: {k: cert[M][k] for k in cert[M] if k != 'distances'} for M in cert}, indent=1))
