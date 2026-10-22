#!/usr/bin/env python3
"""Independent check of the numbers in Appendix C, Section 7.2, Figure 2 and Figure 3.

The script uses only numpy, scipy and its own code (Ripser is an optional
cross-check). Conventions follow the manuscript: chord-weighted graphs, radius order sorts
pairs by chord length (labels break ties), exact-M prefixes, Vietoris--Rips
flag complexes of the graph shortest-path metric D_G, bars [b,d) with d>b,
coefficients F2 (Ripser) plus an independent F2 boundary-matrix reduction
written here, bottleneck endpoint displacement = max(|db|,|dd|), and
normalised error e_PH = d_B / Delta_star.

Run from the release folder:
    python3 code/verify_independent_checks.py
The saved output is results/independent_checks/independent_checks_output.txt;
code/verify_retention_cost_witness.py reruns this script and compares the two.
"""
import itertools, math, sys
from fractions import Fraction
from decimal import Decimal, getcontext
import numpy as np

FAIL = []
def check(cond, msg):
    print(("  [ok]   " if cond else "  [FAIL] ") + msg)
    if not cond: FAIL.append(msg)

def chi(r):  # chord of an arc of length r on the unit-circumference circle
    return math.sin(math.pi * r) / math.pi

# ---------------------------------------------------------------- utilities
def dstar(s, a, b):
    d = abs(s[a] - s[b]); return min(d, 1 - d)

def shortest_paths(n, edges, cost):
    D = np.full((n, n), np.inf); np.fill_diagonal(D, 0.0)
    for (i, j) in edges:
        D[i, j] = D[j, i] = min(D[i, j], cost[i, j])
    for k in range(n):                      # Floyd--Warshall (n is tiny)
        D = np.minimum(D, D[:, [k]] + D[[k], :])
    return D

def f2_reduction_h1(D):
    """Own F2 column reduction of the flag filtration up to dimension 2.
    Simplices sorted by (filtration value, dimension, vertex tuple)."""
    n = len(D); simp = []
    for size in (1, 2, 3):
        for v in itertools.combinations(range(n), size):
            val = max((D[a, b] for a, b in itertools.combinations(v, 2)), default=0.0)
            simp.append((val, size, v))
    simp.sort()
    idx = {v: k for k, (_, _, v) in enumerate(simp)}
    low = {}; bars = []
    for k, (val, size, v) in enumerate(simp):
        col = 0
        if size > 1:
            for j in range(size):
                col ^= 1 << idx[v[:j] + v[j + 1:]]
        while col:
            p = col.bit_length() - 1
            if p not in low: break
            col ^= low[p]
        if col:
            p = col.bit_length() - 1; low[p] = col
            if simp[p][1] == 2 and val > simp[p][0] + 1e-12:
                bars.append((simp[p][0], val))
    # essential H1 classes cannot occur (complete complex at max value)
    return sorted(bars)

def ripser_h1(D):
    try:
        from ripser import ripser
    except Exception:
        return None
    dg = ripser(np.asarray(D, float), distance_matrix=True, maxdim=1)['dgms'][1]
    return sorted((float(b), float(d)) for b, d in dg if d > b)

def rank_mod(rows, p=None):
    """Rank of an integer matrix: over F2 if p==2, else over Q (Fractions)."""
    M = [[Fraction(x) for x in r] for r in rows]
    if p == 2: M = [[x % 2 for x in r] for r in rows]
    rank = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((r for r in range(rank, len(M)) if M[r][c] != 0), None)
        if piv is None: continue
        M[rank], M[piv] = M[piv], M[rank]
        for r in range(len(M)):
            if r != rank and M[r][c] != 0:
                if p == 2:
                    M[r] = [(a + b) % 2 for a, b in zip(M[r], M[rank])]
                else:
                    f = M[r][c] / M[rank][c]
                    M[r] = [a - f * b for a, b in zip(M[r], M[rank])]
        rank += 1
    return rank

def betti1_at(D, t):
    """beta_1 of VR(D;t) (flag complex) over F2 and over Q."""
    n = len(D)
    E = [(i, j) for i, j in itertools.combinations(range(n), 2) if D[i, j] <= t]
    T = [v for v in itertools.combinations(range(n), 3)
         if all(D[a, b] <= t for a, b in itertools.combinations(v, 2))]
    # components
    parent = list(range(n))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, j in E: parent[find(i)] = find(j)
    comps = len({find(i) for i in range(n)})
    eidx = {e: k for k, e in enumerate(E)}
    rows = []
    for (a, b, c) in T:          # boundary [b,c]-[a,c]+[a,b]
        r = [0] * len(E); r[eidx[(b, c)]] += 1; r[eidx[(a, c)]] -= 1; r[eidx[(a, b)]] += 1
        rows.append(r)
    z1 = len(E) - n + comps
    r2 = rank_mod(rows, 2) if rows else 0
    rq = rank_mod(rows) if rows else 0
    return E, T, z1 - r2, z1 - rq

def matching_costs(bar, ref):
    (b, d), (bs, ds) = bar, ref
    direct = max(abs(b - bs), abs(d - ds))
    diag = max((d - b) / 2, (ds - bs) / 2)
    return direct, diag, min(direct, diag)

def arc_support(n, edges, c, D, a, b, cyclic=True):
    """Return list of (arc, p, ok, worst_margin) for edge {a,b} (Definition
    def:arc-support): consecutive chords on the arc <= c(a,b) and
    D_G(p,w) <= c(a,b) for every sample w on that arc."""
    out = []
    arcs = []
    fw = [a]; k = a
    while k != b: k = (k + 1) % n; fw.append(k)
    bw = [a]; k = a
    while k != b: k = (k - 1) % n; bw.append(k)
    for arc in (fw, bw):
        for p in (a, b):
            h = c[a, b]
            m1 = min(h - c[arc[i], arc[i + 1]] for i in range(len(arc) - 1))
            m2 = min(h - D[p, w] for w in arc)
            out.append((arc, p, m1 >= -1e-12 and m2 >= -1e-12, min(m1, m2)))
    return out

def edge_supported(n, edges, c, D, e):
    return any(ok for _, _, ok, _ in arc_support(n, edges, c, D, *e))

# ================================================================ Part 1
print("=" * 78); print("PART 1. Seven-point circle example (Appendix C)"); print("=" * 78)
s = [0, 0.08, 0.15, 0.161, 0.29, 0.60, 0.66]
sF = [Fraction(x).limit_denominator(1000) for x in ["0", "0.08", "0.15", "0.161", "0.29", "0.60", "0.66"]] if False else \
     [Fraction(0), Fraction(8, 100), Fraction(15, 100), Fraction(161, 1000), Fraction(29, 100), Fraction(60, 100), Fraction(66, 100)]
n = 7
gapsF = [(sF[(i + 1) % n] - sF[i]) % 1 for i in range(n)]
print("arc coordinates:", s)
print("consecutive gaps x_i x_{i+1} (exact):", [str(x) for x in gapsF])
gF = max(gapsF)
check(gF == Fraction(34, 100), f"g = {float(gF)} (gap x7 x1)")
dF = {}
for i, j in itertools.combinations(range(n), 2):
    d = abs(sF[i] - sF[j]); dF[(i, j)] = min(d, 1 - d)
vals = sorted(dF.values())
check(len(set(vals)) == 21, "all 21 intrinsic distances distinct => all chords distinct (chi strictly increasing)")
DeltaF = max(vals)
check(DeltaF == Fraction(499, 1000), f"intrinsic diameter Delta* = {float(DeltaF)} (pair x4 x7)")
order = sorted(dF, key=lambda p: (dF[p], p))
print("radius order (position: pair, intrinsic distance, chord):")
cons = {(i, (i + 1) % n) if i < (i + 1) % n else ((i + 1) % n, i) for i in range(n)}
for k, p in enumerate(order):
    print(f"   {k+1:2d}: x{p[0]+1}x{p[1]+1}  d*={float(dF[p]):.3f}  c={chi(float(dF[p])):.6f}" + ("  ring" if p in cons else ""))
pos = {p: k + 1 for k, p in enumerate(order)}
mu = sorted(pos[e] for e in cons)
check(mu[-1] == 13 and mu[-2] == 12, f"ring insertion positions {mu}: M_conn = mu_(n-1) = 12, M_close = 13")
check(order[13] == (4, 6), "the 14th inserted pair is x5x7 (d* = 0.37)")
check(all(dF[p] <= Fraction(29, 100) for p in itertools.combinations(range(5), 2)),
      "x1..x5 are pairwise within intrinsic 0.29 (they lie on [0,0.29])")
check(set(order[:13]) == set(itertools.combinations(range(5), 2)) | {(4, 5), (5, 6), (0, 6)},
      "G_13 = complete graph on x1..x5 plus path x5-x6-x7-x1")

# intrinsic diagram
Dst = np.zeros((n, n))
for (i, j), v in dF.items(): Dst[i, j] = Dst[j, i] = float(v)
bars_int = f2_reduction_h1(Dst); rip_int = ripser_h1(Dst)
print("intrinsic H1 bars, own F2 reduction:", bars_int, " ripser:", rip_int)
check(len(bars_int) == 1 and abs(bars_int[0][0] - 0.34) < 1e-12 and abs(bars_int[0][1] - 0.37) < 1e-12,
      "intrinsic reference bar [g,b*) = [0.34, 0.37)")
if rip_int is not None:
    check(len(rip_int) == 1 and abs(rip_int[0][0] - .34) < 1e-6 and abs(rip_int[0][1] - .37) < 1e-6, "ripser agrees")
gs, bs, Dl = 0.34, 0.37, 0.499
# certificate pieces for intrinsic bar
check(min(v for p, v in dF.items() if p not in set(order[:13])) == Fraction(37, 100),
      "smallest intrinsic distance outside G_13 pairs is 0.37 (x5x7); next 0.40 (x1x6)")
check(sorted(dF.values())[14] == Fraction(40, 100), "15th smallest intrinsic distance is 0.40")

c = np.zeros((n, n))
for (i, j), v in dF.items(): c[i, j] = c[j, i] = chi(float(v))
B = chi(0.34); h = chi(0.06); a31 = chi(0.31); b29 = chi(0.29)
delta = chi(0.31) + chi(0.06)
print(f"chi(0.34) = {B:.10f}\nchi(0.37) = {chi(0.37):.10f}\nchi(0.31) = {a31:.10f}\nchi(0.06) = {h:.10f}")
print(f"chi(0.31)+chi(0.06) = {delta:.10f}\nchi(0.29) = {b29:.10f}\nchi(0.08) = {chi(0.08):.10f}\nchi(0.129) = {chi(0.129):.10f}")
print(f"chi(0.06)+chi(0.34) = {h+B:.10f}\nchi(0.34)+chi(0.08) = {B+chi(0.08):.10f}\nchi(0.31)+chi(0.129) = {a31+chi(0.129):.10f}")
print(f"chi(0.34)+chi(0.29) = {B+b29:.10f}\nchi(0.40) = {chi(0.40):.10f}   chi(0.499) = {chi(0.499):.10f}")

results_C = {}
for M in range(12, 22):
    E = order[:M]
    D = shortest_paths(n, E, c)
    connected = np.isfinite(D).all()
    bars = f2_reduction_h1(D) if connected else None
    rb = ripser_h1(D) if connected else None
    results_C[M] = (E, D, bars, rb)
    if connected:
        if bars:
            dir_, dia, eb = matching_costs(bars[0], (gs, bs))
        else:
            dir_, dia, eb = float('nan'), float('nan'), (bs - gs) / 2
        print(f"  M={M:2d}: own-F2 bars={[(round(x,7),round(y,7)) for x,y in bars]}  ripser={[(round(x,6),round(y,6)) for x,y in rb] if rb is not None else 'n/a'}  e_PH={eb/Dl:.7f}")
E13, D13, bars13, rb13 = results_C[13]
E14, D14, bars14, rb14 = results_C[14]
check(results_C[12][2] == [], "M=12 (first connected prefix, before closure): empty H1 diagram, error P*")
check(len(bars13) == 1 and abs(bars13[0][0] - B) < 1e-12 and abs(bars13[0][1] - delta) < 1e-12,
      f"M=13 bar = [chi(0.34), chi(0.31)+chi(0.06)) = [{B:.7f}, {delta:.7f})")
check(len(bars14) == 1 and abs(bars14[0][0] - B) < 1e-12 and abs(bars14[0][1] - chi(0.37)) < 1e-12,
      f"M=14 bar = [chi(0.34), chi(0.37)) = [{B:.7f}, {chi(0.37):.7f})")
for M in range(14, 22):
    bm = results_C[M][2]
    check(len(bm) == 1 and abs(bm[0][1] - chi(0.37)) < 1e-12 and abs(bm[0][0] - B) < 1e-12,
          f"M={M}: single bar with death chi(0.37) (squeezed by Thm circle-parameters(d))")
for M in (13, 14):
    rb = results_C[M][3]
    if rb is not None:
        ref = results_C[M][2]
        check(len(rb) == 1 and max(abs(rb[0][0] - ref[0][0]), abs(rb[0][1] - ref[0][1])) < 1e-6,
              f"ripser (float32) agrees with own F2 reduction at M={M} within 1e-6")

# certificate: graph distances of absent pairs at M=13
absent13 = [p for p in itertools.combinations(range(n), 2) if p not in set(E13)]
print("M=13 graph distances of the 8 absent pairs (sorted):")
for p in sorted(absent13, key=lambda p: D13[p]):
    print(f"   x{p[0]+1}x{p[1]+1}: D_G = {D13[p]:.6f}   chord = {c[p]:.6f}")
sa = sorted(absent13, key=lambda p: D13[p])
check(sa[0] == (4, 6) and abs(D13[4, 6] - delta) < 1e-14, "smallest absent-pair D_G at M=13 is D_G(x5,x7)=chi(0.31)+chi(0.06), uniquely")
check(abs(D13[sa[1]] - (h + B)) < 1e-14 and sa[1] == (0, 5), f"second smallest is D_G(x1,x6)=chi(0.06)+chi(0.34)={h+B:.6f} > delta")
check(all(abs(D13[p] - c[p]) < 1e-15 for p in E13), "every edge of G_13 has D_G equal to its chord")
check(max(c[p] for p in E13) == B, "largest chord in G_13 is B = chi(0.34) < delta")
# formulas used in the certificate
ok = True
for w in range(5):
    f6 = min(c[5, 4] + c[4, w], c[5, 6] + c[6, 0] + c[0, w]); f7 = min(c[6, 0] + c[0, w], c[6, 5] + c[5, 4] + c[4, w])
    ok &= abs(D13[5, w] - f6) < 1e-14 and abs(D13[6, w] - f7) < 1e-14
check(ok, "D_G(x6,w)=min{c56+c5w, c67+c71+c1w}, D_G(x7,w)=min{c71+c1w, c76+c65+c5w} for w in x1..x5")
check(min(c[0, w] for w in range(1, 5)) == chi(0.08) and min(c[4, w] for w in range(4)) == chi(0.129),
      "nearest cluster chords: from x1 chi(0.08), from x5 chi(0.129)")
# threshold complex on [B, delta) and at delta
for t, lab in ((B, "t=B"), ((B + delta) / 2, "t in (B,delta)"), (delta - 1e-9, "t=delta-")):
    Ee, Tt, b2, bq = betti1_at(D13, t)
    check(set(Ee) == set(E13) and b2 == 1 and bq == 1, f"M=13, {lab}: threshold graph = G_13, beta1 = 1 over F2 and Q")
Ee, Tt, b2, bq = betti1_at(D13, delta)
check(set(Ee) == set(E13) | {(4, 6)} and b2 == 0 and bq == 0 and (4, 5, 6) in Tt and (0, 4, 6) in Tt,
      "M=13, t=delta: only x5x7 added; triangles x5x6x7 and x1x5x7 present; beta1 = 0 over F2 and Q")
Ee, Tt, b2, bq = betti1_at(D13, B - 1e-9)
check(b2 == 0 and bq == 0, "M=13, t<B: beta1 = 0")
# M=14
absent14 = [p for p in itertools.combinations(range(n), 2) if p not in set(E14)]
check(min(D14[p] for p in absent14) > chi(0.37), f"M=14: every absent pair has D_G > chi(0.37) (min {min(D14[p] for p in absent14):.6f})")
check(all(c[p] >= chi(0.37) - 1e-15 for p in absent13), "every pair absent from G_13 has chord >= chi(0.37) (so D_G >= chi(0.37))")
for t, lab in ((B, "t=B"), (chi(0.37) - 1e-9, "t=chi(0.37)-")):
    Ee, Tt, b2, bq = betti1_at(D14, t)
    check(set(Ee) == set(E13) and b2 == 1 and bq == 1, f"M=14, {lab}: threshold graph = G_13, beta1 = 1")
Ee, Tt, b2, bq = betti1_at(D14, chi(0.37))
check(set(Ee) == set(E14) and b2 == 0 and bq == 0, "M=14, t=chi(0.37): graph G_14, beta1 = 0 over F2 and Q")
# domination (collapse) certificate at the death thresholds
for M, Dm, td in ((13, D13, delta), (14, D14, chi(0.37))):
    nb = lambda v: {u for u in range(n) if Dm[u, v] <= td + 1e-15}
    check(nb(6) <= nb(4), f"M={M}, t=death: closed neighbourhood of x7 {sorted(u+1 for u in nb(6))} lies in that of x5 (x7 dominated by x5)")
    check(nb(5) - {6} <= nb(4), f"M={M}, t=death: after deleting x7, x6 is dominated by x5 (closed nbhd {sorted(u+1 for u in nb(5)-{6})})")
    check(all(Dm[u, v] <= td for u, v in itertools.combinations(range(5), 2)), f"M={M}: x1..x5 form a clique (contractible simplex)")
# triangle side lengths in the certificate
print(f"  triangle x5x6x7 at M=13: sides {c[4,5]:.6f}, {c[5,6]:.6f}, {D13[4,6]:.6f}")
print(f"  triangle x1x5x7 at M=13: sides {c[0,4]:.6f}, {c[0,6]:.6f}, {D13[4,6]:.6f}")
# arc support
for M in (13, 14):
    E, D = results_C[M][0], results_C[M][1]
    sup = all(edge_supported(n, E, c, D, e) for e in E)
    check(sup, f"G_{M} is arc-supported (every edge passes Definition def:arc-support)")
check(not all(edge_supported(n, results_C[12][0], c, results_C[12][1], e) for e in results_C[12][0]) or True,
      "(arc support is only defined for graphs containing the ring; M=12 not tested)")
# matchings
print("bottleneck matchings (raw / normalised by Delta*=0.499):")
P_star = (bs - gs) / (2 * Dl)
for M, bar in ((13, bars13[0]), (14, bars14[0])):
    dr, dg_, eb = matching_costs(bar, (gs, bs))
    print(f"  M={M}: direct = max(|{gs}-{bar[0]:.7f}|, |{bs}-{bar[1]:.7f}|) = {dr:.8f} -> {dr/Dl:.7f};"
          f" diagonal = max({(bar[1]-bar[0])/2:.8f}, {(bs-gs)/2:.8f}) = {dg_:.8f} -> {dg_/Dl:.7f}; e_PH = {eb/Dl:.7f}")
    results_C[M] = results_C[M] + ((dr / Dl, dg_ / Dl, eb / Dl),)
d13, g13, e13 = results_C[13][4]; d14, g14, e14 = results_C[14][4]
check(abs(e13 - 0.0440643) < 5e-8, f"e_PH(M=13) = {e13:.10f} ~ 0.0440643 (manuscript value)")
check(abs(e14 - 0.0300601) < 5e-8 and abs(e14 - P_star) < 1e-15, f"e_PH(M=14) = {e14:.10f} = P* = {P_star:.10f} ~ 0.0300601 (manuscript value)")
check(abs(d13 - 0.122371) < 5e-7, f"direct cost M=13 = {d13:.10f} ~ 0.122371 (manuscript value)")
check(abs(d14 - 0.156051) < 5e-7, f"direct cost M=14 = {d14:.10f} ~ 0.156051 (manuscript value)")
print(f"  raw quantities: g-chi(g) = {gs-B:.8f}; b*-delta = {bs-delta:.8f}; b*-chi(b*) = {bs-chi(0.37):.8f};"
      f" (delta-chi(g))/2 = {(delta-B)/2:.8f}; (chi(b*)-chi(g))/2 = {(chi(0.37)-B)/2:.8f}; (b*-g)/2 = {(bs-gs)/2:.8f}")
for M, dd in ((13, delta), (14, chi(0.37))):
    print(f"  table row M={M}: d={dd:.6f}  |0.37-d|={abs(bs-dd):.6f}  d-chi(0.34)={dd-B:.6f}  0.03/0.998={0.03/0.998:.7f}  2/pi-1/2={2/math.pi-0.5:.4f}")
check((delta - B) / 2 <= 0.021989 < gs - B and max(0.015, (delta - B) / 2) == (delta - B) / 2,
      f"App C hypothesis paragraph: diagonal raw cost max(0.015,(d-chi(0.34))/2) <= (delta-chi(0.34))/2 = {(delta-B)/2:.7f} <= 0.021989 < 0.061063")
ceiling = max(gs - B, bs - chi(0.37)) / Dl
print(f"  cost ceiling (Corollary cor:cost-ceiling) = {ceiling:.10f}")
check(abs(ceiling - d14) < 1e-15, "direct cost at M=14 equals the cost-ceiling bound")
check(e13 > P_star, "e_PH at closure exceeds P* (worse than the empty diagram at M=12)")
# cross-check with release helper if present
try:
    import os
    _here = os.path.dirname(os.path.abspath(__file__))
    for _c in (_here, os.path.join(os.getcwd(), "code")):
        if os.path.exists(os.path.join(_c, "finite_diagram.py")):
            sys.path.insert(0, _c); break
    import finite_diagram as fd
    for M in (13, 14):
        bars_rel, _ = fd.reduction_h1(results_C[M][1])
        md = fd.matching_distance(np.array(results_C[M][2]), np.array([[gs, bs]]))
        check(np.allclose(bars_rel, results_C[M][2], atol=1e-12) and abs(md / Dl - results_C[M][4][2]) < 1e-12,
              f"release finite_diagram.py agrees at M={M} (bars and bottleneck)")
except Exception as ex:
    print("  release helper not run:", ex)
g0 = 2 / math.pi - 0.5
check(gs > g0, f"small-gap hypothesis g < 2/pi-1/2 = {g0:.6f} fails (g = 0.34)")
check(gs > 1 / 3, "g > 1/3: positivity of the intrinsic bar is sample-specific (Lemma circle-bar's sufficient condition fails)")
check(gs - B > (bs - gs) / 2, f"key inequality of the corollary's proof fails: g-chi(g)={gs-B:.6f} > (b*-g)/2={(bs-gs)/2:.6f}")
lam = 2 / math.pi
check(not (B < lam / 3), f"geometric survival bound B < 2/(3pi) = {lam/3:.6f} fails (B = {B:.6f}); survival via chi(g)<=t<chi(b*) instead")

# ================================================================ Part 2
print("\n" + "=" * 78); print("PART 2. Independent-key model (Section 7.2)"); print("=" * 78)
def T_OR(m):
    if m == 0: return Fraction(1)
    return Fraction(2, 2 ** m) - 2 * sum(Fraction(math.comb(m + q - 1, q), 3 ** (m + q)) for q in range(m))
def T_MU(m):
    if m == 0: return Fraction(1)
    return Fraction(2, 2 ** m) - T_OR(m)
# direct check against joint pmf Pr(L=l,R=r) = C(l+r,l) 3^{-(l+r+1)}
def joint(l, r): return Fraction(math.comb(l + r, l), 3 ** (l + r + 1))
K = 80
tot = sum(joint(l, r) for l in range(K) for r in range(K))
print(f"  joint pmf mass on [0,{K})^2 = 1 - {float(1-tot):.3e}")
for m in range(1, 9):
    direct_or = sum(joint(l, r) for l in range(m, K) for r in range(m, K))
    direct_mu = 1 - sum(joint(l, r) for l in range(m) for r in range(m))
    trunc = 2 * Fraction(1, 2 ** K)   # Pr(L>=K)+Pr(R>=K)
    check(abs(direct_or - T_OR(m)) <= trunc and abs(direct_mu - T_MU(m)) <= trunc,
          f"m={m}: T_OR={float(T_OR(m)):.12f}, T_mutual={float(T_MU(m)):.12f} match direct joint-pmf sums")
check(all(Fraction(1, 2 ** m) == sum(Fraction(1, 2 ** (l + 1)) for l in range(m, 400)) + Fraction(1, 2 ** 400) for m in (1, 3, 7)),
      "Pr(L=l)=2^{-(l+1)}, so Pr(L>=m)=2^{-m}")
# integral form by quadrature
try:
    from scipy.integrate import quad
    mx = 0
    for m in range(1, 13):
        f = lambda t: 2 * t ** (m - 1) * math.exp(-2 * t) / math.factorial(m - 1) * (1 - math.exp(-t) * sum(t ** q / math.factorial(q) for q in range(m)))
        val = quad(f, 0, np.inf, epsabs=1e-15, epsrel=1e-13, limit=200)[0]
        mx = max(mx, abs(val - float(T_OR(m))))
    check(mx < 1e-12, f"integral representation equals closed form for m=1..12 (max diff {mx:.2e})")
except Exception as ex:
    print("  quadrature skipped:", ex)
# binomial positive form
ok = all(T_OR(m) == Fraction(2, 2 ** m) * sum(Fraction(math.comb(2 * m - 1, j) * 2 ** j, 3 ** (2 * m - 1)) for j in range(m)) for m in range(1, 30))
check(ok, "positive form T_OR(m) = 2^{1-m} Pr{Bin(2m-1,2/3) <= m-1} (m=1..29, exact)")
check(T_OR(1) == Fraction(1, 3) and T_MU(1) == Fraction(2, 3) and T_OR(2) == Fraction(7, 54) and T_MU(2) == Fraction(10, 27),
      "T_OR(1)=1/3, T_mutual(1)=2/3, T_OR(2)=7/54, T_mutual(2)=10/27")
check(all(T_OR(m) > 0 and T_MU(m) > 0 and T_OR(m) < T_MU(m) for m in range(1, 40)), "0 < T_OR(m) < T_mutual(m) for m=1..39")
# tie probability
getcontext().prec = 80
def tie_prob(T, n, mmax=400):
    Td = [Decimal(T(m).numerator) / Decimal(T(m).denominator) for m in range(mmax + 2)]
    s_ = Decimal(0)
    for m in range(1, mmax + 1):
        pm = Td[m] - Td[m + 1]
        s_ += pm * (1 - Td[m]) ** (n - 1)
    trunc = n * Td[mmax + 1]
    return 1 - n * s_, trunc
Tor_cache = {}; Tmu_cache = {}
def Tor(m):
    if m not in Tor_cache: Tor_cache[m] = T_OR(m)
    return Tor_cache[m]
def Tmu(m):
    if m not in Tmu_cache: Tmu_cache[m] = Fraction(2, 2 ** m) - Tor(m) if m else Fraction(1)
    return Tmu_cache[m]
po, tro = tie_prob(Tor, 512, 200); pmu, trm = tie_prob(Tmu, 512, 200)
print(f"  n=512: Pr(max OR mark tied)     = {str(po)[:22]}  (truncation bound {float(tro):.1e})")
print(f"  n=512: Pr(max mutual mark tied) = {str(pmu)[:22]}  (truncation bound {float(trm):.1e})")
check(round(float(po), 6) == 0.328768, f"OR tie probability rounds to 0.328768 (manuscript value); digits {str(po)[:12]}")
check(round(float(pmu), 6) == 0.274923, f"mutual tie probability rounds to 0.274923 (manuscript value); digits {str(pmu)[:12]} -- NOTE: truncation would give 0.274922")
check(round(float(po), 3) == 0.329 and round(float(pmu), 3) == 0.275, "rounded 0.329 and 0.275 as in Section 7.2")
for nn in (128, 4096, 32768):
    a_, _ = tie_prob(Tor, nn, 200); b_, _ = tie_prob(Tmu, nn, 200)
    print(f"  n={nn}: OR {float(a_):.6f}, mutual {float(b_):.6f}")
# Monte Carlo sanity check of the tie formula under the model
rng = np.random.default_rng(20261020)
reps, nn = 4000, 512
T = rng.exponential(1.0, size=(reps, nn)); L = rng.poisson(T); R = rng.poisson(T)
for name, Z, p in (("OR", np.minimum(L, R), float(po)), ("mutual", np.maximum(L, R), float(pmu))):
    mx_ = Z.max(axis=1); ties = (Z == mx_[:, None]).sum(axis=1) >= 2
    est = ties.mean(); se = math.sqrt(p * (1 - p) / reps)
    check(abs(est - p) < 4 * se, f"Monte Carlo ({reps} reps) {name} tie freq {est:.4f} vs exact {p:.4f} (4 SE = {4*se:.4f})")
# asymptotics
for m in (10, 50, 200, 400):
    r1 = float(Fraction(T_OR(m)) / Fraction(1)) / (3 * (4 / 9) ** m / math.sqrt(math.pi * m)) if m <= 200 else None
    if r1 is not None:
        r2 = float(T_MU(m) * 2 ** (m - 1))
        print(f"  m={m}: T_OR / [3(4/9)^m/sqrt(pi m)] = {r1:.6f};  T_mutual / 2^(1-m) = {r2:.12f}")
def ratio_or(m):
    # T_OR(m)/[3(4/9)^m/sqrt(pi m)] computed exactly up to the final float division
    x = T_OR(m) * Fraction(9 ** m, 3 * 4 ** m)
    return float(x) * math.sqrt(math.pi * m)
rs = [(m, ratio_or(m)) for m in (100, 200, 400)]
for m, r in rs: print(f"  m={m}: OR ratio {r:.8f}, m*(1-ratio) = {m*(1-r):.5f}")
check(all(r1 < r2 < 1 for (_, r1), (_, r2) in zip(rs, rs[1:])) and abs(rs[-1][1] - 1) < 2e-2,
      "OR tail ratio increases to 1 with 1-ratio = O(1/m) (slow, as expected from Stirling corrections)")
check(all(Fraction(j, 2 * (2 * m - j)) <= Fraction(1, 2) for m in range(1, 60) for j in range(1, m)),
      "binomial mass ratio a_{m,j-1}/a_{m,j} = j/(2(2m-j)) <= 1/2 for 1<=j<=m-1")
ok = all(Fraction(math.comb(2*m-1, j-1) * 2**(j-1), 1) / Fraction(math.comb(2*m-1, j) * 2**j, 1) == Fraction(j, 2*(2*m-j)) for m in range(2, 30) for j in range(1, m))
check(ok, "ratio formula exact for Bin(2m-1,2/3)")

# ================================================================ Part 3
print("\n" + "=" * 78); print("PART 3. Arc-support figure (Figure 2)"); print("=" * 78)
s6 = [0, 0.051, 0.442, 0.491, 0.543, 0.953]; n6 = 6
c6 = np.zeros((n6, n6))
for i, j in itertools.combinations(range(n6), 2):
    c6[i, j] = c6[j, i] = chi(dstar(s6, i, j))
ring6 = [(i, (i + 1) % n6) if i < (i + 1) % n6 else ((i + 1) % n6, i) for i in range(n6)]
G6 = ring6 + [(0, 3)]
D6 = shortest_paths(n6, G6, c6)
hh = c6[0, 3]
print(f"  h = c(x1,x4) = chi(0.491) = {hh:.8f}")
check(0.31818 < hh < 0.31819, "h in (0.31818, 0.31819)")
for (i, j), bound in (((0, 2), 0.35060), ((0, 4), 0.35249), ((3, 1), 0.34862), ((3, 5), 0.35742)):
    check(D6[i, j] > bound, f"D_G(x{i+1},x{j+1}) = {D6[i,j]:.6f} > {bound}")
print(f"  smallest of the four failed support distances: {min(D6[0,2],D6[0,4],D6[3,1],D6[3,5]):.6f}")
check(min(D6[0, 2], D6[0, 4], D6[3, 1], D6[3, 5]) > 0.348 > hh, "all four exceed 0.348 > h")
res = arc_support(n6, G6, c6, D6, 0, 3)
check(not any(ok for *_, ok, _ in res), "x1x4 fails arc support on both arcs from both endpoints")
check(all(edge_supported(n6, G6, c6, D6, e) for e in ring6), "ring edges pass")
Rr = np.zeros((n6, n6))
DR = shortest_paths(n6, ring6, c6)
ok = all(abs(D6[i, j] - min(DR[i, j], DR[i, 0] + hh + DR[3, j], DR[i, 3] + hh + DR[0, j])) < 1e-14 for i in range(n6) for j in range(n6))
check(ok, "Remark 2 formula D_G(i,j)=min{R_ij, R_i1+h+R_4j, R_i4+h+R_1j} holds for all pairs")
E33, T33, b2, bq = betti1_at(D6, 0.33)
print("  Rips edges at t=0.33:", [f"x{i+1}x{j+1}({D6[i,j]:.4f})" for i, j in E33])
print("  triangles at t=0.33:", [f"x{a+1}x{b+1}x{c_+1}" for a, b, c_ in T33])
check(len(E33) == 9 and set(T33) == {(0, 1, 5), (2, 3, 4)}, "t=0.33: nine edges; exactly the triangles x1x2x6 and x3x4x5")
check(set(E33) - set(G6) == {(1, 5), (2, 4)}, "the two Rips-only pairs at t=0.33 are x2x6 and x3x5 (not graph edges)")
check(b2 == 2 and bq == 2, "beta1 = 2 at t=0.33 over F2 and Q")
check(all(not any(set(e) <= set(tr) for tr in T33) for e in [(0, 3), (1, 2), (4, 5)]), "joining edges x1x4, x2x3, x5x6 lie in no triangle")
nxt = sorted(D6[i, j] for i, j in itertools.combinations(range(n6), 2) if D6[i, j] > 0.33)[0]
lst = sorted(D6[i, j] for i, j in itertools.combinations(range(n6), 2) if D6[i, j] <= 0.33)[-1]
print(f"  threshold graph is constant on [{lst:.6f}, {nxt:.6f}); 0.33 lies inside")
print(f"  ring chords: x2x3 = chi(0.391) = {c6[1,2]:.6f}, x5x6 = chi(0.41) = {c6[4,5]:.6f}")
angs = {s_: (90 - 360 * s_) % 360 for s_ in s6}
print("  inset label angles 90-360s (mod 360):", {k: round(v if v <= 180 else v - 360, 2) for k, v in angs.items()})
check(all(abs(((90 - 360 * s_ + 180) % 360) - 180 - a_) < 1e-9 for s_, a_ in ((0.051, 71.64), (0.442, -69.12), (0.491, -86.76), (0.543, -105.48), (0.953, 106.92))),
      "hard-coded inset label angles 71.64, -69.12, -86.76, -105.48, 106.92 equal 90-360s")
# supported panel (drawn configuration): five points w0..w4 on an arc
sw = [0, 0.07, 0.15, 0.23, 0.30]
cw = lambda i, j: chi(abs(sw[i] - sw[j]))
cab = cw(0, 4)
path = [sum(cw(k, k + 1) for k in range(i)) for i in range(5)]
print(f"  supported panel: arc coords {sw}; c(a,b) = chi(0.30) = {cab:.6f}")
print("   ring-path distances from p=a=w0:", [round(x, 6) for x in path])
check(all(cw(k, k + 1) <= cab for k in range(4)), "panel (a): every consecutive chord <= c(a,b)")
check(all(min(path[i], cab + (path[4] - path[i])) <= cab + 1e-15 for i in range(5)),
      "panel (a): D_G(p,w_i) <= c(a,b) for all i, already within ring + ab, so the dashed spokes pw2, pw3 are Rips pairs at t=c(a,b)")
check(path[2] > cw(0, 2) and path[3] > cw(0, 3), "panel (a): spokes pw2, pw3 are not graph edges (their D_G are ring paths, not chords)")

# ================================================================ Part 4
print("\n" + "=" * 78); print("PART 4. Six-point ellipse (Figure 3, Example 2)"); print("=" * 78)
th = [30, 85, 150, 255, 290, 335]
X = np.array([[2 * math.cos(math.radians(t)), math.sin(math.radians(t))] for t in th])
for i, p in enumerate(X): print(f"  x{i} = ({p[0]:.6f}, {p[1]:.6f})")
ne = 6
ce = np.array([[np.linalg.norm(X[i] - X[j]) for j in range(ne)] for i in range(ne)])
pairs = list(itertools.combinations(range(ne), 2))
chords = sorted(ce[p] for p in pairs)
mg = min(b - a for a, b in zip(chords, chords[1:]))
check(mg > 0.03, f"fifteen chords distinct; smallest gap {mg:.4f}")
ringE = [(i, (i + 1) % ne) if i < (i + 1) % ne else ((i + 1) % ne, i) for i in range(ne)]
# ranks
def rank(i, j):
    others = sorted((k for k in range(ne) if k != i), key=lambda k: (ce[i, k], k))
    return 1 + others.index(j)
key_or = {p: min(rank(*p), rank(p[1], p[0])) for p in pairs}
key_mu = {p: max(rank(*p), rank(p[1], p[0])) for p in pairs}
ord_rad = sorted(pairs, key=lambda p: (ce[p], p))
ord_or = sorted(pairs, key=lambda p: (key_or[p], ce[p], p))
ord_mu = sorted(pairs, key=lambda p: (key_mu[p], ce[p], p))
Gset = set(ringE) | {(0, 4), (1, 4)}
for name, o in (("radius", ord_rad), ("OR", ord_or), ("mutual", ord_mu)):
    clos = max(o.index(e) + 1 for e in ringE)
    print(f"  {name}: closure at insertion {clos}; 8th prefix = {sorted(o[:8])}")
    check(set(o[:8]) == Gset, f"{name}: eighth prefix is R + x0x4 + x1x4")
    check(clos == {"radius": 7, "OR": 6, "mutual": 6}[name], f"{name}: closure position {clos}")
Ge = sorted(Gset)
De = shortest_paths(ne, Ge, ce)
print(f"  c04 = {ce[0,4]:.6f}, c14 = {ce[1,4]:.6f}, c15 = {ce[1,5]:.6f}, c12 = {ce[1,2]:.6f}")
check(1.7807 < ce[0, 4] < 1.7808, "1.7807 < c04 < 1.7808")
check(2.0018 < ce[1, 4] < 2.0019, "2.0018 < c14 < 2.0019")
check(ce[1, 5] > 2.1672 > ce[1, 4], f"c15 = {ce[1,5]:.5f} > 2.1672 > c14")
B2 = max(ce[e] for e in ringE)
check(abs(B2 - ce[1, 2]) < 1e-15 and 1.9698 < B2 < 1.9699, f"B = c12 = {B2:.6f} in (1.9698, 1.9699)")
# support of x0x4 via x4,x5,x0 and x1x4 via x4,x5,x0,x1 from x4
check(max(ce[4, 5], ce[5, 0], De[4, 5]) < 1.242 < ce[0, 4], f"x0x4: consecutive costs {ce[4,5]:.4f},{ce[5,0]:.4f} and D(x4,x5) below 1.242")
check(max(ce[4, 5], ce[5, 0], ce[0, 1]) < 1.635 and max(De[4, 5], De[4, 0]) < 1.781 < ce[1, 4],
      f"x1x4: consecutive costs <= {max(ce[4,5],ce[5,0],ce[0,1]):.4f} < 1.635; interior distances {De[4,5]:.4f},{De[4,0]:.4f} < 1.781")
check(all(edge_supported(ne, Ge, ce, De, e) for e in Ge), "G is arc-supported (full check)")
E21, T21, b2, bq = betti1_at(De, 2.1)
check(set(E21) == Gset, "t=2.1: threshold graph is exactly G (8 edges)")
check(set(T21) == {(0, 1, 4), (0, 4, 5)}, "t=2.1: exactly the triangles x0x1x4 and x0x4x5")
om = min(De[p] for p in pairs if p not in Gset)
check(om > 2.443, f"omitted pairs have D_G > 2.443 (min {om:.4f})")
check(b2 == 1 and bq == 1, "beta1 = 1 at t=2.1 over F2 and Q")
check(all(not any(set(e) <= set(tr) for tr in T21) for e in [(1, 2), (2, 3), (3, 4)]), "x1x2, x2x3, x3x4 in no triangle (forced block)")
check(B2 < 2.1, "t=2.1 > B, inside the positive bar")
bars_e = f2_reduction_h1(De); rbe = ripser_h1(De)
print(f"  H1 bars of D_G (own F2): {bars_e}; ripser: {rbe}")
check(len(bars_e) == 1 and bars_e[0][0] < 2.1 < bars_e[0][1], "single H1 bar containing t=2.1")
# arc lengths / minor arc
from scipy.integrate import quad
sp = lambda u: math.sqrt(4 * math.sin(u) ** 2 + math.cos(u) ** 2)
Lp = quad(sp, 0, 2 * math.pi, epsabs=1e-13)[0]
arc41 = quad(sp, math.radians(290), math.radians(445), epsabs=1e-13)[0]
print(f"  perimeter L = {Lp:.6f}; arc x4->x5->x0->x1 = {arc41:.6f} (L/2 = {Lp/2:.6f})")
check(arc41 < Lp / 2, "x4->x5->x0->x1 is the minor arc from x4 to x1")
gaps_e = [quad(sp, math.radians(th[i]), math.radians(th[(i + 1) % 6] + (360 if i == 5 else 0)))[0] for i in range(6)]
check(max(gaps_e) < 2.554 and max(gaps_e) / Lp < 1 / 3, f"largest arc gap {max(gaps_e):.4f} < 2.554; normalized {max(gaps_e)/Lp:.4f} < 1/3")
check(Lp >= 4 * math.sqrt(5), "L >= 4 sqrt 5")
print("\nSUMMARY:", "ALL CHECKS PASSED" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}")
sys.exit(1 if FAIL else 0)
