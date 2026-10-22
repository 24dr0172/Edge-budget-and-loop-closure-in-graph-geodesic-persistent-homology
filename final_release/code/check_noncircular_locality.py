"""
Numerical check of the locality proposition on two non-circular curves
(Appendix A of the manuscript). Additional evidence: no theorem depends on it.

Curves, normalized to length one:
  ellipse : (2 cos t, sin t)
  peanut  : r = 1 + 0.4 cos 2t, (r cos t, r sin t)   (non-convex, with a waist)
A round circle is run as a control for the code only (on by default; --no-control skips it).

n = 512 points uniform in arc length, 20 samples per curve; radius, OR-kNN and
mutual-kNN orders with the manuscript's pair ordering (key, chord, endpoint labels).

r_loc = min(nu, beta) of Appendix A is evaluated on a grid (nu_beta below). It is a
numerical value, not a certified one: nu is found from sign changes of
(gamma(s)-gamma(t)).T(t) on a grid and refined by linear interpolation, so a tangential
(double) zero would be missed. Two grid sizes are printed for convergence.

Per sample and order:
  loc_ok        the proposition's locality hypothesis for short-hops-first holds:
                radius: B < r_loc; kNN: with k_c the largest consecutive key, every
                vertex's distance to its k_c-th nearest neighbour is < r_loc
  shh_ok        short-hops-first holds for every edge inserted through ring closure
                (Definition 3, checked edge by edge, independent of loc_ok)
  conn_close    closure at connectivity: exact prefix for radius (M_conn == M_close);
                integer k for OR and mutual (k_conn == k_close), computed directly
  closed_sup    the prefix at M_close is arc-supported (Definition 4 with that prefix's
                distances)
  M_end         end of the INITIAL UNINTERRUPTED RUN of arc-supported prefixes that
                starts at M_close: the largest M such that every prefix M_close..M is
                arc-supported (0 if the prefix at M_close is not). Arc support need
                not be monotone in M beyond that point: an edge that fails at its
                insertion can become supported after later edges, and no claim is made
                about prefixes after the run ends.
  fail_chord    chord length of the first edge that breaks the run (nan if the run
                reaches the complete graph)
  bars_*        ALL H1 bars returned by Ripser (float32 endpoints, no lifetime cutoff,
                the release's convention) at M_close and at M_end; ring_* counts bars
                born at B; other_max_* is the longest other bar (0 if none)

Within the run, a newly inserted edge is checked against the distances of the prefix
that contains it. Edges already supported stay supported when edges are added
(distances only shrink), so the prefix is supported exactly when every new edge was.

Usage (from the release folder; about 10 minutes on a multi-core machine):
    python3 code/check_noncircular_locality.py [outdir] [--no-control] [--samples=k] [--workers=k]
Defaults to results/noncircular_locality/ including the round-circle control, the saved configuration.
Seeds are fixed, so a rerun reproduces the saved CSV.
"""
import csv
import json
import os
import sys
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from ripser import ripser

ROOT = Path(__file__).resolve().parents[1]
N = 512
SAMPLES = 20


def curve(kind, theta):
    if kind == "ellipse":
        return np.c_[2 * np.cos(theta), np.sin(theta)]
    if kind == "peanut":
        r = 1 + 0.4 * np.cos(2 * theta)
        return np.c_[r * np.cos(theta), r * np.sin(theta)]
    return np.c_[np.cos(theta), np.sin(theta)]


def curve_d(kind, theta):
    if kind == "ellipse":
        return np.c_[-2 * np.sin(theta), np.cos(theta)]
    if kind == "peanut":
        r = 1 + 0.4 * np.cos(2 * theta)
        dr = -0.8 * np.sin(2 * theta)
        return np.c_[dr * np.cos(theta) - r * np.sin(theta), dr * np.sin(theta) + r * np.cos(theta)]
    return np.c_[-np.sin(theta), np.cos(theta)]


def arclength_table(kind):
    t = np.linspace(0, 2 * np.pi, 400001)
    Q = curve(kind, t)
    cum = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))]
    return t, cum


def nu_beta(kind, grid):
    """Grid evaluation of nu and beta of Appendix A for the unit-length curve."""
    t, cum = arclength_table(kind)
    L = cum[-1]
    s = np.linspace(0, L, grid, endpoint=False)
    th = np.interp(s, cum, t)
    P = curve(kind, th) / L
    T = curve_d(kind, th)
    T /= np.linalg.norm(T, axis=1)[:, None]
    P2 = np.vstack([P, P[0]])
    nu = np.inf
    for j in range(grid):
        h = (P - P[j]) @ T[j]
        h2 = np.r_[h, h[0]]
        cr = np.flatnonzero(np.sign(h2[:-1]) * np.sign(h2[1:]) < 0)
        cr = cr[(cr != j) & (cr != (j - 1) % grid)]       # sign change at s = t itself
        for c in cr:
            lam = h2[c] / (h2[c] - h2[c + 1])
            pt = P2[c] + lam * (P2[c + 1] - P2[c])
            nu = min(nu, float(np.linalg.norm(pt - P[j])))
    beta = min(float(np.linalg.norm(P[i] - P[(i + grid // 2) % grid])) for i in range(grid))
    return nu, beta


def sample(kind, rng):
    t, cum = arclength_table(kind)
    L = cum[-1]
    s = np.sort(rng.random(N))
    theta = np.interp(s * L, cum, t)
    P = curve(kind, theta) / L
    return s, P


def orderings(C):
    n = len(C)
    iu, ju = np.triu_indices(n, 1)
    R = np.zeros((n, n), int)
    for i in range(n):
        idx = np.argsort(C[i], kind="stable")
        idx = idx[idx != i]
        R[i, idx] = np.arange(1, n)
    c = C[iu, ju]
    out = {"radius": np.lexsort((ju, iu, c))}
    keys = {"OR": np.minimum(R[iu, ju], R[ju, iu]), "mutual": np.maximum(R[iu, ju], R[ju, iu])}
    for name, key in keys.items():
        out[name] = np.lexsort((ju, iu, c, key))
    return iu, ju, out, keys


def minor_arc(a, b, s, n):
    """Samples a=w_0,...,w_m=b along the shorter arc."""
    if s[b] - s[a] < 0.5:
        return np.arange(a, b + 1)
    return np.r_[np.arange(a, -1, -1), np.arange(n - 1, b - 1, -1)]


def short_hops_first(C, pos, iu, ju, order, m_close, s):
    n = len(C)
    for e in order[:m_close + 1]:
        a, b = int(iu[e]), int(ju[e])
        seq = minor_arc(a, b, s, n)
        if len(seq) <= 2:
            continue
        pab, cab = pos[a, b], C[a, b]
        if np.any(pos[seq[:-1], seq[1:]] >= pab):
            return False
        cc = C[seq[:-1], seq[1:]]
        cum = np.r_[0.0, np.cumsum(cc)]
        mid = seq[1:-1]
        t = np.arange(1, len(seq) - 1)
        left = ((pos[a, mid] < pab) & (C[a, mid] <= cab)) | (cum[t] <= cab)
        right = ((pos[mid, b] < pab) & (C[mid, b] <= cab)) | (cum[-1] - cum[t] <= cab)
        if not (left.all() and right.all()):
            return False
    return True


def edge_supported(D, C, rs, a, b, n):
    cab = C[a, b]
    arcs = ((np.arange(a, b + 1), rs[a:b]),
            (np.r_[np.arange(b, n), np.arange(0, a + 1)], np.r_[rs[b:], rs[:a]]))
    for seq, steps in arcs:
        if steps.max() > cab:
            continue
        if D[a, seq].max() <= cab or D[b, seq].max() <= cab:
            return True
    return False


def graph_dist(C, edges_a, edges_b):
    n = len(C)
    # A dense csgraph conversion treats weights below about 1e-8 as null.
    # Preserve every positive chord, including near-coincident sample pairs.
    weights = C[edges_a, edges_b]
    W = csr_matrix((np.r_[weights, weights],
                    (np.r_[edges_a, edges_b], np.r_[edges_b, edges_a])), shape=(n, n))
    return shortest_path(W, method="D", directed=False)


def h1_bars(D):
    """All H1 bars Ripser returns (it omits zero-length pairs); no extra cutoff."""
    dg = ripser(D, distance_matrix=True, maxdim=1)["dgms"][1]
    return [(float(b), float(d)) for b, d in dg if d > b]


def bar_summary(bars, B):
    ring = [x for x in bars if abs(x[0] - B) <= 1e-6 * max(1.0, B)]
    other = [x for x in bars if x not in ring]
    return dict(n=len(bars), ring=len(ring),
                other_max=max((d - b for b, d in other), default=0.0))


def components_at(order, iu, ju, n):
    """Position (0-based) of the edge that first connects the graph."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    comps = n
    for p, e in enumerate(order):
        ra, rb = find(int(iu[e])), find(int(ju[e]))
        if ra != rb:
            parent[ra] = rb
            comps -= 1
            if comps == 1:
                return p
    return len(order) - 1


def run_sample(kind, seed, r_loc):
    rng = np.random.default_rng(seed)
    s, P = sample(kind, rng)
    n = N
    C = np.linalg.norm(P[:, None] - P[None], axis=2)
    gaps = np.diff(np.r_[s, s[0] + 1])
    iu, ju, orders, keys = orderings(C)
    rs = np.array([C[i, (i + 1) % n] for i in range(n)])
    B = rs.max()
    ring = [(i, i + 1) for i in range(n - 1)] + [(0, n - 1)]
    pid = np.zeros((n, n), int)
    pid[iu, ju] = np.arange(len(iu))
    Csort = np.sort(C, axis=1)                      # column 0 is the vertex itself
    row = {"curve": kind, "seed": seed, "g": float(gaps.max()), "B": float(B)}
    for name, order in orders.items():
        pos = np.zeros((n, n), int)
        pos[iu[order], ju[order]] = np.arange(len(order))
        pos[ju[order], iu[order]] = np.arange(len(order))
        m_close = max(pos[i, j] for i, j in ring)
        m_conn = components_at(order, iu, ju, n)
        r = {"M_close": int(m_close) + 1, "M_conn": int(m_conn) + 1}
        if name == "radius":
            r["conn_close"] = bool(m_conn == m_close)
            r["loc_ok"] = bool(B < r_loc)
        else:
            k = keys[name]
            kc = int(k[order[m_conn]])
            kcl = int(max(k[pid[i, j]] for i, j in ring))
            r["conn_close"] = bool(kc == kcl)
            r["k_conn"], r["k_close"] = kc, kcl
            r["loc_ok"] = bool(Csort[:, kcl].max() < r_loc)   # k_c-th nearest neighbour
        r["shh_ok"] = bool(short_hops_first(C, pos, iu, ju, order, m_close, s))
        ea, eb = iu[order[:m_close + 1]], ju[order[:m_close + 1]]
        D = graph_dist(C, ea, eb)
        sup = all(edge_supported(D, C, rs, int(a), int(b), n) for a, b in zip(ea, eb))
        r["closed_sup"] = bool(sup)
        bc = bar_summary(h1_bars(D), B)
        r["bars_close"], r["ring_close"], r["other_max_close"] = bc["n"], bc["ring"], bc["other_max"]
        end, fail_chord = m_close, float("nan")
        if sup:
            total = len(order)
            for m in range(m_close + 1, total):
                a, b = int(iu[order[m]]), int(ju[order[m]])
                w = C[a, b]
                if w < D[a, b]:
                    # only rows Ia (closer via a) and columns Jb (closer via b) can improve
                    da, db = D[:, a].copy(), D[:, b].copy()
                    Ia, Jb = np.flatnonzero(da + w < db), np.flatnonzero(db + w < da)
                    blk = np.minimum(D[np.ix_(Ia, Jb)], da[Ia][:, None] + w + db[Jb][None, :])
                    D[np.ix_(Ia, Jb)] = blk
                    D[np.ix_(Jb, Ia)] = blk.T
                if not edge_supported(D, C, rs, a, b, n):
                    fail_chord = float(w)
                    break
                end = m
            # D already contains the edge that broke the run; rebuild at M_end
            D = graph_dist(C, iu[order[:end + 1]], ju[order[:end + 1]])
            be = bar_summary(h1_bars(D), B)
            r["bars_end"], r["ring_end"], r["other_max_end"] = be["n"], be["ring"], be["other_max"]
            r["M_end"] = int(end) + 1
            r["fail_chord"] = fail_chord
        else:
            r["bars_end"], r["ring_end"], r["other_max_end"], r["M_end"] = -1, 0, float("nan"), 0
            r["fail_chord"] = float("nan")
        for kx, v in r.items():
            row[f"{name}_{kx}"] = v
    return row


def _task(args):
    return run_sample(*args)


def main():
    global SAMPLES
    workers = None
    arguments = list(sys.argv[1:])
    # A separate --workers value is an option argument, not an output path.
    if "--workers" in arguments:
        index = arguments.index("--workers")
        if index+1 >= len(arguments):
            raise SystemExit("--workers requires a positive integer")
        del arguments[index:index+2]
    args = [a for a in arguments if not a.startswith("--")]
    out = args[0] if args else str(ROOT / "results/noncircular_locality")
    Path(out).mkdir(parents=True, exist_ok=True)
    for a in sys.argv[1:]:
        if a.startswith("--samples="):
            SAMPLES = int(a.split("=")[1])
        if a.startswith("--workers="):
            workers = int(a.split("=")[1])
    if "--workers" in sys.argv:
        index = sys.argv.index("--workers")
        if index+1 >= len(sys.argv):
            raise SystemExit("--workers requires a positive integer")
        workers = int(sys.argv[index+1])
    if SAMPLES < 1 or (workers is not None and workers < 1):
        raise SystemExit("samples and workers must be positive")
    jobs = [("ellipse", 20350101), ("peanut", 20350201), ("circle", 20350301)]
    if "--no-control" in sys.argv:
        jobs = jobs[:2]
    loc = {}
    for kind, _ in jobs:
        v = {g: nu_beta(kind, g) for g in (4000, 8000)}
        loc[kind] = dict(grid_values={str(g): dict(nu=x[0], beta=x[1]) for g, x in v.items()},
                         r_loc=min(min(x) for x in v.values()))
        print(kind, "nu, beta by grid:", v, flush=True)
    tasks = [(kind, seed0 + j, loc[kind]["r_loc"]) for kind, seed0 in jobs for j in range(SAMPLES)]
    from multiprocessing import Pool
    rows = []
    pool_size = min(len(tasks), workers) if workers is not None else max(1, min(len(tasks), (os.cpu_count() or 2) // 2))
    with Pool(pool_size) as pool:
        for row in pool.imap_unordered(_task, tasks):
            rows.append(row)
            print(len(rows), "/", len(tasks), row["curve"], row["seed"], flush=True)
    rows.sort(key=lambda r: (r["curve"], r["seed"]))
    with open(f"{out}/locality_check_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    summary = {"r_loc": loc}
    for kind in dict.fromkeys(r["curve"] for r in rows):
        rr = [r for r in rows if r["curve"] == kind]
        for name in ("radius", "OR", "mutual"):
            fc = np.array([r[f"{name}_fail_chord"] for r in rr], float)
            have = ~np.isnan(fc)
            summary[f"{kind}/{name}"] = {
                "samples": len(rr),
                "locality_hypothesis_holds": sum(r[f"{name}_loc_ok"] for r in rr),
                "shh_ok": sum(r[f"{name}_shh_ok"] for r in rr),
                "closure_at_connectivity": sum(r[f"{name}_conn_close"] for r in rr),
                "closed_prefix_supported": sum(r[f"{name}_closed_sup"] for r in rr),
                "runs_reaching_complete_graph": sum(r[f"{name}_M_end"] == N*(N-1)//2 for r in rr),
                "fail_chord_min": float(fc[have].min()) if have.any() else None,
                "fail_chord_median": float(np.median(fc[have])) if have.any() else None,
                "fail_chord_max": float(fc[have].max()) if have.any() else None,
                "exactly_one_bar_at_close": sum(r[f"{name}_bars_close"] == 1 for r in rr),
                "exactly_one_bar_at_end": sum(r[f"{name}_bars_end"] == 1 for r in rr),
                "ring_bar_at_end": sum(r[f"{name}_ring_end"] == 1 for r in rr),
                "max_other_bar_lifetime_at_end": float(max(r[f"{name}_other_max_end"] for r in rr)),
                "max_other_bar_lifetime_at_close": float(max(r[f"{name}_other_max_close"] for r in rr)),
            }
    with open(f"{out}/locality_check_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
