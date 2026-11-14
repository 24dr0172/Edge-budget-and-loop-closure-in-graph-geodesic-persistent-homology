"""Winding-triangle survival check and certified first-nonempty range on the peanut clouds.

This reuses the 100 paired peanut clouds and 300 method records of the
varying-waist study; it draws no new clouds and computes no persistence.

For a graph G containing the ring, W_G is the smallest threshold at which the
graph-metric Rips complex contains a winding triangle (three samples whose
three cyclic arc gaps are all below one half) or an intrinsically antipodal
pair. A ring class born at B survives at least until W_G, so W_G > B
certifies survival. W_E is the same quantity for the chord metric, which is
the complete graph's metric.

For each record the script computes, at exact ring closure:
  W_G, W_E, the exact intrinsic death b* (the smallest largest gap of a
  winding triple), lambda_G and the edgewise survival test lambda_G b* > B;
  the bound of the sharper fixed-graph certificate,
  min{max{g-B, b*-max{B,W_G}}, (b*-B)/2} / Delta*, for arc-supported graphs;
  short-hops-first through closure (Definition 2, recomputed);
  the certified range end U*: the largest exact prefix M >= M_close with
  W_{G_M} > B. W_{G_M} cannot increase as edges are added, so U* is found by
  bisection; U* is the complete graph exactly when W_E > B.
It also recomputes W_G and the edgewise test at the later prefix M_shortcut of
the peanut study, the first prefix after closure with an edge spanning a
quarter of the curve. It then compares the hypotheses with the saved
first-nonempty selections.
"""
from pathlib import Path
import argparse, csv, hashlib, json, os, sys
sys.dont_write_bytecode = True
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np
from graph_methods import build_orders, first_connected, closure, graph_metric, intrinsic, certificate

ROOT = Path(__file__).resolve().parents[1]
FOLDER = 'results/winding_check'
METHODS = ('R2_global', 'OR_kNN', 'mutual_kNN')
ALPHAS = (0., .4, .6, .8, .9)
TOLERANCE = 5e-8
FIELDS = ('alpha', 'seed', 'method', 'n', 'g', 'B', 'b_star', 'b_star_saved', 'diameter', 'P_star',
          'M_conn', 'M_close', 'pairs', 'lambda_G', 'edgewise_survival', 'W_G', 'W_G_triple',
          'W_G_over_lambda_b_star', 'W_E', 'tau_E', 'tau_G', 'closure_supported', 'closure_bars',
          'W_G_gt_B', 'W_E_gt_B', 'tau_E_gt_B', 'antipodal_pairs', 'short_hops_first',
          'certified', 'U_star', 'U_star_is_complete', 'U_star_probes', 'M_first', 'first_at_closure',
          'first_early', 'theorem_applies', 'closure_error', 'euclidean_error', 'U_G', 'winding_bound',
          'winding_bound_applicable', 'M_shortcut', 'lambda_G_shortcut', 'W_G_shortcut',
          'W_G_gt_B_shortcut', 'edgewise_survival_shortcut')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def winding_threshold(D, s):
    """W for metric D on samples with sorted arc positions s; also the minimizing triple.

    A triple i<j<k winds when s_j-s_i, s_k-s_j and 1-(s_k-s_i) are all below 1/2.
    Intrinsically antipodal pairs (separation exactly 1/2) are also included.
    """
    n = len(s)
    best, triple = np.inf, None
    for i in range(n-2):
        j1 = int(np.searchsorted(s, s[i]+.5, side='left'))      # s_j - s_i < 1/2
        k0 = int(np.searchsorted(s, s[i]+.5, side='right'))     # s_k - s_i > 1/2
        J = np.arange(i+1, min(j1, n))
        K = np.arange(max(k0, i+2), n)
        if not len(J) or not len(K):
            continue
        ok = (K[None, :] > J[:, None]) & (s[K][None, :]-s[J][:, None] < .5)
        value = np.maximum(np.maximum(D[i, J][:, None], D[i, K][None, :]), D[np.ix_(J, K)])
        value = np.where(ok, value, np.inf)
        flat = int(np.argmin(value))
        if value.flat[flat] < best:
            a, b = divmod(flat, len(K))
            best, triple = float(value.flat[flat]), (i, int(J[a]), int(K[b]))
    sep = np.abs(s[:, None]-s[None, :])
    antipodal = np.argwhere(np.triu(np.minimum(sep, 1-sep) == .5, 1))
    for i, j in antipodal:
        if D[i, j] < best:
            best, triple = float(D[i, j]), (int(i), int(j), -1)
    return best, triple, len(antipodal)


def intrinsic_death(s):
    """Exact intrinsic death: min over winding triples of the largest gap, capped at 1/2."""
    gaps = np.minimum(np.abs(s[:, None]-s[None, :]), 1.)
    forward = (s[None, :]-s[:, None]) % 1.
    value, _, _ = winding_threshold(np.where(forward <= .5, forward, 1-forward), s)
    return min(value, .5), gaps


def prefix_metric(order, iu, ju, w, m, n, C):
    return C if m == len(order) else graph_metric(order, iu, ju, w, m, n)


def certified_end(order, iu, ju, w, n, close, B, C, s, W_close, W_E):
    """Largest prefix M >= close with W_{G_M} > B, or -1 if closure itself is not certified."""
    if not W_close > B:
        return -1, 0
    total = len(order)
    if W_E > B:
        return total, 0
    lo, hi, probes = close, total, 0
    while hi-lo > 1:
        mid = (lo+hi)//2
        value, _, _ = winding_threshold(prefix_metric(order, iu, ju, w, mid, n, C), s)
        probes += 1
        if value > B:
            lo = mid
        else:
            hi = mid
    return lo, probes


def one(job):
    alpha, seed, trials, first = job
    with np.load(ROOT/'results/peanut_neck/checkpoints'/f'a{alpha:g}_seed{seed}.npz', allow_pickle=False) as saved:
        s, xy, reference = saved['s'], saved['xy'], saved['reference_bars']
    n = len(s)
    C, iu, ju, w, orders, _ = build_orders(xy)
    ref = intrinsic(s)
    diameter = float(ref.max())
    b_star = float(intrinsic_death(s)[0])
    g = float(np.diff(np.r_[s, s[0]+1]).max())
    B = float(C[np.arange(n), (np.arange(n)+1) % n].max())
    W_E, _, antipodal = winding_threshold(C, s)
    rows = []
    for method in METHODS:
        t, f = trials[method], first[method]
        order = orders[method]
        conn, close = first_connected(order, iu, ju, n), closure(order, iu, ju, n)
        assert (conn, close) == (int(t['M_conn']), int(t['M_close'])) == (int(f['M_conn']), int(f['M_close']))
        selected = order[:close]
        lam = float(np.min(w[selected]/ref[iu[selected], ju[selected]]))
        D = graph_metric(order, iu, ju, w, close, n)
        W_G, triple, _ = winding_threshold(D, s)
        shf = bool(certificate(s, C, order, iu, ju)[1])
        U, probes = certified_end(order, iu, ju, w, n, close, B, C, s, W_G, W_E)
        b_saved = float(reference[0, 1])
        P = (b_saved-float(reference[0, 0]))/2
        supported = t['closure_supported'] == 'True'
        direct = max(g-B, b_saved-max(B, W_G))
        bound = min(direct, (b_saved-B)/2)/diameter
        row = dict(alpha=alpha, seed=seed, method=method, n=n, g=g, B=B, b_star=b_star, b_star_saved=b_saved,
                   diameter=diameter, P_star=P/diameter, M_conn=conn, M_close=close, pairs=len(order),
                   lambda_G=lam, edgewise_survival=bool(B < lam*b_saved), W_G=W_G,
                   W_G_triple='-'.join(map(str, triple)), W_G_over_lambda_b_star=W_G/(lam*b_star),
                   W_E=W_E, tau_E=float(t['euclidean_ring_death']), tau_G=float(t['closure_ring_death']),
                   closure_supported=supported, closure_bars=int(t['closure_bars']),
                   W_G_gt_B=bool(W_G > B), W_E_gt_B=bool(W_E > B),
                   tau_E_gt_B=bool(float(t['euclidean_ring_death']) > B+TOLERANCE),
                   antipodal_pairs=antipodal, short_hops_first=shf, certified=bool(W_G > B),
                   U_star=U, U_star_is_complete=bool(U == len(order)), U_star_probes=probes,
                   M_first=int(f['M_first']), first_at_closure=int(f['M_first']) == close,
                   first_early=int(f['M_first']) < close, theorem_applies=bool(shf and W_G > B),
                   closure_error=float(t['closure_error']), euclidean_error=float(t['euclidean_error']),
                   U_G=max(g-B, (1-lam)*b_saved)/diameter, winding_bound=bound, winding_bound_applicable=supported)
        # Later in the scan: the first prefix after closure with an edge spanning a quarter of the curve.
        late = int(t['M_shortcut'])
        lam_late = float(np.min(w[order[:late]]/ref[iu[order[:late]], ju[order[:late]]]))
        W_late, _, _ = winding_threshold(graph_metric(order, iu, ju, w, late, n), s)
        row.update(M_shortcut=late, lambda_G_shortcut=lam_late, W_G_shortcut=W_late,
                   W_G_gt_B_shortcut=bool(W_late > B), edgewise_survival_shortcut=bool(B < lam_late*b_saved))
        assert abs(lam-float(t['closure_min_metric_ratio'])) < 5e-12
        rows.append(row)
    return rows


def summarize(rows):
    groups = []
    for alpha in ALPHAS:
        for method in METHODS:
            q = [r for r in rows if r['alpha'] == alpha and r['method'] == method]
            applies = [r for r in q if r['theorem_applies']]
            groups.append(dict(alpha=alpha, method=method, samples=len(q),
                supported=sum(r['closure_supported'] for r in q),
                edgewise_survival=sum(r['edgewise_survival'] for r in q),
                W_G_gt_B=sum(r['W_G_gt_B'] for r in q),
                W_G_gt_B_shortcut=sum(r['W_G_gt_B_shortcut'] for r in q),
                edgewise_survival_shortcut=sum(r['edgewise_survival_shortcut'] for r in q), W_E_gt_B=sum(r['W_E_gt_B'] for r in q),
                tau_E_gt_B=sum(r['tau_E_gt_B'] for r in q),
                short_hops_first=sum(r['short_hops_first'] for r in q),
                theorem_applies=len(applies),
                theorem_applies_and_agrees=sum(r['first_at_closure'] for r in applies),
                certified_to_complete=sum(r['U_star_is_complete'] for r in q),
                first_at_closure=sum(r['first_at_closure'] for r in q),
                early=sum(r['first_early'] for r in q),
                early_with_short_hops_first=sum(r['first_early'] and r['short_hops_first'] for r in q),
                U_star_min=min((r['U_star'] for r in q if r['certified']), default=None),
                U_star_max=max((r['U_star'] for r in q if r['certified']), default=None),
                W_G_over_lambda_b_star_min=min(r['W_G_over_lambda_b_star'] for r in q),
                winding_bound_min=min((r['winding_bound'] for r in q if r['winding_bound_applicable']), default=None),
                winding_bound_max=max((r['winding_bound'] for r in q if r['winding_bound_applicable']), default=None)))
    early = [dict(alpha=r['alpha'], seed=r['seed'], method=r['method'], M_first=r['M_first'], M_close=r['M_close'],
                  short_hops_first=r['short_hops_first'], W_G_gt_B=r['W_G_gt_B'], U_star=r['U_star'])
             for r in rows if r['first_early']]
    checks = dict(
        W_G_at_least_lambda_b_star=all(r['W_G'] >= r['lambda_G']*r['b_star']*(1-1e-12) for r in rows),
        W_G_at_least_W_E=all(r['W_G'] >= r['W_E']*(1-1e-12) for r in rows),
        edgewise_survival_implies_W_G=all(r['W_G_gt_B'] for r in rows if r['edgewise_survival']),
        edgewise_survival_implies_W_G_shortcut=all(r['W_G_gt_B_shortcut'] for r in rows if r['edgewise_survival_shortcut']),
        W_G_shortcut_at_least_lambda_b_star=all(r['W_G_shortcut'] >= r['lambda_G_shortcut']*r['b_star']*(1-1e-12) for r in rows),
        early_stops_fail_short_hops_first=all(not r['short_hops_first'] for r in rows if r['first_early']),
        theorem_cases_agree=all(r['first_at_closure'] for r in rows if r['theorem_applies']),
        b_star_matches_saved=max(abs(r['b_star']-r['b_star_saved']) for r in rows),
        tau_G_at_least_W_G_when_certified=min((r['tau_G']-r['W_G'] for r in rows if r['W_G_gt_B'] and r['closure_bars'] == 1 and r['closure_supported']), default=None),
        tau_E_at_least_W_E_when_certified=min((r['tau_E']-r['W_E'] for r in rows if r['W_E_gt_B']), default=None),
        winding_bound_holds=all(r['closure_error'] <= r['winding_bound']+TOLERANCE/r['diameter'] for r in rows if r['winding_bound_applicable']),
        winding_bound_not_above_U_G=all(r['winding_bound'] <= r['U_G']+1e-15 for r in rows if r['winding_bound_applicable'] and r['edgewise_survival']),
        antipodal_pairs=sum(r['antipodal_pairs'] for r in rows)//3)
    return dict(status='completed', records=len(rows), clouds=len(rows)//3, groups=groups, early_stops=early, checks=checks)


def load_inputs(root):
    trials = list(csv.DictReader((root/'results/peanut_neck/trials.csv').open()))
    first = list(csv.DictReader((root/'results/first_nonempty/trials.csv').open()))
    jobs = []
    for alpha in ALPHAS:
        for seed in sorted({int(r['seed']) for r in trials}):
            pick = lambda rows: {r['method']: r for r in rows if float(r['alpha']) == alpha and int(r['seed']) == seed}
            jobs.append((alpha, seed, pick(trials), pick(first)))
    return jobs


def write(root, rows):
    folder = root/FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    with (folder/'records.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(FIELDS))
        writer.writeheader(); writer.writerows(rows)
    summary = summarize(rows)
    (folder/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    with (folder/'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary['groups'][0]))
        writer.writeheader(); writer.writerows(summary['groups'])
    return summary


def read_records(folder):
    rows = list(csv.DictReader((folder/'records.csv').open()))
    for r in rows:
        for k, v in r.items():
            if v in ('True', 'False'):
                r[k] = v == 'True'
            elif k in ('method', 'W_G_triple'):
                pass
            elif k in ('seed', 'n', 'M_conn', 'M_close', 'pairs', 'antipodal_pairs', 'U_star', 'U_star_probes',
                       'M_first', 'closure_bars', 'M_shortcut'):
                r[k] = int(v)
            else:
                r[k] = float(v)
    return rows


def verify(root, require, replay=False):
    root = Path(root)
    folder = root/FOLDER
    config = json.loads((folder/'config.json').read_text())
    expected_script = config.get('verification_script_sha256', config['script_sha256'])
    require('winding_check_script_hash', expected_script == sha(root/'code/winding_check.py'))
    require('winding_check_sources', all(config['source_sha256'][k] == sha(root/k) for k in config['source_sha256']))
    rows = read_records(folder)
    require('winding_check_scope', len(rows) == 300 and len({(r['alpha'], r['seed']) for r in rows}) == 100 and
            list(rows[0]) == list(FIELDS))
    saved = json.loads((folder/'summary.json').read_text())
    recomputed = summarize(rows)
    require('winding_check_summary_replay', json.loads(json.dumps(recomputed)) == saved)
    c = saved['checks']
    require('winding_check_implications', c['W_G_at_least_lambda_b_star'] and c['W_G_at_least_W_E'] and
            c['edgewise_survival_implies_W_G'] and c['edgewise_survival_implies_W_G_shortcut'] and
            c['W_G_shortcut_at_least_lambda_b_star'] and c['theorem_cases_agree'] and c['winding_bound_holds'] and
            c['winding_bound_not_above_U_G'] and c['b_star_matches_saved'] < 1e-7, c)
    trials = list(csv.DictReader((root/'results/peanut_neck/trials.csv').open()))
    first = list(csv.DictReader((root/'results/first_nonempty/trials.csv').open()))
    key = lambda r: (float(r['alpha']), int(r['seed']), r['method'])
    tmap, fmap = {key(r): r for r in trials}, {key(r): r for r in first}
    require('winding_check_joins', all(int(tmap[key(r)]['M_close']) == r['M_close'] and
            int(fmap[key(r)]['M_first']) == r['M_first'] and
            int(tmap[key(r)]['M_shortcut']) == r['M_shortcut'] and
            (tmap[key(r)]['short_hops_first'] == 'True') == r['short_hops_first'] for r in rows))
    require('winding_check_range_logic', all((r['U_star'] == -1) == (not r['certified']) and
            (r['U_star_is_complete'] == (r['U_star'] == r['pairs'])) and
            (not r['W_E_gt_B'] or r['U_star_is_complete'] or not r['certified']) and
            (r['U_star'] == -1 or r['M_close'] <= r['U_star'] <= r['pairs']) for r in rows))
    replayed = 0
    if replay:
        jobs = {(j[0], j[1]): j for j in load_inputs(root)}
        chosen = []
        for alpha in ALPHAS:
            group = [r for r in rows if r['alpha'] == alpha]
            pick = next((r for r in group if r['certified'] and not r['U_star_is_complete']), group[0])
            chosen.append((alpha, pick['seed']))
        for alpha, seed in chosen:
            again = one(jobs[(alpha, seed)])
            for row in again:
                record = next(r for r in rows if r['alpha'] == alpha and r['seed'] == seed and r['method'] == row['method'])
                require(f"winding_check_replay_{alpha:g}_{seed}_{row['method']}",
                        all((abs(row[k]-record[k]) <= 1e-12*max(1., abs(record[k]))) if isinstance(record[k], float) else
                            (str(row[k]) == str(record[k])) for k in FIELDS))
                replayed += 1
                # The certified end is the last prefix with W > B: check the next prefix directly.
                if row['certified'] and not row['U_star_is_complete']:
                    with np.load(root/'results/peanut_neck/checkpoints'/f'a{alpha:g}_seed{seed}.npz') as z:
                        s, xy = z['s'], z['xy']
                    C, iu, ju, w, orders, _ = build_orders(xy)
                    after, _, _ = winding_threshold(graph_metric(orders[row['method']], iu, ju, w, row['U_star']+1, len(s)), s)
                    at, _, _ = winding_threshold(graph_metric(orders[row['method']], iu, ju, w, row['U_star'], len(s)), s)
                    require(f"winding_check_range_end_{alpha:g}_{seed}_{row['method']}", at > row['B'] >= after)
    groups = saved['groups']
    return dict(status='passed', records=len(rows), replayed_records=replayed,
                W_G_gt_B=sum(g['W_G_gt_B'] for g in groups),
                W_G_gt_B_shortcut=sum(g['W_G_gt_B_shortcut'] for g in groups), edgewise_survival=sum(g['edgewise_survival'] for g in groups),
                theorem_applies=sum(g['theorem_applies'] for g in groups),
                early_stops=len(saved['early_stops']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--verify', action='store_true', help='Check saved records; generate nothing.')
    parser.add_argument('--replay', action='store_true', help='With --verify, recompute one cloud per width.')
    args = parser.parse_args()
    if args.verify:
        def require(name, condition, detail=None):
            if not bool(condition): raise AssertionError((name, detail))
        print(json.dumps(verify(ROOT, require, args.replay), indent=2)); return
    jobs = load_inputs(ROOT)
    from multiprocessing import Pool
    rows = []
    with Pool(max(1, args.workers)) as pool:
        for i, batch in enumerate(pool.imap_unordered(one, jobs), 1):
            rows.extend(batch)
            if i % 10 == 0:
                print(f'Winding check: {i}/{len(jobs)} clouds', flush=True)
    rows.sort(key=lambda r: (r['alpha'], r['seed'], METHODS.index(r['method'])))
    sources = ['results/peanut_neck/trials.csv', 'results/first_nonempty/trials.csv']
    config = dict(n=512, widths=list(ALPHAS), seeds=20, methods=list(METHODS), new_clouds=False,
                  persistence_computed=False,
                  winding_triple='three samples whose cyclic arc gaps are all below 1/2',
                  W_G='smallest graph-metric threshold containing a winding triangle or an intrinsically antipodal pair',
                  b_star='exact intrinsic death: smallest largest gap over winding triples (capped at 1/2)',
                  bound='min{max{g-B, b*-max{B,W_G}}, (b*-B)/2} / Delta*, for arc-supported closure graphs, with the saved b*',
                  later_prefix='M_shortcut from the peanut study: the first prefix after closure containing an edge whose endpoints are at least a quarter of the curve apart; W_G and lambda_G are recomputed there',
                  certified_range='largest exact prefix M >= M_close with W_{G_M} > B, found by bisection',
                  tau_E_tolerance=TOLERANCE,
                  source_sha256={k: sha(ROOT/k) for k in sources},
                  script_sha256=sha(__file__),
                  helper_sha256={k: sha(ROOT/'code'/k) for k in ('graph_methods.py',)})
    (ROOT/FOLDER).mkdir(parents=True, exist_ok=True)
    (ROOT/FOLDER/'config.json').write_text(json.dumps(config, indent=2)+'\n')
    summary = write(ROOT, rows)
    print(json.dumps({k: v for k, v in summary.items() if k != 'groups'}, indent=2))


if __name__ == '__main__':
    main()
