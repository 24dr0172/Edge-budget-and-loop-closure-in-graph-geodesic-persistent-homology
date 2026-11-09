"""Check the exact circle optimum on the saved clean-circle study.

The 60 uniform, cosine and trough samples at n=512 are reused with the radius,
OR and mutual orders (180 series). No new samples are drawn. For each series
the script reads the closure error E_c and the empty-diagram cost P_star and
checks the exact optimum min over connected prefixes = min{E_c, P_star}:
  * every connected pre-closure cell on the saved grid has error P_star;
  * every saved post-closure prefix has error at least min{E_c, P_star};
  * the complete graph (the chord metric, one diagram per sample) has error at
    most P_star.
It also evaluates the strict condition g < psi(b*), psi(b) = 2 sin(pi b)/pi - b,
with the exact intrinsic death b*, and the uniform threshold g < g0 ~ 0.1513.
Diagrams of the complete graph use unrestricted Ripser over F2.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, os, sys
sys.dont_write_bytecode = True
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FOLDER = 'results/exact_optimum'
METHODS = ('R2_global', 'OR_kNN', 'mutual_kNN')
TOLERANCE = 5e-8          # normalized error, as in the clean-circle study
G0 = 0.15127097021666     # root certified in results/rational_certificates/certificates.json
FIELDS = ('seed', 'profile', 'method', 'M_close', 'g', 'b_star', 'b_star_saved', 'diameter', 'P_star',
          'closure_error', 'closure_at_most_P_star', 'psi_b_star', 'strict_condition', 'g_below_g0',
          'saved_post_closure_prefixes', 'post_closure_error_max', 'post_closure_below_P_star',
          'post_closure_minimum_at_closure', 'preclosure_connected_cells', 'preclosure_error_deviation',
          'complete_graph_error', 'complete_at_most_P_star', 'optimum', 'optimum_attained_by')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def psi(b):
    return 2*math.sin(math.pi*b)/math.pi-b


def complete_errors(job):
    from ripser import ripser
    from finite_diagram import matching_distance
    from winding_check import intrinsic_death
    import graph_methods as m
    profile, seed, root = job
    with np.load(Path(root)/'results/useful_diagram/checkpoints'/f'{profile}_{seed}.npz', allow_pickle=False) as z:
        s, xy, reference = z['latent'], z['observed'], z['reference']
    chord = np.linalg.norm(xy[:, None, :]-xy[None, :, :], axis=2)
    bars = np.asarray(ripser(chord, distance_matrix=True, maxdim=1)['dgms'][1], float).reshape(-1, 2)
    bars = bars[bars[:, 1] > bars[:, 0]]
    diameter = float(m.intrinsic(s).max())
    return (profile, seed), dict(complete_bars=len(bars),
        complete_error=matching_distance(bars, reference)/diameter, b_star=float(intrinsic_death(np.asarray(s))[0]))


def build(root, complete):
    root = Path(root)
    prefixes = pd.read_csv(root/'results/useful_diagram/saved_prefixes.csv')
    cells = pd.read_csv(root/'results/clean_closure/cells.csv')
    rows = []
    for (seed, profile, method), g in prefixes.groupby(['seed', 'profile', 'method'], sort=True):
        close = g[g.M == g.M_close]
        assert len(close) == 1
        c = close.iloc[0]
        post = g[g.M > g.M_close]
        P = float(c.P_star)
        E = float(c.e_PH)
        pre = cells[(cells.seed == seed) & (cells.profile == profile) & (cells.method == method) &
                    cells.connected & (cells.M < c.M_close)]
        info = complete[(profile, int(seed))]
        b = info['b_star']
        optimum = min(E, P)
        rows.append(dict(seed=int(seed), profile=profile, method=method, M_close=int(c.M_close), g=float(c.g),
            b_star=b, b_star_saved=float(c.true_death), diameter=float(c.diameter), P_star=P,
            closure_error=E, closure_at_most_P_star=bool(E <= P), psi_b_star=psi(b),
            strict_condition=bool(c.g < psi(b)), g_below_g0=bool(c.g < G0),
            saved_post_closure_prefixes=len(post),
            post_closure_error_max=float(post.e_PH.max()) if len(post) else None,
            post_closure_below_P_star=bool((post.e_PH < P).all()),
            post_closure_minimum_at_closure=bool((post.e_PH >= E-TOLERANCE).all()),
            preclosure_connected_cells=len(pre),
            preclosure_error_deviation=float(np.abs(pre.p1-P).max()) if len(pre) else 0.,
            complete_graph_error=info['complete_error'],
            complete_at_most_P_star=bool(info['complete_error'] <= P+TOLERANCE),
            optimum=optimum, optimum_attained_by='closure' if E <= P else 'complete graph'))
    return rows


def summarize(rows):
    count = lambda key: sum(bool(r[key]) for r in rows)
    groups = []
    for method in METHODS:
        for profile in ('uniform', 'strong_cosine', 'localized_trough'):
            q = [r for r in rows if r['method'] == method and r['profile'] == profile]
            groups.append(dict(method=method, profile=profile, series=len(q),
                closure_at_most_P_star=sum(r['closure_at_most_P_star'] for r in q),
                strict_condition=sum(r['strict_condition'] for r in q),
                closure_error_max=max(r['closure_error'] for r in q),
                closure_over_P_star_max=max(r['closure_error']/r['P_star'] for r in q)))
    return dict(status='passed' if all(r['closure_at_most_P_star'] and r['complete_at_most_P_star'] and
                                       r['post_closure_minimum_at_closure'] for r in rows) else 'failed',
                series=len(rows), samples=len({(r['seed'], r['profile']) for r in rows}),
                closure_at_most_P_star=count('closure_at_most_P_star'),
                strict_condition=count('strict_condition'), g_below_g0=count('g_below_g0'),
                post_closure_below_P_star=count('post_closure_below_P_star'),
                post_closure_minimum_at_closure=count('post_closure_minimum_at_closure'),
                complete_at_most_P_star=count('complete_at_most_P_star'),
                saved_post_closure_prefixes=sum(r['saved_post_closure_prefixes'] for r in rows),
                preclosure_connected_cells=sum(r['preclosure_connected_cells'] for r in rows),
                preclosure_error_deviation_max=max(r['preclosure_error_deviation'] for r in rows),
                largest_gap=max(r['g'] for r in rows),
                smallest_strict_margin=min(r['psi_b_star']-r['g'] for r in rows),
                closure_over_P_star_max=max(r['closure_error']/r['P_star'] for r in rows),
                b_star_difference_max=max(abs(r['b_star']-r['b_star_saved']) for r in rows),
                groups=groups)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    from multiprocessing import Pool
    clouds = pd.read_csv(ROOT/'results/clean_closure/clouds.csv')
    jobs = [(p, int(s), str(ROOT)) for s, p in zip(clouds.seed, clouds.profile)]
    with Pool(max(1, args.workers)) as pool:
        complete = dict(pool.map(complete_errors, jobs))
    rows = build(ROOT, complete)
    folder = ROOT/FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    with (folder/'records.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(FIELDS)); writer.writeheader(); writer.writerows(rows)
    summary = summarize(rows)
    (folder/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    sources = ['results/useful_diagram/saved_prefixes.csv', 'results/clean_closure/cells.csv',
               'results/clean_closure/clouds.csv']
    config = dict(samples=60, orders=list(METHODS), new_samples=False, tolerance=TOLERANCE, g0=G0,
                  complete_graph='chord metric; unrestricted Ripser over F2; every bar with death > birth kept',
                  b_star='exact intrinsic death from winding triples (winding_check.intrinsic_death)',
                  source_sha256={k: sha(ROOT/k) for k in sources}, script_sha256=sha(__file__))
    (folder/'config.json').write_text(json.dumps(config, indent=2)+'\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'groups'}, indent=2))


def verify(root, require, replay=False):
    root = Path(root)
    folder = root/FOLDER
    config = json.loads((folder/'config.json').read_text())
    require('exact_optimum_script_hash', config.get('verification_script_sha256',config['script_sha256']) == sha(root/'code/exact_optimum_check.py'))
    require('exact_optimum_sources', all(config['source_sha256'][k] == sha(root/k) for k in config['source_sha256']))
    saved = pd.read_csv(folder/'records.csv', float_precision='round_trip')
    require('exact_optimum_scope', len(saved) == 180 and list(saved.columns) == list(FIELDS))
    complete = {(r.profile, int(r.seed)): dict(complete_error=float(r.complete_graph_error), b_star=float(r.b_star))
                for r in saved.itertuples()}
    rows = build(root, complete)
    for row, (_, record) in zip(rows, saved.iterrows()):
        tag = f"{row['profile']}_{row['seed']}_{row['method']}"
        require('exact_optimum_rebuild_'+tag, all(
            (pd.isna(record[k]) and row[k] is None) or
            (abs(float(record[k])-float(row[k])) <= 1e-15*max(1, abs(float(row[k]))) if isinstance(row[k], float) else str(record[k]) == str(row[k]))
            for k in FIELDS))
    summary = json.loads((folder/'summary.json').read_text())
    require('exact_optimum_summary_replay', json.loads(json.dumps(summarize(rows))) == summary)
    require('exact_optimum_outcome', summary['status'] == 'passed' and summary['closure_at_most_P_star'] == 180,
            {k: summary[k] for k in ('closure_at_most_P_star', 'strict_condition', 'complete_at_most_P_star')})
    if replay:
        for profile in ('uniform', 'strong_cosine', 'localized_trough'):
            seed = int(saved[saved.profile == profile].seed.min())
            _, again = complete_errors((profile, seed, str(root)))
            require(f'exact_optimum_complete_replay_{profile}_{seed}',
                    abs(again['complete_error']-complete[(profile, seed)]['complete_error']) < 1e-15 and
                    again['b_star'] == complete[(profile, seed)]['b_star'])
    return dict(status='passed', series=180, closure_at_most_P_star=summary['closure_at_most_P_star'],
                strict_condition=summary['strict_condition'], complete_at_most_P_star=summary['complete_at_most_P_star'])


if __name__ == '__main__':
    if '--verify' in sys.argv:
        def require(name, condition, detail=None):
            if not bool(condition): raise AssertionError((name, detail))
        print(json.dumps(verify(ROOT, require, '--replay' in sys.argv), indent=2))
    else:
        main()
