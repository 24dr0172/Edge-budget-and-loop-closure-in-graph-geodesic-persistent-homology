"""Noise-free peanut curves with varying waist width, with a circle control.

The same uniform arc coordinates are used at every width within a seed.
Each prescribed order is evaluated at exact first connectivity, exact ring
closure, and the first subsequent prefix containing an edge whose endpoints
are at least a quarter of the circumference apart. Complete H1 diagrams are
computed with unrestricted Ripser, retaining every positive returned bar.
Euclidean failure is evaluated against the intrinsic reference diagram and
at the common threshold 1/4; it does not mean that every Euclidean bar is absent.
"""
from pathlib import Path
from functools import lru_cache
import argparse, csv, hashlib, json, os, sys
sys.dont_write_bytecode = True
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np
from ripser import ripser
from graph_methods import build_orders, first_connected, closure, graph_metric, intrinsic, certificate
from run_studies import arc_support
from finite_diagram import matching_distance

ROOT = Path(__file__).resolve().parents[1]
ALPHAS = (0., .4, .6, .8, .9)
METHODS = ('R2_global', 'OR_kNN', 'mutual_kNN')
BASE_SEED, N, SAMPLES, GRID, PROBE = 20361021, 512, 20, 400001, .25

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

@lru_cache(maxsize=5)
def arc_table(alpha):
    theta = np.linspace(0., 2*np.pi, GRID)
    radius = 1 + alpha*np.cos(2*theta)
    xy = np.c_[radius*np.cos(theta), radius*np.sin(theta)]
    cumulative = np.r_[0., np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
    return theta, cumulative

def sample(seed, alpha):
    s = np.sort(np.random.default_rng(seed).random(N))
    theta, cumulative = arc_table(alpha)
    length = float(cumulative[-1])
    angles = np.interp(s*length, cumulative, theta)
    radius = 1 + alpha*np.cos(2*angles)
    xy = np.c_[radius*np.cos(angles), radius*np.sin(angles)]/length
    return s, xy, length

def persistence(distance):
    result = ripser(distance, distance_matrix=True, maxdim=1, do_cocycles=True)
    bars = np.asarray(result['dgms'][1], float).reshape(-1, 2)
    keep = bars[:, 1] > bars[:, 0]
    return bars[keep], [co for co, yes in zip(result['cocycles'][1], keep) if yes]

def ring_death(bars, cocycles, birth, n):
    edges = {tuple(sorted((i, (i+1) % n))) for i in range(n)}
    deaths = [float(birth)]
    for (b, d), co in zip(bars, cocycles):
        if d > birth and sum(int(v) for i, j, v in co
                             if tuple(sorted((int(i), int(j)))) in edges) % 2:
            deaths.append(float(d))
    return max(deaths)

def alive(bars):
    return int(((bars[:, 0] <= PROBE) & (bars[:, 1] > PROBE)).sum())

def reference_matching(bars, reference, diameter):
    """Compare dropping the reference bar with matching it to each observed bar."""
    half = (bars[:, 1]-bars[:, 0])/2 if len(bars) else np.array([])
    diagonal = max(float(np.diff(reference[0])[0]/2), max(half, default=0.))
    direct = []
    for i, bar in enumerate(bars):
        direct.append(max(float(np.abs(bar-reference[0]).max()),
                          max(np.delete(half, i), default=0.)))
    best_direct = min(direct, default=float('inf'))
    return matching_distance(bars, reference)/diameter, bool(diagonal+5e-8*diameter < best_direct)

def one_seed(job):
    seed, output = job
    s = np.sort(np.random.default_rng(seed).random(N))
    reference_metric = intrinsic(s)
    reference, _ = persistence(reference_metric)
    assert len(reference) == 1
    diameter = float(reference_metric.max())
    empty_cost = float((reference[0, 1]-reference[0, 0])/(2*diameter))
    mask = ~np.eye(N, dtype=bool)
    rows = []
    for alpha in ALPHAS:
        s, xy, length = sample(seed, alpha)
        C, iu, ju, weights, orders, keys = build_orders(xy)
        B = float(C[np.arange(N), (np.arange(N)+1) % N].max())
        gap = float(np.diff(np.r_[s, s[0]+1]).max())
        euclidean, euclidean_co = persistence(C)
        e_euclidean, discard = reference_matching(euclidean, reference, diameter)
        tau_e = ring_death(euclidean, euclidean_co, B, N)
        eccentricity = C.max(axis=1)
        cone_vertex = int(eccentricity.argmin())
        payload = dict(s=s, xy=xy, reference_bars=reference, euclidean_bars=euclidean)
        for i, co in enumerate(euclidean_co):
            payload[f'euclidean_cocycle_{i}'] = co
        for method in METHODS:
            order = orders[method]
            m_conn = first_connected(order, iu, ju, N)
            m_close = closure(order, iu, ju, N)
            first_long = int(np.flatnonzero(reference_metric[iu[order], ju[order]] >= .25)[0])+1
            m_shortcut = max(m_close+1, first_long)
            edge = int(order[first_long-1])
            row = dict(alpha=alpha, waist_width=2*(1-alpha)/length, n=N, seed=seed,
                       method=method, g=gap, B=B, perimeter=length, P_star=empty_cost,
                       M_conn=m_conn, M_close=m_close, M_first_long_edge=first_long,
                       M_shortcut=m_shortcut, shortcut_i=int(iu[edge]), shortcut_j=int(ju[edge]),
                       shortcut_arc=float(reference_metric[iu[edge], ju[edge]]),
                       shortcut_chord=float(weights[edge]),
                       short_hops_first=bool(certificate(s, C, order, iu, ju)[1]),
                       intrinsic_birth=float(reference[0, 0]), intrinsic_death=float(reference[0, 1]),
                       intrinsic_alive=alive(reference), euclidean_bars=len(euclidean),
                       euclidean_last_death=float(euclidean[:, 1].max()) if len(euclidean) else 0.,
                       euclidean_ring_death=tau_e, euclidean_alive=alive(euclidean),
                       euclidean_error=e_euclidean, euclidean_reference_discarded=discard,
                       euclidean_cone_vertex=cone_vertex,
                       euclidean_cone_radius=float(eccentricity[cone_vertex]))
            for event, budget in [('connection', m_conn), ('closure', m_close), ('shortcut', m_shortcut)]:
                D = graph_metric(order, iu, ju, weights, budget, N)
                assert np.isfinite(D).all()
                bars, co = persistence(D)
                payload[f'{method}_{event}_bars'] = bars
                loss, _ = reference_matching(bars, reference, diameter)
                row.update({event+'_bars':len(bars), event+'_alive':alive(bars),
                            event+'_error':loss,
                            event+'_last_death':float(bars[:, 1].max()) if len(bars) else 0.})
                if event != 'connection':
                    supported, bad = arc_support(D, C, iu, ju, order[:budget], tol=0.)
                    row.update({event+'_supported':supported, event+'_unsupported_edges':bad,
                                event+'_ring_death':ring_death(bars, co, B, N)})
                if event == 'closure':
                    row['closure_birth'] = float(bars[:, 0].min()) if len(bars) else 0.
                    row['closure_min_metric_ratio'] = float(np.min(D[mask]/reference_metric[mask]))
            rows.append(row)
        stem = f'a{alpha:g}_seed{seed}'
        np.savez_compressed(output/'checkpoints'/f'{stem}.npz', **payload)
    return rows

def summarize(rows):
    groups = []
    for alpha in ALPHAS:
        for method in METHODS:
            group = [v for v in rows if v['alpha'] == alpha and v['method'] == method]
            rec = dict(alpha=alpha, method=method, samples=len(group),
                       waist_width=group[0]['waist_width'],
                       short_hops_first=sum(v['short_hops_first'] for v in group),
                       closure_supported=sum(v['closure_supported'] for v in group),
                       shortcut_supported=sum(v['shortcut_supported'] for v in group),
                       euclidean_reference_discarded=sum(v['euclidean_reference_discarded'] for v in group),
                       euclidean_cone_at_probe=sum(v['euclidean_cone_radius'] < PROBE for v in group))
            for event in ['connection', 'closure', 'shortcut', 'euclidean']:
                rec.update({event+'_empty':sum(v[event+'_bars'] == 0 for v in group),
                            event+'_single':sum(v[event+'_bars'] == 1 for v in group),
                            event+'_multiple':sum(v[event+'_bars'] > 1 for v in group),
                            event+'_alive_at_probe':sum(v[event+'_alive'] > 0 for v in group),
                            event+'_error_min':min(v[event+'_error'] for v in group),
                            event+'_error_max':max(v[event+'_error'] for v in group)})
            rec['closure_winding_certificate'] = sum(v['closure_min_metric_ratio'] > 3*PROBE for v in group)
            groups.append(rec)
    return dict(status='completed', independent_seeds=len(rows)//(len(ALPHAS)*len(METHODS)),
                geometric_clouds=len(rows)//len(METHODS), method_records=len(rows),
                probe=PROBE, groups=groups)

def plot(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder = root/'results/peanut_neck'
    with (folder/'trials.csv').open() as f:
        rows = list(csv.DictReader(f))
    fig = plt.figure(figsize=(7.3, 4.5))
    grid = fig.add_gridspec(2, 4, height_ratios=[.65, 1.8], hspace=.35, wspace=.4)
    for j, alpha in enumerate(ALPHAS[1:]):
        ax = fig.add_subplot(grid[0, j])
        theta = np.linspace(0, 2*np.pi, 3000)
        rad = 1+alpha*np.cos(2*theta)
        length = arc_table(alpha)[1][-1]
        ax.plot(rad*np.cos(theta)/length, rad*np.sin(theta)/length, color='0.3', lw=1)
        half = (1-alpha)/length
        ax.plot([0, 0], [-half, half], color='#b2182b', lw=2)
        ax.set_aspect('equal')
        ax.set_title(f'a = {alpha:g}, w = {2*half:.3f}', fontsize=9)
        ax.set_xlim(-.205, .205); ax.set_ylim(-.15, .15); ax.axis('off')
    ax = fig.add_subplot(grid[1, :2])
    data = np.load(folder/'checkpoints'/f'a0.8_seed{BASE_SEED}.npz', allow_pickle=False)
    labels = [('Intrinsic', 'reference_bars'), ('Connectivity', 'R2_global_connection_bars'),
              ('Closure', 'R2_global_closure_bars'), ('Euclidean', 'euclidean_bars'),
              ('Shortcut', 'R2_global_shortcut_bars')]
    colors = ['0.25', '0.5', '#0072B2', '#D55E00', '#009E73']
    for j, ((name, key), color) in enumerate(zip(labels, colors)):
        bars = data[key]
        if not len(bars):
            ax.text(.02, j, 'empty', va='center', fontsize=8, color='0.4')
        for k, (b, d) in enumerate(bars):
            y = j + (k-(len(bars)-1)/2)*.19
            ax.plot([b, d], [y, y], color=color, lw=2)
    ax.axvline(PROBE, ls='--', color='0.3', lw=.8)
    ax.set_yticks(range(len(labels)), [x[0] for x in labels], fontsize=9)
    ax.invert_yaxis(); ax.set_xlim(0, .35)
    ax.set_xlabel('Rips threshold', fontsize=9)
    ax.set_title('Full barcodes: a = 0.8, radius order', fontsize=9)
    ax.tick_params(labelsize=8)
    ax = fig.add_subplot(grid[1, 2:])
    for event, label, color, marker in [('closure','Closure','#0072B2','o'),
                                      ('euclidean','Euclidean','#D55E00','s'),
                                      ('shortcut','Shortcut','#009E73','^')]:
        xx, yy, low, high = [], [], [], []
        for alpha in ALPHAS[1:]:
            group = [v for v in rows if float(v['alpha']) == alpha and
                     (event != 'euclidean' or v['method'] == 'R2_global')]
            vals = [float(v[event+'_error'])/float(v['P_star']) for v in group]
            xx.append(float(group[0]['waist_width'])); yy.append(float(np.median(vals)))
            low.append(min(vals)); high.append(max(vals))
        order = np.argsort(xx); xx, yy, low, high = [np.asarray(v)[order] for v in [xx, yy, low, high]]
        ax.plot(xx, yy, marker=marker, color=color, label=label, lw=1.4, ms=4)
        ax.fill_between(xx, low, high, color=color, alpha=.15)
    ax.axhline(1, ls='--', color='0.4', lw=.8)
    ax.set_xlabel('Waist width w (length-one curve)', fontsize=9)
    ax.set_ylabel(r'Diagram error / $P_\star$', fontsize=9)
    ax.set_title('Medians and full observed ranges', fontsize=9)
    ax.set_ylim(-.03, 1.07); ax.tick_params(labelsize=8)
    ax.legend(fontsize=8, frameon=False, loc='center right')
    fig.subplots_adjust(left=.12, right=.985, bottom=.10, top=.92)
    dest = root/'results/additional_figures/peanut_neck_study'
    fig.savefig(dest.with_suffix('.pdf'), metadata={'CreationDate':None, 'ModDate':None})
    fig.savefig(dest.with_suffix('.png'), dpi=300)
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--samples', type=int, default=SAMPLES)
    parser.add_argument('--plot-only', action='store_true')
    args = parser.parse_args()
    if args.plot_only:
        plot(ROOT); return
    out = ROOT/'results/peanut_neck'
    (out/'checkpoints').mkdir(parents=True, exist_ok=True)
    config = dict(n=N, samples_per_width=args.samples, alphas=list(ALPHAS),
                  curve='r=1+a*cos(2*theta), normalized to length one',
                  waist_width='2*(1-a)/perimeter', base_seed=BASE_SEED,
                  profile='uniform in arc length', arc_length_grid=GRID,
                  methods=list(METHODS), probe=PROBE, noise=False,
                  coupling='Within a seed, every width and every method uses the same uniform arc coordinates; seeds are independent.',
                  events=['exact first connectivity', 'exact ring closure',
                          'first prefix after closure containing an edge with intrinsic separation at least 1/4'],
                  diagrams='Unrestricted Ripser, F2, all positive returned bars, no lifetime cutoff.',
                  source_sha256=sha(__file__),
                  helper_sha256={name:sha(ROOT/'code'/name) for name in
                                 ['graph_methods.py','run_studies.py','finite_diagram.py']},
                  example=dict(alpha=.8, seed=BASE_SEED, method='R2_global'))
    (out/'config.json').write_text(json.dumps(config, indent=2)+'\n')
    from multiprocessing import Pool
    rows = []
    with Pool(max(1, min(args.workers, args.samples))) as pool:
        for count, batch in enumerate(pool.imap_unordered(one_seed,
                     [(BASE_SEED+j, out) for j in range(args.samples)]), 1):
            rows.extend(batch)
            print(f'Peanut neck study: {count}/{args.samples} paired seeds', flush=True)
    rows.sort(key=lambda v:(v['alpha'],v['seed'],METHODS.index(v['method'])))
    with (out/'trials.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    summary = summarize(rows)
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    plot(ROOT)
    print(json.dumps(summary, indent=2), flush=True)

if __name__ == '__main__':
    main()
