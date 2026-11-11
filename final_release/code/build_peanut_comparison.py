"""Summarize saved peanut diagrams and reproduce the manuscript figure.

No samples or persistence diagrams are generated. Selected edges are reconstructed
to check lambda_G independently of the saved all-pairs metric ratio. The edgewise
bound U_G needs support and lambda_G b* > B. The winding bound replaces lambda_G b*
by the winding threshold W_G (results/winding_check) and needs support only.
Figure 4 uses OR exact closure and U_G. The CSV/JSON summary also includes winding-based checks as additional evidence.
"""
from pathlib import Path
import argparse, csv, hashlib, json, shutil
import numpy as np
from graph_methods import build_orders, closure, intrinsic

ROOT = Path(__file__).resolve().parents[1]
METHODS = ('R2_global', 'OR_kNN', 'mutual_kNN')


def winding_records(root):
    rows = csv.DictReader((root/'results/winding_check/records.csv').open())
    return {(float(r['alpha']), int(r['seed']), r['method']): r for r in rows}


def summarize(root):
    folder = root/'results/peanut_neck'
    rows = list(csv.DictReader((folder/'trials.csv').open()))
    winding = winding_records(root)
    detailed, checks = [], 0
    for alpha in (0., .4, .6, .8, .9):
        for seed in sorted({int(r['seed']) for r in rows}):
            with np.load(folder/'checkpoints'/f'a{alpha:g}_seed{seed}.npz') as saved:
                C, iu, ju, costs, orders, _ = build_orders(saved['xy'])
                reference = intrinsic(saved['s'])
                diameter = float(reference.max())
                bstar = float(saved['reference_bars'][0, 1])
                for method in METHODS:
                    r = next(r for r in rows if float(r['alpha']) == alpha and
                             int(r['seed']) == seed and r['method'] == method)
                    budget = closure(orders[method], iu, ju, len(C))
                    selected = orders[method][:budget]
                    lam = float(np.min(costs[selected]/reference[iu[selected], ju[selected]]))
                    assert budget == int(r['M_close'])
                    assert abs(lam-float(r['closure_min_metric_ratio'])) < 5e-12
                    bound = max(float(r['g'])-float(r['B']), (1-lam)*bstar)/diameter
                    support = r['closure_supported'] == 'True'
                    survives = float(r['B']) < lam*bstar
                    eligible = support and survives
                    bars = saved[f'{method}_closure_bars']
                    w = winding[(alpha, seed, method)]
                    W_G = float(w['W_G'])
                    winding_bound = min(max(float(r['g'])-float(r['B']), bstar-max(float(r['B']), W_G)),
                                        (bstar-float(r['B']))/2)/diameter
                    assert int(w['M_close']) == budget and abs(float(w['winding_bound'])-winding_bound) < 1e-15
                    if support:
                        assert float(r['closure_error']) <= winding_bound+5e-8/diameter
                    if eligible:
                        assert len(bars) == 1
                        assert abs(float(bars[0, 0])-float(r['B'])) < 5e-8
                        assert lam*bstar-5e-8 <= float(bars[0, 1]) <= bstar+5e-8
                        assert float(r['closure_error']) <= bound+5e-8/diameter
                    detailed.append(dict(alpha=alpha, seed=seed, method=method,
                        M_close=budget, lambda_G=lam, U_G=bound,
                        support=support, survival=survives, bound_applicable=eligible,
                        strict_bound_improvement=eligible and bound<float(r['euclidean_error']),
                        euclidean_error=float(r['euclidean_error']),
                        closure_error=float(r['closure_error']),
                        euclidean_bars=int(r['euclidean_bars']), closure_bars=len(bars),
                        W_G=W_G, winding_survival=W_G > float(r['B']), winding_bound=winding_bound,
                        winding_bound_applicable=support,
                        winding_bound_improvement=support and winding_bound<float(r['euclidean_error'])))
                    checks += 1
    groups = []
    for alpha in (0., .4, .6, .8, .9):
        for method in METHODS:
            q = [r for r in detailed if r['alpha'] == alpha and r['method'] == method]
            valid = [r['U_G'] for r in q if r['bound_applicable']]
            winding = [r['winding_bound'] for r in q if r['winding_bound_applicable']]
            groups.append(dict(alpha=alpha, method=method, samples=len(q),
                supported=sum(r['support'] for r in q),
                closure_single=sum(r['closure_bars']==1 for r in q),
                closure_multiple=sum(r['closure_bars']>1 for r in q),
                euclidean_multiple=sum(r['euclidean_bars']>1 for r in q),
                bound_applicable=len(valid),
                strict_bound_improvement=sum(r['strict_bound_improvement'] for r in q),
                actual_improvement=sum(r['closure_error']<r['euclidean_error'] for r in q),
                euclidean_error_min=min(r['euclidean_error'] for r in q),
                euclidean_error_max=max(r['euclidean_error'] for r in q),
                closure_error_min=min(r['closure_error'] for r in q),
                closure_error_max=max(r['closure_error'] for r in q),
                U_G_min=min(valid) if valid else None, U_G_max=max(valid) if valid else None,
                edgewise_survival=sum(r['survival'] for r in q),
                winding_survival=sum(r['winding_survival'] for r in q),
                winding_bound_improvement=sum(r['winding_bound_improvement'] for r in q),
                winding_bound_min=min(winding) if winding else None,
                winding_bound_max=max(winding) if winding else None))
    return dict(status='passed', scope='Frozen n=512 peanut cohort; exact closure; no new diagrams.',
                selected_edge_ratio_checks=checks, independent_seeds=20,
                methods=list(METHODS), groups=groups, records=detailed,
                trials_sha256=hashlib.sha256((folder/'trials.csv').read_bytes()).hexdigest(),
                winding_sha256=hashlib.sha256((root/'results/winding_check/records.csv').read_bytes()).hexdigest(),
                tolerance_note='Saved Ripser endpoints: 5e-8 raw units; edge ratio comparison: 5e-12.')


NAMES = {'R2_global':'Radius', 'OR_kNN':'OR', 'mutual_kNN':'Mutual'}


def table_rows(report):
    """Table 2: one row per noncircular waist; counts retain the three orders."""
    lines = []
    for alpha in (.4, .6, .8, .9):
        groups = [next(g for g in report['groups']
                       if g['alpha'] == alpha and g['method'] == method)
                  for method in METHODS]
        lo = min(g['closure_error_min'] for g in groups)
        hi = max(g['closure_error_max'] for g in groups)
        closure = f"{1000*lo:.3g}$--${1000*hi:.3g}"
        euclidean = {(g['euclidean_error_min'], g['euclidean_error_max']) for g in groups}
        assert len(euclidean) == 1
        e_lo, e_hi = euclidean.pop()
        error = f"{e_lo:.5f}$--${e_hi:.5f}" if alpha == .4 else r"P_\star"
        counts = '/'.join(str(g['supported']) for g in groups)
        single = '/'.join(str(g['closure_single']) for g in groups)
        bars = [int(r['euclidean_bars']) for r in report['records'] if r['alpha'] == alpha]
        b_lo, b_hi = min(bars), max(bars)
        bar_range = str(b_lo) if b_lo == b_hi else f"{b_lo}$--${b_hi}"
        lines.append(f"${alpha:g}$ & {counts} & {single} & ${bar_range}$ & ${error}$ & ${closure}$" + "\\\\")
    return lines


def draw(root, report, manuscript=None):
    """Draw the manuscript's peanut comparison from saved records."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
                         'legend.fontsize': 8})
    rows = list(csv.DictReader((root/'results/peanut_neck/trials.csv').open()))
    fig = plt.figure(figsize=(7.2, 5.3))
    grid = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.25])
    for i, a in enumerate((.4, .9)):
        ax = fig.add_subplot(grid[0, i])
        row = next(r for r in rows if float(r['alpha']) == a)
        theta = np.linspace(0, 2*np.pi, 1601)
        radius = 1 + a*np.cos(2*theta)
        perimeter = float(row['perimeter'])
        ax.plot(radius*np.cos(theta)/perimeter, radius*np.sin(theta)/perimeter, color='.25')
        ax.set_aspect('equal'); ax.set_xlim(-.20, .20); ax.set_ylim(-.13, .13)
        ax.set_title(f'Peanut geometry: a = {a:g}'); ax.set_xlabel('x'); ax.set_ylabel('y')
    with np.load(root/'results/peanut_neck/checkpoints/a0.4_seed20361021.npz', allow_pickle=False) as saved:
        ax = fig.add_subplot(grid[1, 0])
        for j, (key, color) in enumerate([('reference_bars', '.25'),
                ('euclidean_bars', '#b35c00'), ('OR_kNN_closure_bars', '#006da8')]):
            bars = saved[key]
            for k, (birth, death) in enumerate(bars):
                y = j + (k - (len(bars)-1)/2)*.22
                ax.plot([birth, death], [y, y], lw=2, color=color)
        ax.set_yticks([0, 1, 2], ['Arc-length metric', 'Euclidean', 'Symmetric kNN'])
        ax.set_ylim(2.5, -.5); ax.set_xlim(0, .35)
        ax.set_xlabel('Rips threshold'); ax.set_title(r'One point set at $a=0.4$')
    ax = fig.add_subplot(grid[1, 1])
    for field, label, color, marker in [('euclidean_error', 'Euclidean', '#b35c00', 'o'),
            ('closure_error', 'Symmetric kNN', '#006da8', 's'), ('U_G', 'Bound $U_G$', '0.35', '^')]:
        xx, yy, low, high = [], [], [], []
        for a in (0., .4, .6, .8, .9):
            group = [r for r in report['records'] if r['method'] == 'OR_kNN' and r['alpha'] == a]
            if field == 'U_G' and not all(r['bound_applicable'] for r in group):
                continue
            values = np.array([r[field] for r in group]); median = float(np.median(values))
            xx.append(a); yy.append(median); low.append(median-float(values.min())); high.append(float(values.max())-median)
        ax.errorbar(xx, yy, yerr=[low, high], label=label, color=color, marker=marker,
                    lw=1.1, ms=4, capsize=2)
    ax.set_yscale('log'); ax.set_ylim(3e-5, .7); ax.set_xticks([0, .4, .6, .8, .9])
    ax.set_xlabel('Peanut parameter a'); ax.set_ylabel('Normalized diagram error')
    ax.set_title('20 paired point sets')
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(0, .81))
    fig.tight_layout()
    png = root/'manuscript_figures/peanut_comparison.png'
    pdf = root/'results/figure_exports/peanut_comparison.pdf'
    png.parent.mkdir(parents=True, exist_ok=True); pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=250); fig.savefig(pdf, metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)
    if manuscript:
        shutil.copy2(png, manuscript/'figures/peanut_comparison.png')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manuscript-dir',type=Path)
    parser.add_argument('--preserve-figure',action='store_true',help='Rebuild records while retaining the supplied frozen Figure 4 files.')
    args=parser.parse_args()
    if args.preserve_figure and args.manuscript_dir:
        parser.error('--preserve-figure cannot update a manuscript figure.')
    if args.preserve_figure and not all((ROOT/p).is_file() for p in ('manuscript_figures/peanut_comparison.png','results/figure_exports/peanut_comparison.pdf')):
        parser.error('The frozen Figure 4 PNG and PDF are required with --preserve-figure.')
    report=summarize(ROOT)
    dest=ROOT/'results/peanut_neck/closure_error_bound.json'
    dest.write_text(json.dumps(report,indent=2)+'\n')
    with (dest.with_suffix('.csv')).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(report['groups'][0]));writer.writeheader();writer.writerows(report['groups'])
    if not args.preserve_figure:
        draw(ROOT,report,args.manuscript_dir)
    action="retained the supplied Figure 4" if args.preserve_figure else "redrew Figure 4"
    print(f"Checked {report['selected_edge_ratio_checks']} saved closure graphs; wrote paired bounds and {action}.")


if __name__=='__main__': main()
