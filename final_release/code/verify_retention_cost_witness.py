"""Check Table 3, the additional post-closure plot and Appendix C against the saved records.

Run from anywhere; paths are relative to the release root (the parent of `code/`).
No measurements are generated. The checks cover:
  * results/retention_table/: every count, percentage and Wilson endpoint of the
    24-group table and the 48 fixed-k cells, recomputed from the joined raw rank
    (`extra_k/rank_sweep/samples.csv`) and bridge-free
    (`comparisons/bridgeless/sweep_samples.csv`) records, with all-attempt
    denominators, and every displayed cell of `retention_three_layers.csv`;
  * results/additional_figures/post_closure_error_checks.json: the plotted-series inventory, grid and
    complete-graph maxima, P* band and finite ceilings, recomputed from
    `useful_diagram/saved_prefixes.csv` and `clean_closure/clouds.csv`;
  * results/seven_point_check/seven_point_numbers.json against the release witness
    `useful_diagram/rational_counterexample.json`;
  * a rerun of code/verify_independent_checks.py (Appendix C, the geometry of
    Figures 2 and 3 and the n=512 independent-key values), compared with its
    saved output.
`python3 code/verify_retention_cost_witness.py --write` rewrites that saved output first.
"""
from pathlib import Path
import argparse, json, math, re, subprocess, sys
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
Z = 1.959963984540054
SIZES = (128, 512, 4096, 32768)
EXPECTED = {128: 1500, 512: 1500, 4096: 1500, 32768: 500}
PROFILES = ('uniform', 'strong_cosine', 'localized_trough')
METHODS = ('OR_kNN', 'mutual_kNN')
INDEPENDENT_OUTPUT = 'results/independent_checks/independent_checks_output.txt'


def wilson(s, n):
    p = s / n
    den = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / den
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / den
    return c - h, c + h


def joined_records(root):
    data = Path(root) / 'results'
    key = ['profile', 'seed', 'method', 'n']
    rank = pd.read_csv(data / 'extra_k/rank_sweep/samples.csv')
    bridge = pd.read_csv(data / 'comparisons/bridgeless/sweep_samples.csv')
    rank, bridge = rank[rank.n.isin(SIZES)], bridge[bridge.n.isin(SIZES)]
    joined = rank.merge(bridge, on=key, suffixes=('_rank', '_br'), validate='one_to_one')
    return rank, bridge, joined


def practical_counts(root, require):
    rank, bridge, J = joined_records(root)
    require('retention_table_join', len(J) == len(rank) == len(bridge) == 30000 and
            (J.k_conn_rank == J.k_conn_br).all() and (J.k_close_rank == J.k_close_br).all(),
            {'joined_records': len(J)})
    J = J.rename(columns={'k_conn_rank': 'k_conn', 'k_close_rank': 'k_close'})
    require('retention_table_rank_order', ((J.k_conn <= J.k_br) & (J.k_br <= J.k_close)).all())
    require('retention_table_gap_below_one_third', (J.largest_gap < 1 / 3).all(), float(J.largest_gap.max()))
    out = {}
    for n in SIZES:
        seeds = None
        for p in PROFILES:
            for m in METHODS:
                g = J[(J.profile == p) & (J.method == m) & (J.n == n)]
                require(f'retention_table_denominator_{p}_{m}_{n}', len(g) == EXPECTED[n] and g.seed.is_unique)
                seeds = set(g.seed) if seeds is None else seeds
                require(f'retention_table_shared_seeds_{p}_{m}_{n}', set(g.seed) == seeds)
                out[p, m, n] = dict(
                    samples=len(g),
                    k_conn=int((g.k_close == g.k_conn).sum()), k_br=int((g.k_close == g.k_br).sum()),
                    k_conn_plus3=int((g.k_close <= g.k_conn + 3).sum()),
                    k10=int((g.k_close <= 10).sum()), k15=int((g.k_close <= 15).sum()),
                    connected_k10=int((g.k_conn <= 10).sum()), connected_k15=int((g.k_conn <= 15).sum()))
    return out


def verify_retention_table(root, require):
    root = Path(root); folder = root / 'results/retention_table'
    counts = practical_counts(root, require)
    checks = json.loads((folder / 'retention_checks.json').read_text())
    require('retention_table_generator_checks', all(c['passed'] for c in checks['checks'].values()) and
            checks['joined_records'] == 30000, {'checks': len(checks['checks'])})
    saved = pd.read_csv(folder / 'retention_counts.csv')
    require('retention_table_counts_rows', len(saved) == 24)
    for r in saved.itertuples():
        c = counts[r.profile, r.method, r.n]
        for field in ('k_conn', 'k_br', 'k_conn_plus3', 'k10', 'k15'):
            s = c[field]; lo, hi = wilson(s, c['samples'])
            require(f'retention_table_cell_{r.profile}_{r.method}_{r.n}_{field}',
                    getattr(r, 'closed_at_' + field) == s and r.samples == c['samples'] and
                    abs(getattr(r, field + '_pct') - 100 * s / c['samples']) < 1e-11 and
                    abs(getattr(r, field + '_wilson_lo_pct') - 100 * lo) < 1e-9 and
                    abs(getattr(r, field + '_wilson_hi_pct') - 100 * hi) < 1e-9)
        require(f'retention_table_connected_{r.profile}_{r.method}_{r.n}',
                r.connected_at_k10 == c['connected_k10'] and r.connected_at_k15 == c['connected_k15'])
    # Saved bridge-free summary: the same counts and Wilson endpoints.
    summary = pd.read_csv(root / 'results/comparisons/bridgeless/sweep_summary.csv')
    for r in summary[summary.n.isin(SIZES)].itertuples():
        c = counts[r.profile, r.method, r.n]
        for event, field in (('connection', 'k_conn'), ('bridgeless', 'k_br'), ('plus3', 'k_conn_plus3')):
            lo, hi = wilson(c[field], c['samples'])
            require(f'retention_table_sweep_summary_{r.profile}_{r.method}_{r.n}_{event}',
                    getattr(r, event + '_closed') == c[field] and
                    abs(getattr(r, event + '_wilson_lo') - lo) < 1e-12 and abs(getattr(r, event + '_wilson_hi') - hi) < 1e-12)
    fixed = pd.read_csv(folder / 'retention_fixed_k.csv')
    require('retention_fixed_k_rows', len(fixed) == 48 and set(fixed.k) == {10, 15})
    for r in fixed.itertuples():
        c = counts[r.profile, r.method, r.n]; s = c[f'k{r.k}']; lo, hi = wilson(s, r.samples)
        require(f'retention_fixed_k_{r.profile}_{r.method}_{r.n}_{r.k}',
                r.samples == c['samples'] and r.closed == s and r.connected == c[f'connected_k{r.k}'] and
                r.closed <= r.connected <= r.samples and abs(r.closed_pct - 100 * s / r.samples) < 1e-11 and
                abs(r.wilson_lo_pct - 100 * lo) < 1e-9 and abs(r.wilson_hi_pct - 100 * hi) < 1e-9)
    # The uniform rows of the rank-sweep fixed-k summary agree.
    audit = pd.read_csv(root / 'results/extra_k/rank_sweep/fixed_k_summary.csv')
    for r in audit.itertuples():
        require(f'retention_fixed_k_uniform_summary_{r.method}_{r.k}_{r.n}', counts[r.profile, r.method, r.n][f'k{r.k}'] == r.closed)
    displayed=pd.read_csv(folder/'retention_three_layers.csv',dtype=str)
    require('retention_table_display_rows',len(displayed)==72)
    seen=set()
    for row in displayed.itertuples():
        key=(row.profile,row.method,int(row.n),row.layer)
        c=counts[key[:3]];lo,hi=wilson(c[row.layer],c['samples'])
        require('retention_display_'+str(key),key not in seen and
                (row.percent,row.wilson_lo_percent,row.wilson_hi_percent)==
                (f"{100*c[row.layer]/c['samples']:.1f}",f"{100*lo:.1f}",f"{100*hi:.1f}"))
        seen.add(key)
    require('retention_display_inventory',seen=={(*key,layer) for key in counts for layer in ('k_conn','k_br','k_conn_plus3')})
    return dict(status='passed',groups=24,joined_records=30000,fixed_k_cells=48,
                scope='All counts, Wilson intervals, fixed-k outcomes and CSV display values recomputed from joined raw records; all-attempt denominators.')


def chi(x):
    return math.sin(math.pi * x) / math.pi


def verify_cost_figure(root, require):
    root = Path(root); data = root / 'results'
    saved = json.loads((root / 'results/additional_figures/post_closure_error_checks.json').read_text())
    summary = saved['summary']
    require('post_closure_error_figure_generator_checks', all(c['passed'] for c in saved['checks'].values()), {'checks': len(saved['checks'])})
    P = pd.read_csv(data / 'useful_diagram/saved_prefixes.csv')
    complete = P.M == P.n * (P.n - 1) // 2
    series = P.groupby(['profile', 'seed', 'method'])
    require('post_closure_error_figure_inventory', len(P) == summary['post_closure_evaluations'] == 3556 and
            series.ngroups == 180 and int((P.M == P.M_close).sum()) == summary['exact_closures'] == 180 and
            int(complete.sum()) == summary['complete_graph_excluded'] == 180 and
            int((~complete).sum()) == summary['plotted_post_closure'])
    grid = P[P.is_saved_budget & ~complete]
    require('post_closure_error_figure_grid_max', float(grid.e_PH.max()) == summary['grid_max'] == float(P[~complete].e_PH.max()),
            summary['grid_max'])
    require('post_closure_error_figure_complete_max', float(P[complete].e_PH.max()) == summary['complete_graph_max'], summary['complete_graph_max'])
    clouds = pd.read_csv(data / 'clean_closure/clouds.csv')
    require('post_closure_error_figure_P_star_band', np.isclose(P.P_star.min(), summary['P_star_range'][0], rtol=0, atol=1e-12) and
            np.isclose(P.P_star.max(), summary['P_star_range'][1], rtol=0, atol=1e-12) and len(clouds) == 60,
            summary['P_star_range'])
    ceiling = np.maximum(P.true_birth - np.sin(np.pi * P.true_birth) / np.pi,
                         P.true_death - np.sin(np.pi * P.true_death) / np.pi) / P.diameter
    require('post_closure_error_figure_finite_ceiling_range', abs(ceiling.min() - summary['finite_ceiling_range'][0]) < 1e-6 and
            abs(ceiling.max() - summary['finite_ceiling_range'][1]) < 1e-6 and (P.e_PH <= ceiling + 1e-6).all(),
            [float(ceiling.min()), float(ceiling.max())])
    require('post_closure_error_figure_asymptotic_reference', abs(summary['asymptotic_reference'] - (2 / 3 - math.sqrt(3) / math.pi)) < 1e-15)
    steps = P.sort_values('M').groupby(['profile', 'seed', 'method']).e_PH.diff().dropna()
    require('post_closure_error_figure_nondecreasing', (steps >= -summary['tolerance']).all() and summary['tolerance'] == 5e-8,
            float(steps.min()))
    from pypdf import PdfReader
    require('post_closure_error_figure_pdf', len(PdfReader(root / 'results/additional_figures/post_closure_error.pdf').pages) == 1)
    return dict(status='passed', series=180, evaluations=3556, grid_max=summary['grid_max'],
                complete_graph_max=summary['complete_graph_max'], P_star_range=summary['P_star_range'],
                finite_ceiling_range=summary['finite_ceiling_range'])


def verify_seven_point(root, require):
    root = Path(root); data = root / 'results'
    numbers = json.loads((data / 'seven_point_check/seven_point_numbers.json').read_text())
    witness = json.loads((data / 'useful_diagram/rational_counterexample.json').read_text())
    require('seven_point_generator_checks', all(c['passed'] for c in numbers['checks'].values()),
            {'checks': len(numbers['checks'])})
    radius = witness['methods']['R2_global']
    require('seven_point_coordinates', [float(x) for x in witness['coordinates']] == [float(x) for x in numbers['coordinates']])
    require('seven_point_reference', abs(numbers['g'] - 0.34) < 1e-15 and
            np.allclose(numbers['reference_bar'], witness['reference_diagram'][0], rtol=0, atol=1e-12) and
            numbers['intrinsic_diameter'] == radius['before']['diameter'] == 0.499 and
            numbers['M_close_radius'] == radius['before']['M_close'] == 13)
    for stage, M in (('before', 'M13'), ('after', 'M14')):
        row = radius[stage]
        require(f'seven_point_bar_{M}', np.allclose(numbers['bars'][M], [row['birth'], row['death']], rtol=0, atol=1e-12))
        require(f'seven_point_error_{M}', abs(numbers['matchings'][M]['bottleneck_normalized'] - row['e_PH']) < 1e-12)
    require('seven_point_direct_costs',
            abs(numbers['matchings']['M13']['direct_normalized'] - 0.122371) < 5e-7 and
            abs(numbers['matchings']['M14']['direct_normalized'] - 0.156051) < 5e-7)
    for method in ('OR_kNN', 'mutual_kNN'):
        require(f'seven_point_knn_{method}', numbers['knn'][method]['M_close'] == witness['methods'][method]['before']['M_close'])
    return dict(status='passed', g=numbers['g'], errors=[radius['before']['e_PH'], radius['after']['e_PH']],
                direct_costs=[numbers['matchings']['M13']['direct_normalized'], numbers['matchings']['M14']['direct_normalized']])


def run_independent_checks(root):
    result = subprocess.run([sys.executable, '-B', str(Path(root) / 'code/verify_independent_checks.py')], cwd=root,
                            capture_output=True, text=True, check=False)
    return result.returncode, result.stdout


def verify_independent_checks_replay(root, require):
    root = Path(root)
    saved = (root / INDEPENDENT_OUTPUT).read_text()
    code, output = run_independent_checks(root)
    require('independent_checks_replay', code == 0 and output == saved and 'SUMMARY: ALL CHECKS PASSED' in saved and
            '[FAIL]' not in saved, {'ok_lines': saved.count('[ok]')})
    return dict(status='passed', checks=saved.count('[ok]'), output=INDEPENDENT_OUTPUT)


def verify_figure_sources(root, require):
    folder = Path(root) / 'results/additional_figures'
    names = ['arc_support_figure', 'ellipse_figure']
    require('additional_drawing_exports',all((folder/(n+ext)).is_file() for n in names for ext in ('.pdf','.png')))
    require('post_closure_error_figure_assets', all((folder / name).is_file() for name in ('post_closure_error.pdf', 'post_closure_error_preview.png', 'post_closure_error_grey.png', 'post_closure_error_checks.json')))
    return dict(status='passed',drawings=names,exports='PDF and 300-dpi PNG; native document sources are kept outside the release')


def verify_retention_cost_witness(root, require):
    return dict(status='passed',
                retention_table=verify_retention_table(root, require),
                cost_figure=verify_cost_figure(root, require),
                seven_point=verify_seven_point(root, require),
                independent_checks=verify_independent_checks_replay(root, require),
                figure_sources=verify_figure_sources(root, require),
                scope='Table 3, the additional post-closure plot and Appendix C checked from saved records; no persistence diagrams are computed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--write', action='store_true', help='Rewrite the saved output of verify_independent_checks.py first.')
    args = parser.parse_args()
    if args.write:
        code, output = run_independent_checks(ROOT)
        if code != 0:
            raise SystemExit('verify_independent_checks.py failed.')
        (ROOT / INDEPENDENT_OUTPUT).parent.mkdir(parents=True, exist_ok=True)
        (ROOT / INDEPENDENT_OUTPUT).write_text(output)
    checks = []
    def require(name, condition, detail=None):
        if not bool(condition):
            raise AssertionError(name)
        checks.append(name)
    report = verify_retention_cost_witness(ROOT, require)
    print(f'Table 3, the additional post-closure plot and Appendix C verified: {len(checks)} checks.')
