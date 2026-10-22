"""Check the saved two-curve locality evidence behind Appendix A of the manuscript.

Reads results/noncircular_locality/ only. It recomputes every summary count from the
per-sample CSV, compares the saved summary, replays the grid thresholds nu and beta at a
coarser grid, and returns the numbers the manuscript quotes. It generates no persistence
diagrams. The thresholds are numerical grid values, not certified bounds.
"""
from pathlib import Path
import csv
import json
import sys
sys.dont_write_bytecode = True

ORDERS = ('radius', 'OR', 'mutual')
CURVES = ('ellipse', 'peanut', 'circle')


def _b(value):
    return value == 'True'


def verify_noncircular(root, require, replay_thresholds=True):
    folder = Path(root)/'results/noncircular_locality'
    with (folder/'locality_check_results.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    summary = json.loads((folder/'locality_check_summary.json').read_text())
    require('noncircular_inventory', len(rows) == 60 and
            all(sum(r['curve'] == c for r in rows) == 20 for c in CURVES) and
            len({(r['curve'], r['seed']) for r in rows}) == 60)
    r_loc = {c: summary['r_loc'][c]['r_loc'] for c in CURVES}
    for c in CURVES:
        grids = summary['r_loc'][c]['grid_values']
        require(f'noncircular_grids_agree_{c}', all(
            abs(grids['4000'][k]-grids['8000'][k]) < 1e-6 for k in ('nu', 'beta')))
    require('noncircular_circle_threshold', abs(r_loc['circle']-1/3.141592653589793) < 1e-6)
    require('noncircular_peanut_waist', abs(summary['r_loc']['peanut']['grid_values']['8000']['nu'] -
                                           summary['r_loc']['peanut']['grid_values']['8000']['beta']) < 1e-6)
    quoted = {}
    for curve in CURVES:
        group = [r for r in rows if r['curve'] == curve]
        for order in ORDERS:
            col = lambda name: [r[f'{order}_{name}'] for r in group]
            fail = [float(x) for x in col('fail_chord') if x not in ('', 'nan')]
            saved = summary[f'{curve}/{order}']
            require(f'noncircular_locality_{curve}_{order}', all(_b(x) for x in col('loc_ok')) and
                    saved['locality_hypothesis_holds'] == 20)
            require(f'noncircular_shh_support_{curve}_{order}',
                    all(_b(x) for x in col('shh_ok')) and all(_b(x) for x in col('closed_sup')) and
                    saved['shh_ok'] == saved['closed_prefix_supported'] == 20)
            # exactly one returned bar, born at B, at closure and at the end of the run
            require(f'noncircular_one_bar_{curve}_{order}',
                    all(int(x) == 1 for x in col('bars_close')) and all(int(x) == 1 for x in col('bars_end')) and
                    all(int(x) == 1 for x in col('ring_close')) and all(int(x) == 1 for x in col('ring_end')) and
                    all(float(x) == 0.0 for x in col('other_max_close')) and
                    all(float(x) == 0.0 for x in col('other_max_end')) and
                    saved['exactly_one_bar_at_close'] == saved['exactly_one_bar_at_end'] == 20 and
                    saved['max_other_bar_lifetime_at_end'] == 0.0)
            require(f'noncircular_conn_close_{curve}_{order}',
                    sum(_b(x) for x in col('conn_close')) == saved['closure_at_connectivity'])
            if curve == 'circle':
                require(f'noncircular_circle_complete_{order}', not fail and
                        saved['runs_reaching_complete_graph'] == 20 and
                        all(int(x) == 130816 for x in col('M_end')))
            else:
                require(f'noncircular_fail_chord_{curve}_{order}', len(fail) == 20 and
                        saved['runs_reaching_complete_graph'] == 0 and
                        abs(min(fail)-saved['fail_chord_min']) < 1e-15 and
                        abs(max(fail)-saved['fail_chord_max']) < 1e-15 and
                        min(fail) >= r_loc[curve]-1e-9)
                quoted.setdefault(curve, []).extend(fail)
    stats = {}
    for curve in ('ellipse', 'peanut'):
        low, high = min(quoted[curve]), max(quoted[curve])
        stats[curve] = dict(first_failure_chord_min=low, first_failure_chord_max=high,
                            r_loc=r_loc[curve], percent_above_min=100*(low/r_loc[curve]-1),
                            percent_above_max=100*(high/r_loc[curve]-1))
    if replay_thresholds:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from check_noncircular_locality import nu_beta
        for curve in ('ellipse', 'peanut'):
            nu, beta = nu_beta(curve, 2000)
            require(f'noncircular_threshold_replay_{curve}', abs(min(nu, beta)-r_loc[curve]) < 1e-5)
    return dict(status='passed', cases_non_circular=120, cases_circle_control=60,
                curves=stats, r_loc={c: r_loc[c] for c in CURVES},
                scope='Saved per-sample records recomputed; thresholds are grid values, not certified; no diagram recomputed.')


if __name__ == '__main__':
    checks = []

    def require(name, condition, detail=None):
        if not condition:
            raise AssertionError(name)
        checks.append(name)
    root = Path(__file__).resolve().parents[1]
    print(json.dumps(verify_noncircular(root, require), indent=2))
    print(f'{len(checks)} checks passed.')
