"""Check numbers reported in the manuscript against the saved records.

No measurements or diagrams are generated. Run directly for a compact report;
every verifier mode also calls verify_reported_claims().
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd


def verify_reported_claims(root, require):
    root = Path(root)
    circle = pd.read_csv(root / 'results/useful_diagram/saved_prefixes.csv')
    direct = np.maximum(circle.true_birth - circle.B,
                        circle.true_death - circle.death) / circle.diameter
    capped = np.minimum(direct, circle.P_star)
    discrepancy = float(np.max(np.abs(capped - circle.clipped_error)))
    require('reported_circle_capped_matching', len(circle) == 3556 and
            discrepancy < 1e-7, discrepancy)
    require('reported_circle_capped_column',
            np.max(np.abs(circle.clipped_error - np.minimum(circle.e_PH, circle.P_star))) < 1e-12)

    ellipse = pd.read_csv(root / 'results/ellipse/cells.csv')
    order = pd.read_csv(root / 'results/curve_order/curve_order_checks.csv')
    eo = order[order.curve == 'ellipse']
    failed = eo[eo.certified_while_alive]
    supported_failed = failed[failed.arc_supported == True]
    require('reported_ellipse_curve_order_scope', len(eo) == 3516 and
            len(failed) == 77 and len(supported_failed) == 17 and
            (supported_failed.bars == 1).all())
    # Match graph bars with independent error records by explicit sample keys.
    keys = ['seed', 'profile', 'method', 'M']
    joined = ellipse.merge(eo[keys + ['birth', 'death']], on=keys, validate='one_to_one')
    supported = joined[joined.arc_supported].copy()
    diameter = (supported.b_star - supported.g) / (2 * supported.P_star)
    direct = np.maximum(supported.g - supported.birth,
                        supported.b_star - supported.death) / diameter
    delta = float(np.max(np.abs(np.minimum(direct, supported.P_star) -
                                np.minimum(supported.e_PH, supported.P_star))))
    require('reported_ellipse_capped_matching', len(supported) > 0 and delta < 2e-7, delta)
    decreases = []
    for key, group in supported.groupby(['seed', 'profile', 'method']):
        g = group.sort_values('M')
        loss = np.minimum(g.e_PH, g.P_star).to_numpy()
        if len(loss) > 1:
            decreases.append(float(np.min(np.diff(loss))))
    require('reported_ellipse_capped_order', min(decreases) > -2e-7, min(decreases))

    cert = pd.read_csv(root / 'results/order_checks/certificates.csv')
    natural = cert[cert.method == 'R1_local']
    permuted = pd.read_csv(root / 'results/order_checks/shuffled_control.csv')
    require('reported_cknn_control', len(natural) == 60 and natural.holds.sum() == 59 and
            len(permuted) == 60 and permuted.short_hops_first_passes.sum() == 0 and
            (permuted.first_connected_H1_bars > 0).sum() == 6 and
            (permuted.M_first_connected < permuted.M_close).all())
    require('reported_cknn_precision', json.loads(
            (root / 'results/order_checks/precision_check.json').read_text())['status'] == 'passed')

    summary = pd.read_csv(root / 'results/noise/selector_summary.csv')
    selectors = pd.read_csv(root / 'results/noise/selectors.csv')
    noise = []
    for row in summary.itertuples():
        data = selectors[(selectors.model == row.model) & (selectors.delta == row.delta) &
                         (selectors.method == row.method) & (selectors.selector == row.selector)]
        require(f'reported_noise_group_{row.model}_{row.delta}_{row.method}_{row.selector}',
                len(data) == row.samples == 20 and data.seed.is_unique and
                int(data.before_closure.sum()) == row.before_closure and
                abs(float(data.e_PH.median()) - row.median_e_PH) < 1e-12)
    for delta in sorted(summary.delta.unique()):
        record = dict(delta=float(delta), evaluations_per_rule=120)
        for rule in ('first_nonempty', 'longest_lifetime'):
            data = summary[(summary.delta == delta) & (summary.selector == rule)]
            record[rule] = dict(median_range=[float(data.median_e_PH.min()),
                                             float(data.median_e_PH.max())],
                                before_closure=int(data.before_closure.sum()))
        noise.append(record)
    euclidean = pd.read_csv(root / 'results/comparisons/euclidean/summary.csv')
    baselines = euclidean[euclidean.delta == 0.01].median_e_PH
    require('reported_noise_euclidean_scale', len(baselines) == 2 and
            round(float(baselines.min()), 3) == 0.143 and round(float(baselines.max()), 3) == 0.145)

    iid = pd.read_csv(root / 'results/extra_k/rank_sweep/iid_predictions.csv')
    empirical = pd.read_csv(root / 'results/rank_table.csv')
    model = iid[(iid.method == 'mutual_kNN') & iid.n.isin([128, 512, 4096, 32768])]
    observed = empirical[(empirical.profile == 'uniform') & (empirical.method == 'mutual_kNN')]
    gap = model.merge(observed, on=['n', 'method'], validate='one_to_one')
    difference = gap.closure_percent - 100 * gap.cdf_0
    require('reported_model_difference', len(gap) == 4 and
            (difference >= 1.05).all() and (difference < 1.75).all() and
            (100 * gap.cdf_0 >= gap.wilson_lo_percent).all() and
            (100 * gap.cdf_0 <= gap.wilson_hi_percent).all(), difference.tolist())
    return dict(status='passed', circle_capped_checks=len(circle),
                supported_ellipse_capped_checks=len(supported),
                ellipse_curve_order_failures=77, supported_curve_order_failures=17,
                scope='Curve-order failures do not test alternative vertex orders.',
                noise_groups=len(summary), noise=noise)


if __name__ == '__main__':
    def require(name, condition, detail=None):
        if not condition:
            raise AssertionError((name, detail))
    print(json.dumps(verify_reported_claims(Path(__file__).resolve().parents[1], require), indent=2))
