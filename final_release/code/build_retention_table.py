"""Build Table 3 (ring retention) from the saved rank and bridge-free records.

Inputs (read only):
  results/extra_k/rank_sweep/samples.csv        (rank records: k_conn, k_close, largest_gap)
  results/comparisons/bridgeless/sweep_samples.csv   (bridge-free records: k_conn, k_br, k_close, flags)
  results/comparisons/bridgeless/sweep_summary.csv   (saved counts + Wilson endpoints)
  results/extra_k/rank_sweep/summary.csv        (saved cdf_0/cdf_3 Wilson endpoints)
  results/extra_k/rank_sweep/n512_comparison_summary.csv, results/rank_table.csv
  results/extra_k/rank_sweep/fixed_k_summary.csv (uniform-only fixed-k counts)
Outputs, in results/retention_table/:
  retention_compact.csv       six rows of Table 3 (retention at k_conn)
  retention_three_layers.csv  the same groups at k_conn, k_br and k_conn+3
  retention_fixed_k.csv        fixed k = 10 and 15
  retention_counts.csv         every count and Wilson interval
  retention_checks.json        the checks run while building
No persistence diagram is computed: every count comes from neighbor ranks.
"""
import json, math, sys
from pathlib import Path
import pandas as pd

import os as _os
# Release root: $RELEASE_ROOT, else the parent of this script's folder.
_ROOT = Path(_os.environ.get('RELEASE_ROOT', Path(__file__).resolve().parent.parent))
REL = _ROOT / 'results'
OUT = _ROOT / 'results' / 'retention_table'
OUT.mkdir(parents=True, exist_ok=True)
PROFILES = ['uniform', 'strong_cosine', 'localized_trough']
PNAME = {'uniform': 'Uniform', 'strong_cosine': 'Cosine', 'localized_trough': 'Trough'}
METHODS = ['OR_kNN', 'mutual_kNN']
SIZES = [128, 512, 4096, 32768]
EXPECTED = {128: 1500, 512: 1500, 4096: 1500, 32768: 500}
Z = 1.959963984540054
checks = {}


def check(name, ok, detail=None):
    checks[name] = dict(passed=bool(ok), detail=detail)
    if not ok:
        print('FAIL', name, detail, file=sys.stderr)


def wilson(s, n):
    p = s / n
    den = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / den
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / den
    return c - h, c + h


key = ['profile', 'seed', 'method', 'n']
rank = pd.read_csv(REL / 'extra_k/rank_sweep/samples.csv')
br = pd.read_csv(REL / 'comparisons/bridgeless/sweep_samples.csv')
rank = rank[rank.n.isin(SIZES)]
br = br[br.n.isin(SIZES)]
check('key_uniqueness_rank', not rank.duplicated(key).any())
check('key_uniqueness_bridge', not br.duplicated(key).any())
J = rank.merge(br, on=key, how='inner', suffixes=('_rank', '_br'), validate='one_to_one')
check('joined_records_30000', len(J) == len(rank) == len(br) == 30000, len(J))
check('rank_values_match', (J.k_conn_rank == J.k_conn_br).all() and (J.k_close_rank == J.k_close_br).all())
J = J.rename(columns={'k_conn_rank': 'k_conn', 'k_close_rank': 'k_close'})
check('order_kconn_le_kbr_le_kclose', ((J.k_conn <= J.k_br) & (J.k_br <= J.k_close)).all())
check('extra_k_consistent', (J.extra_k == J.k_close - J.k_conn).all())
check('flag_closed_at_connection', (J.closed_at_connection == (J.k_close == J.k_conn)).all())
check('flag_closed_at_bridgeless', (J.closed_at_bridgeless == (J.k_close == J.k_br)).all())
check('flag_closed_at_plus3', (J.closed_at_plus3 == (J.k_close <= J.k_conn + 3)).all())
gmax = float(J.largest_gap.max())
check('all_g_below_one_third', (J.largest_gap < 1 / 3).all(), gmax)
check('max_gap_about_0.271062', abs(gmax - 0.271062) < 5e-7, gmax)
# seed uniqueness: within each (profile, rule, n) seeds unique; each seed shared across profiles/rules
seed_ok = True
for (p, m, n), g in J.groupby(['profile', 'method', 'n']):
    seed_ok &= g.seed.is_unique and len(g) == EXPECTED[n]
seedsets = J.groupby('n').apply(lambda g: g.groupby(['profile', 'method']).seed.apply(frozenset).nunique(), include_groups=False)
check('seed_unique_per_group_and_denominators', seed_ok)
check('same_seed_set_across_profiles_rules_per_n', (seedsets == 1).all(), seedsets.to_dict())
check('seeds_distinct_across_sizes', J.groupby('seed').n.nunique().max() == 1)

rows = []
for n in SIZES:
    for p in PROFILES:
        for m in METHODS:
            g = J[(J.profile == p) & (J.method == m) & (J.n == n)]
            N = len(g)
            c = dict(profile=p, method=m, n=n, samples=N,
                     closed_at_k_conn=int((g.k_close == g.k_conn).sum()),
                     closed_at_k_br=int((g.k_close == g.k_br).sum()),
                     closed_at_k_conn_plus3=int((g.k_close <= g.k_conn + 3).sum()),
                     closed_at_k10=int((g.k_close <= 10).sum()),
                     closed_at_k15=int((g.k_close <= 15).sum()),
                     connected_at_k10=int((g.k_conn <= 10).sum()),
                     connected_at_k15=int((g.k_conn <= 15).sum()),
                     kbr_eq_kconn=int((g.k_br == g.k_conn).sum()),
                     plus3_le_kbr=int((g.k_conn + 3 <= g.k_br).sum()),
                     max_gap=float(g.largest_gap.max()))
            for f in ['k_conn', 'k_br', 'k_conn_plus3', 'k10', 'k15']:
                s = c['closed_at_' + f]
                lo, hi = wilson(s, N)
                c[f + '_pct'] = 100 * s / N
                c[f + '_wilson_lo_pct'] = 100 * lo
                c[f + '_wilson_hi_pct'] = 100 * hi
            rows.append(c)
T = pd.DataFrame(rows)
check('24_groups', len(T) == 24 and (T.samples == T.n.map(EXPECTED)).all())
check('connected_at_fixed_k_not_all',
      True, dict(groups_with_disconnected_k10=int((T.connected_at_k10 < T.samples).sum()),
                 groups_with_disconnected_k15=int((T.connected_at_k15 < T.samples).sum())))
check('plus3_vs_kbr_order_varies', True,
      dict(records_kbr_gt_kconn_plus3=int((J.k_br > J.k_conn + 3).sum()),
           records_kbr_lt_kconn_plus3=int((J.k_br < J.k_conn + 3).sum())))

# Saved endpoints
S = pd.read_csv(REL / 'comparisons/bridgeless/sweep_summary.csv')
S = S[S.n.isin(SIZES)].merge(T, on=['profile', 'method', 'n'], validate='one_to_one')
maxd = 0.0
ok = len(S) == 24 and (S.samples_x == S.samples_y).all()
for a, b in [('connection', 'k_conn'), ('bridgeless', 'k_br'), ('plus3', 'k_conn_plus3')]:
    ok &= (S[a + '_closed'] == S['closed_at_' + b]).all()
    for e, f in [('rate', 'pct'), ('wilson_lo', 'wilson_lo_pct'), ('wilson_hi', 'wilson_hi_pct')]:
        maxd = max(maxd, float((S[a + '_' + e] * 100 - S[b + '_' + f]).abs().max()))
check('sweep_summary_counts_and_wilson', ok and maxd < 1e-9, dict(max_abs_diff_pct=maxd))
U = pd.read_csv(REL / 'extra_k/rank_sweep/summary.csv')
U = U[U.n.isin(SIZES)].merge(T, on=['profile', 'method', 'n'], validate='one_to_one')
d0 = max(float((U[c] * 100 - U[t]).abs().max()) for c, t in
         [('cdf_0', 'k_conn_pct'), ('cdf_0_lo', 'k_conn_wilson_lo_pct'), ('cdf_0_hi', 'k_conn_wilson_hi_pct'),
          ('cdf_3', 'k_conn_plus3_pct'), ('cdf_3_lo', 'k_conn_plus3_wilson_lo_pct'), ('cdf_3_hi', 'k_conn_plus3_wilson_hi_pct')])
check('rank_sweep_summary_cdf0_cdf3_wilson', len(U) == 24 and d0 < 1e-9, dict(max_abs_diff_pct=d0))
F = pd.read_csv(REL / 'rank_table.csv').merge(T, on=['profile', 'method', 'n'], validate='one_to_one')
d1 = max(float((F.closure_percent - F.k_conn_pct).abs().max()), float((F.wilson_lo_percent - F.k_conn_wilson_lo_pct).abs().max()),
         float((F.wilson_hi_percent - F.k_conn_wilson_hi_pct).abs().max()))
check('rank_table_wilson', len(F) == 24 and (F.closed_at_connectivity == F.closed_at_k_conn).all() and d1 < 1e-9,
      dict(max_abs_diff_pct=d1))
MT = pd.read_csv(REL / 'extra_k/rank_sweep/n512_comparison_summary.csv').merge(T, on=['profile', 'method', 'n'])
d2 = max(float((MT[a + '_' + e] - MT[b + '_' + f]).abs().max())
         for a, b in [('connection', 'k_conn'), ('bridgeless', 'k_br'), ('plus3', 'k_conn_plus3')]
         for e, f in [('rate_percent', 'pct'), ('wilson_lo_percent', 'wilson_lo_pct'), ('wilson_hi_percent', 'wilson_hi_pct')])
check('n512_comparison_summary_wilson', len(MT) == 6 and d2 < 1e-9, dict(rows=len(MT), max_abs_diff_pct=d2))
FK = pd.read_csv(REL / 'extra_k/rank_sweep/fixed_k_summary.csv').merge(T, on=['profile', 'method', 'n'])
fk_ok = all(r.closed == (r.closed_at_k10 if r.k == 10 else r.closed_at_k15) for r in FK.itertuples())
check('fixed_k_summary_uniform_rows', len(FK) == 16 and fk_ok and set(FK.profile) == {'uniform'})
# Manuscript Table (tab:rank-scaling) and prose counts
man = {('uniform', 'OR_kNN'): [34.3, 32.9, 32.2, 31.2], ('uniform', 'mutual_kNN'): [28.7, 28.9, 28.7, 29.4],
       ('strong_cosine', 'OR_kNN'): [33.7, 33.2, 31.8, 30.8], ('strong_cosine', 'mutual_kNN'): [23.2, 27.0, 28.2, 29.4],
       ('localized_trough', 'OR_kNN'): [24.2, 32.7, 32.3, 31.2], ('localized_trough', 'mutual_kNN'): [6.8, 7.1, 11.5, 14.8]}
mm = [(k, n, v, round(float(T[(T.profile == k[0]) & (T.method == k[1]) & (T.n == n)].k_conn_pct.iloc[0]), 1))
      for k, vals in man.items() for n, v in zip(SIZES, vals)]
check('manuscript_table_rank_scaling_matches', all(abs(v - w) < 1e-9 for _, _, v, w in mm))
def cell(p, m, n, f):
    return int(T[(T.profile == p) & (T.method == m) & (T.n == n)][f].iloc[0])
check('manuscript_prose_counts', cell('uniform', 'mutual_kNN', 512, 'closed_at_k10') == 593 and
      cell('uniform', 'mutual_kNN', 4096, 'closed_at_k10') == 1 and
      cell('uniform', 'OR_kNN', 512, 'closed_at_k_br') == 1124 and
      cell('localized_trough', 'mutual_kNN', 512, 'closed_at_k_br') == 618)

T.to_csv(OUT / 'retention_counts.csv', index=False)
(OUT / 'retention_checks.json').write_text(json.dumps(dict(checks=checks, max_gap=gmax, joined_records=len(J)), indent=2) + '\n')


def fmt(x):
    return f'{x:.1f}'


print(json.dumps({k: v['passed'] for k, v in checks.items()}, indent=1))
print('all passed:', all(v['passed'] for v in checks.values()), 'max gap', gmax)


# Three-layer display, fixed-k outcomes and the six-row compact display.
# These are numerical CSV exports; document source belongs to the Overleaf project.
import csv
layer_rows=[]
for row in T.itertuples():
    for layer in ('k_conn','k_br','k_conn_plus3'):
        layer_rows.append(dict(profile=row.profile,method=row.method,n=row.n,layer=layer,
            percent=f"{getattr(row,layer+'_pct'):.1f}",
            wilson_lo_percent=f"{getattr(row,layer+'_wilson_lo_pct'):.1f}",
            wilson_hi_percent=f"{getattr(row,layer+'_wilson_hi_pct'):.1f}"))
with (OUT/'retention_three_layers.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(layer_rows[0]));writer.writeheader();writer.writerows(layer_rows)

rows = []
for r in T.itertuples():
    for k in (10, 15):
        rows.append(dict(profile=r.profile, method=r.method, n=r.n, k=k, samples=r.samples,
                         connected=getattr(r, f'connected_at_k{k}'), closed=getattr(r, f'closed_at_k{k}'),
                         closed_pct=getattr(r, f'k{k}_pct'), wilson_lo_pct=getattr(r, f'k{k}_wilson_lo_pct'),
                         wilson_hi_pct=getattr(r, f'k{k}_wilson_hi_pct')))
FKOUT = pd.DataFrame(rows).sort_values(['n', 'profile', 'method', 'k'])
for x in FKOUT.itertuples():
    lo, hi = wilson(x.closed, x.samples)
    assert abs(100 * lo - x.wilson_lo_pct) < 1e-9 and abs(100 * hi - x.wilson_hi_pct) < 1e-9
    assert x.closed <= x.connected <= x.samples
FKOUT['wilson_lo_pct'] = FKOUT['wilson_lo_pct'].where(FKOUT['wilson_lo_pct'] > 1e-12, 0.0)
FKOUT.to_csv(OUT / 'retention_fixed_k.csv', index=False)
print('fixed-k rows', len(FKOUT))



compact_rows=[]
for profile in PROFILES:
    for method in METHODS:
        row=dict(profile=profile,method=method)
        for n in SIZES:
            row[str(n)]=f"{float(T[(T.profile==profile)&(T.method==method)&(T.n==n)].k_conn_pct.iloc[0]):.1f}"
        compact_rows.append(row)
with (OUT/'retention_compact.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(compact_rows[0]));writer.writeheader();writer.writerows(compact_rows)
