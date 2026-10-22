"""Plot Figure 5: diagram error against M - M_close for all 180 clean-circle series.

Inputs (read only):
  results/useful_diagram/saved_prefixes.csv  (3,556 post-closure evaluations incl. 180 closures and 180 complete graphs)
  results/clean_closure/cells.csv            (grid + event budgets for every order; pre-closure measured points)
  results/clean_closure/knn_event_checks.csv (kNN first_connected / pre_coverage_last / own_coverage events)
  results/clean_closure/clouds.csv, closure_selection.csv, covered_monotonicity.csv, number_audit.json
Outputs: post_closure_error.pdf (vector), post_closure_error_preview.png, post_closure_error_grey.png, post_closure_error_checks.json,
Only stored evaluations are plotted; connecting segments are visual guides.
Complete-graph endpoints (M = n(n-1)/2) are verified but excluded from the panels.
"""
import json, math, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import os as _os
# Release root: $RELEASE_ROOT, else the parent of this script's folder.
_ROOT = Path(_os.environ.get('RELEASE_ROOT', Path(__file__).resolve().parent.parent))
REL = _ROOT / 'results'
OUT = _ROOT / 'results/additional_figures'
OUT.mkdir(parents=True, exist_ok=True)
TOL = 5e-8          # release tolerance: code/generate_clean.py (differences < -5e-8 count as decreases)
TOL_TRANS = 1e-7    # release tolerance for saved_transitions (experiment_core.transitions)
METHODS = ['R2_global', 'OR_kNN', 'mutual_kNN']
MTITLE = {'R2_global': 'Radius order', 'OR_kNN': r'OR-$k$NN order', 'mutual_kNN': r'Mutual-$k$NN order'}
PROFILES = ['uniform', 'strong_cosine', 'localized_trough']
PNAME = {'uniform': 'Uniform', 'strong_cosine': 'Cosine', 'localized_trough': 'Trough'}
STYLE = {'uniform': dict(ls='-', marker='o', color='0.15'),
         'strong_cosine': dict(ls='--', marker='s', color='0.40'),
         'localized_trough': dict(ls=':', marker='^', color='0.60')}
checks = {}


def check(name, ok, detail=None):
    checks[name] = dict(passed=bool(ok), detail=detail)
    if not ok:
        print('FAIL', name, detail, file=sys.stderr)


chi = lambda r: np.sin(np.pi * r) / np.pi
key = ['profile', 'seed', 'method', 'M']
P = pd.read_csv(REL / 'useful_diagram/saved_prefixes.csv')
C = pd.read_csv(REL / 'clean_closure/cells.csv')
E = pd.read_csv(REL / 'clean_closure/knn_event_checks.csv')
clouds = pd.read_csv(REL / 'clean_closure/clouds.csv')
NPAIRS = 512 * 511 // 2

check('post_closure_evaluations_3556', len(P) == 3556 and (P.M >= P.M_close).all() and not P.duplicated(key).any(), len(P))
check('series_180', P.groupby(['profile', 'seed', 'method']).ngroups == 180 and
      P.groupby(['profile', 'method']).seed.nunique().eq(20).all())
check('exact_closures_180', int((P.M == P.M_close).sum()) == 180)
check('complete_graphs_180', int((P.M == NPAIRS).sum()) == 180)
check('single_bar_each', (P.bars == 1).all())
P['bound'] = np.maximum(P.g - chi(P.g), P.true_death - chi(P.true_death)) / P.diameter
check('P_star_formula', np.allclose(P.P_star, (P.true_death - P.true_birth) / (2 * P.diameter), atol=1e-12, rtol=0))
check('g_small_gap_condition', P.g.max() < 2 / np.pi - 0.5, float(P.g.max()))
check('error_below_finite_bound_all_3556', (P.e_PH <= P.bound + 1e-6).all(),
      dict(max_ratio=float((P.e_PH / P.bound).max()), min_slack=float((P.bound - P.e_PH).min())))
check('P_star_range', abs(P.P_star.min() - 0.243836) < 5e-7 and abs(P.P_star.max() - 0.324190) < 5e-7,
      [float(P.P_star.min()), float(P.P_star.max())])

# Join saved prefixes with clean_closure records (exact budget keys, not row positions)
Cc = C[C.connected == True]
Ev = E.rename(columns={'landmark': 'budget_label'})
src = pd.concat([Cc[key + ['budget_label', 'p1', 'est_h1_bars']].assign(src='cells'),
                 Ev[key + ['budget_label', 'p1', 'est_h1_bars']].assign(src='knn_events')])
srcu = src.groupby(key).agg(p1=('p1', 'first'), p1_spread=('p1', lambda v: v.max() - v.min()),
                            est_bars=('est_h1_bars', 'max'), labels=('budget_label', lambda v: ','.join(sorted(set(v))))).reset_index()
check('duplicate_budget_records_agree', srcu.p1_spread.max() < 1e-15)
noncomp = P[P.M < NPAIRS]
J = noncomp.merge(srcu, on=key, how='left', validate='one_to_one')
check('every_noncomplete_prefix_joins_clean_record', J.p1.notna().all() and len(J) == 3376,
      dict(joined=int(J.p1.notna().sum()), noncomplete=len(noncomp)))
check('joined_error_values_agree', (J.e_PH - J.p1).abs().max() < 1e-12, float((J.e_PH - J.p1).abs().max()))
check('joined_bars_one', (J.est_bars == 1).all())
grid_post = J[J.labels.str.contains(r'(?:^|,)a', regex=True)]
gm = float(grid_post.e_PH.max())
audit = json.loads((REL / 'number_audit.json').read_text())
check('grid_max_0.043988', abs(gm - 0.043988) < 5e-7 and abs(gm - audit['maximum_postclosure_grid_e_PH']) < 1e-14, gm)
check('noncomplete_max_equals_grid_max', abs(float(noncomp.e_PH.max()) - gm) < 1e-15, float(noncomp.e_PH.max()))
cm = float(P[P.M == NPAIRS].e_PH.max())
check('complete_graph_max_0.115561', abs(cm - 0.115561) < 5e-7, cm)
comp = P[P.M == NPAIRS]
check('complete_graph_error_equals_finite_ceiling_float32', (comp.e_PH - comp.bound).abs().max() < 1e-7,
      dict(max_abs_diff=float((comp.e_PH - comp.bound).abs().max()), max_excess=float((comp.e_PH - comp.bound).max()),
           max_abs_death_minus_chi_bstar=float((comp.death - chi(comp.true_death)).abs().max())))
asym = 2 / 3 - math.sqrt(3) / math.pi
check('complete_max_exceeds_asymptotic_reference', cm > asym, dict(complete_max=cm, asymptotic=asym))

# monotonicity within the release tolerance, per series, incl. complete endpoint
mins = []
bad = 0
for _, g in P.sort_values('M').groupby(['profile', 'seed', 'method']):
    d = np.diff(g.e_PH.values)
    if len(d):
        mins.append(d.min())
        bad += int((d < -TOL).sum())
check('nondecreasing_after_closure_tol_5e-8', bad == 0, dict(decreases=bad, min_step=float(min(mins)), tolerance=TOL))
mono = pd.read_csv(REL / 'clean_closure/covered_monotonicity.csv')
check('saved_covered_monotonicity_180', len(mono) == 180 and (mono.p1_decreases == 0).all())
tr = pd.read_csv(REL / 'useful_diagram/saved_transitions.csv')
check('saved_transitions_no_drop_tol_1e-7', len(tr) == 3376 and (tr.delta >= -TOL_TRANS).all(), float(tr.delta.min()))

# Pre-closure measured evaluations: connected records with M < M_close
mc = P[P.M == P.M_close][['profile', 'seed', 'method', 'M_close', 'P_star']]
pre = src.merge(mc, on=['profile', 'seed', 'method'])
pre = pre[pre.M < pre.M_close].drop_duplicates(key)
check('pre_closure_measured_empty_and_P_star', (pre.est_h1_bars == 0).all() and (pre.p1 - pre.P_star).abs().max() < 1e-12,
      dict(points=len(pre), max_abs=float((pre.p1 - pre.P_star).abs().max())))
check('pre_closure_point_at_M_close_minus_1_every_series',
      pre[pre.M == pre.M_close - 1].groupby(['profile', 'seed', 'method']).ngroups == 180)
check('first_connected_measured_every_series',
      src[src.budget_label == 'first_connected'].groupby(['profile', 'seed', 'method']).ngroups == 180)
check('clouds_gap_max_below_0.090', clouds.g_max.max() <= 0.090, float(clouds.g_max.max()))

bounds = P.groupby(['profile', 'seed']).bound.first()
bmin, bmax = float(bounds.min()), float(bounds.max())
pmin, pmax = float(P.P_star.min()), float(P.P_star.max())

# ---------------- plot ----------------
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7, 'pdf.fonttype': 42,
                     'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm', 'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False})
fig, axes = plt.subplots(1, 3, figsize=(6.5, 3.0), sharey=True)
for ax, meth in zip(axes, METHODS):
    ax.axhspan(pmin, pmax, color='0.88', lw=0, zorder=0)
    ax.axhline(bmax, color='0.1', ls=(0, (5, 1.5, 1, 1.5)), lw=0.9, zorder=1)
    ax.axvline(0, color='0.3', lw=0.5, zorder=1)
    for prof in PROFILES:
        st = STYLE[prof]
        for seed in sorted(P[(P.profile == prof)].seed.unique()):
            post = noncomp[(noncomp.profile == prof) & (noncomp.seed == seed) & (noncomp.method == meth)].sort_values('M')
            pr = pre[(pre.profile == prof) & (pre.seed == seed) & (pre.method == meth)].sort_values('M')
            x0 = post.M_close.iloc[0]
            xs_pre = (pr.M - x0).values; ys_pre = pr.p1.values
            xs = (post.M - x0).values; ys = post.e_PH.values
            # guide segments (measured points only)
            ax.plot(np.r_[xs_pre, xs], np.r_[ys_pre, ys], ls=st['ls'], color=st['color'], lw=0.35, alpha=0.8, zorder=2)
            ax.plot(xs, ys, ls='none', marker=st['marker'], ms=1.6, mew=0, color=st['color'], zorder=3)
            ax.plot(xs_pre, ys_pre, ls='none', marker=st['marker'], ms=2.0, mfc='white', mec=st['color'], mew=0.4, zorder=3)
            ax.plot([0], [ys[0]], ls='none', marker=st['marker'], ms=2.6, mfc='black', mec='black', mew=0.3, zorder=4)
    ax.set_xscale('symlog', linthresh=10, linscale=0.8)
    ax.set_yscale('log')
    ax.set_ylim(3e-5, 0.6)
    ax.set_xlim(-6e4, 6e4)
    ax.set_xticks([-1e4, -1e2, 0, 1e2, 1e4])
    ax.set_xticklabels([r'$-10^4$', r'$-10^2$', '0', r'$10^2$', r'$10^4$'])
    ax.set_xticks([-1e3, -10, 10, 1e3], minor=True)
    ax.tick_params(axis='x', labelsize=6.5)
    ax.set_title(MTITLE[meth], fontsize=8.5)
    ax.set_xlabel(r'$M-M_{\mathrm{close}}$')
axes[0].set_ylabel(r'normalized bottleneck error $e_{\mathrm{PH}}$')
handles = [Line2D([], [], **{**STYLE[p], 'lw': 0.8, 'ms': 3}, label=PNAME[p]) for p in PROFILES]
handles += [Line2D([], [], ls='none', marker='o', mfc='white', mec='0.2', ms=3, label='measured, before closure (empty diagram)'),
            Line2D([], [], ls='none', marker='o', mfc='black', mec='black', ms=3, label='measured exact closure'),
            Patch(facecolor='0.88', label=r'range of $P_\star$ over the 60 samples (not a CI)'),
            Line2D([], [], color='0.1', ls=(0, (5, 1.5, 1, 1.5)), lw=0.9, label='largest sample-specific finite ceiling')]
fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=6.5,
           handlelength=2.2, columnspacing=1.0)
fig.tight_layout(rect=(0, 0.17, 1, 1), w_pad=0.6)
fig.savefig(OUT / 'post_closure_error.pdf', metadata={'CreationDate': None})  # no timestamp: byte-reproducible
fig.savefig(OUT / 'post_closure_error_preview.png', dpi=200)
fig.savefig(OUT / 'post_closure_error_grey.png', dpi=200)
from PIL import Image
Image.open(OUT / 'post_closure_error_grey.png').convert('L').save(OUT / 'post_closure_error_grey.png')

summary = dict(post_closure_evaluations=len(P), plotted_post_closure=len(noncomp), complete_graph_excluded=int((P.M == NPAIRS).sum()),
               exact_closures=int((P.M == P.M_close).sum()), pre_closure_measured_points=len(pre),
               grid_max=gm, complete_graph_max=cm, P_star_range=[pmin, pmax], finite_ceiling_range=[bmin, bmax],
               asymptotic_reference=asym, min_step=float(min(mins)), tolerance=TOL,
               x_range_noncomplete=[int((pre.M - pre.M_close).min()), int((noncomp.M - noncomp.M_close).max())])
(OUT / 'post_closure_error_checks.json').write_text(json.dumps(dict(checks=checks, summary=summary), indent=2, default=float) + '\n')

print(json.dumps({k: v['passed'] for k, v in checks.items()}, indent=1))
print(json.dumps(summary, indent=1, default=float))
print('all passed:', all(v['passed'] for v in checks.values()))
