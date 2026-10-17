"""Rebuild Paper 1 tables from released rows; no graph or persistence computation."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
METHODS=('R2_global','OR_kNN','mutual_kNN')
PROFILES=('uniform','strong_cosine','localized_trough')

def reproduce_rank_table(root):
    """Export all six Table 3 groups and their pointwise Wilson intervals."""
    folder = root/'results/extra_k/candidate_audit'
    source = pd.read_csv(root/'results/phase1/bridgeless/sweep_summary.csv')
    rows = []
    for profile in PROFILES:
        for method in ('OR_kNN', 'mutual_kNN'):
            row = source[(source.profile==profile)&(source.method==method)&(source.n==512)].iloc[0]
            out = dict(profile=profile, method=method, n=512, samples=int(row.samples))
            for event in ('connection', 'bridgeless', 'plus3'):
                for field in ('rate', 'wilson_lo', 'wilson_hi'):
                    out[event+'_'+field+'_percent'] = 100*float(row[event+'_'+field])
            rows.append(out)
    pd.DataFrame(rows).to_csv(folder/'main_text_summary.csv', index=False)
    print('Table 3 export rebuilt: six rows, including bridge-free rates and Wilson intervals.')

def reproduce(root):
    data=root/'results';clean=data/'clean_closure';extra=data/'extra_k';extra.mkdir(exist_ok=True)
    cells=pd.read_csv(clean/'cells.csv')
    events=[('R2_global','first_connected'),('R2_global','pre_completion'),('R2_global','completion'),
            ('OR_kNN','first_connected'),('OR_kNN','smallest_connected_k'),
            ('OR_kNN','pre_completion'),('OR_kNN','completion'),
            ('mutual_kNN','first_connected'),('mutual_kNN','smallest_connected_k'),
            ('mutual_kNN','pre_completion'),('mutual_kNN','completion')]
    knn_events=pd.read_csv(clean/'knn_event_checks.csv')
    rows=[]
    for method,event in events:
        if method in ['OR_kNN','mutual_kNN'] and event in ['pre_completion','completion']:
            landmark={'pre_completion':'pre_coverage_last','completion':'own_coverage'}[event]
            g=knn_events[(knn_events.method==method)&(knn_events.landmark==landmark)]
        else:
            g=cells[(cells.method==method)&(cells.budget_label==event)]
        connected=g[g.connected]
        rows.append(dict(method=method,event=event,attempted=len(g),connected=len(connected),
                         empty=int((connected.est_h1_bars==0).sum()),nonempty=int((connected.est_h1_bars>0).sum())))
    pd.DataFrame(rows).to_csv(clean/'event_table.csv',index=False)
    selected=pd.read_csv(clean/'closure_selection.csv');rows=[]
    for profile in PROFILES:
        for method in METHODS:
            g=selected[(selected.profile==profile)&(selected.method==method)]
            rows.append(dict(profile=profile,method=method,samples=len(g),
                             median_e_inf=g.first_covered_e_inf.median(),median_e_2=g.first_covered_e_2.median(),
                             median_e_PH=g.first_covered_p1.median(),median_first_closed_grid_M=g.first_covered_M.median()))
    pd.DataFrame(rows).to_csv(clean/'first_closed_errors.csv',index=False)
    rows=[];grid=cells[cells.budget_label.str.startswith('a')]
    for (profile,method,M),g in grid.groupby(['profile','method','M'],sort=True):
        # Disconnected attempts stay in the median with infinite loss.
        row=dict(profile=profile,method=method,M=M,a_M=g.a_M.iloc[0],attempted=len(g),connected=int(g.connected.sum()))
        for source,name in [('e_inf','median_e_inf'),('e_2','median_e_2'),('p1','median_e_PH')]:
            loss=np.where(g.connected,g[source],np.inf)
            row[name]=float(np.median(loss))
        rows.append(row)
    pd.DataFrame(rows).to_csv(clean/'budget_summary.csv',index=False)
    ranks=pd.read_csv(clean/'rank_landmarks.csv');samples=ranks[['seed','profile','method','n','g_max','k_conn','k_largest_consecutive_key']].copy()
    samples=samples.rename(columns={'k_largest_consecutive_key':'k_close'})
    samples['extra_k']=samples.k_close-samples.k_conn;samples['miss_at_connectivity']=samples.extra_k>0
    samples.to_csv(extra/'samples.csv',index=False);rows=[]
    for (profile,method),g in samples.groupby(['profile','method'],sort=True):
        rows.append(dict(profile=profile,method=method,n=int(g.n.iloc[0]),samples=len(g),
                         misses=int(g.miss_at_connectivity.sum()),miss_rate=float(g.miss_at_connectivity.mean()),
                         median_k_conn=g.k_conn.median(),median_k_close=g.k_close.median(),
                         median_extra_k=g.extra_k.median(),q25_extra_k=g.extra_k.quantile(.25),
                         q75_extra_k=g.extra_k.quantile(.75),q90_extra_k=g.extra_k.quantile(.9),
                         min_extra_k=int(g.extra_k.min()),max_extra_k=int(g.extra_k.max())))
    pd.DataFrame(rows).to_csv(extra/'summary.csv',index=False)
    (extra/'config.json').write_text(json.dumps(dict(analysis='Paired extra integer k required for closure',
        source='results/clean_closure/rank_landmarks.csv',n=512,independent_streams=20,
        profiles=list(PROFILES),methods=['OR_kNN','mutual_kNN'],sample_method_rows=120,
        new_samples=False,new_persistence_runs=False,
        statistic='Median of paired k_close-k_conn; not difference of marginal medians.',
        miss_definition='k_conn < k_close; reference loops are positive in these clean samples.',
        scope='Descriptive n=512 evidence only; no size-scaling, asymptotic or fixed extra-k guarantee.'),indent=2)+'\n')
    closure_noise=pd.read_csv(data/'matched_noise/closure_comparison.csv');rows=[]
    for (method,event),g in closure_noise.groupby(['method','event'],sort=True):
        rows.append(dict(delta=.01,delta_n=5.12,method=method,event=event,attempted=len(g),
                         connected=int(g.connected.sum()),empty_H1=int(g.empty_H1.sum()),
                         median_longest_over_reference=g.longest_over_reference.median(),
                         long_at_half=int((g.longest_over_reference>=.5).sum())))
    pd.DataFrame(rows).to_csv(data/'matched_noise/closure_summary.csv',index=False)
    from audit_numbers import main as audit_numbers
    from summarize_studies import main as summarize_studies
    audit_numbers()
    summarize_studies(make_figures=False)
    from audit_rank_scaling import audit_rank_scaling
    audit_rank_scaling(root,write=True)
    from phase1_bridgeless import summarize as bridge_summaries
    from phase1_noise import summarize as euclidean_summaries
    from phase1_summaries import main as comparison_summaries
    from phase1_report import main as comparison_guide
    bridge_summaries('clean')
    bridge_summaries('sweep')
    euclidean_summaries()
    comparison_summaries()
    reproduce_rank_table(root)
    comparison_guide()
    print('Tables rebuilt: 11 event rows, 9 closure-error rows, 120 paired k gaps, 4 noise-closure rows.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--rank-table-only',action='store_true',help='Export all six Table 3 groups without rebuilding other outputs')
    args=ap.parse_args()
    if args.rank_table_only:
        reproduce_rank_table(args.release.resolve())
    else:
        reproduce(args.release.resolve())
