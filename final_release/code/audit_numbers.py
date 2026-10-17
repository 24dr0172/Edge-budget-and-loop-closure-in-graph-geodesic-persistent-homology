"""Audit manuscript counts and error extrema against saved measurements."""
from pathlib import Path
import json,sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'final_release/code'))
import numpy as np
import pandas as pd
import paper1_methods as m
PAPER=Path(__file__).resolve().parents[2]
SOURCE=PAPER/'final_release/results'
OUT=PAPER/'final_release/results'


def main():
    samples=pd.read_csv(SOURCE/'extra_k/samples.csv')
    counts=samples.groupby(['profile','method']).extra_k.agg(samples='count',closed=lambda s:int((s==0).sum())).reset_index()
    counts.to_csv(OUT/'profile_closure_counts.csv',index=False)
    cells=pd.read_csv(SOURCE/'clean_closure/cells.csv')
    own=pd.read_csv(SOURCE/'clean_closure/own_closure_cells.csv')
    joined=cells.merge(own[['seed','profile','method','budget_label','M','M_cov','covered']],on=['seed','profile','method','budget_label','M'])
    grid=joined[joined.budget_label.str.startswith('a')&joined.connected&joined.covered]
    series=grid.groupby(['seed','profile','method']).p1.agg(['min','max']);series['increase']=series['max']-series['min']
    witness=series.increase.idxmax()
    knn=pd.read_csv(SOURCE/'clean_closure/knn_event_checks.csv')
    exact=pd.concat([cells[(cells.method=='R2_global')&(cells.budget_label=='completion')],knn[knn.landmark=='own_coverage']],ignore_index=True)
    assert len(exact)==180
    exact_summary=exact.groupby(['profile','method'])[['e_inf','e_2','p1']].median().reset_index()
    exact_summary.to_csv(OUT/'exact_closure_errors.csv',index=False)
    # Replay the local count formula on exactly the original 60 sample streams.
    from check_knn_miss_candidate import neighboring_keys
    replay=[]
    for profile in m.PROFILES:
        for seed in range(20310101,20310121):
            latent,_=m.sample_cloud(seed,profile)
            keys,_,_=neighboring_keys(latent)
            for method in ['OR_kNN','mutual_kNN']:
                top=np.sort(np.partition(keys[method],-2)[-2:])
                row=samples[(samples.seed==seed)&(samples.profile==profile)&(samples.method==method)].iloc[0]
                assert list(top)==[row.k_conn,row.k_close]
                replay.append(dict(seed=seed,profile=profile,method=method,k_conn=int(top[0]),k_close=int(top[1]),closed=bool(top[0]==top[1])))
    pd.DataFrame(replay).to_csv(OUT/'original_seed_rank_replay.csv',index=False)
    large=pd.read_csv(SOURCE/'extra_k/candidate_audit/samples.csv');large=large[large.n==512]
    calibration=[]
    for method in ['OR_kNN','mutual_kNN']:
        g=large[large.method==method].copy();g['closed']=g.extra_k==0
        per_seed=g.groupby('seed').closed.sum().to_numpy()
        pmf=np.bincount(per_seed,minlength=4)/len(per_seed)
        distribution=np.array([1.])
        for _ in range(20):distribution=np.convolve(distribution,pmf)
        observed=int(((samples.method==method)&(samples.extra_k==0)).sum())
        calibration.append(dict(method=method,observed=observed,expected=float(20*per_seed.mean()),
                                probability_at_least_observed=float(distribution[observed:].sum()),
                                paired_seed_triple_pmf=pmf.tolist()))
    minimum=np.inf
    diagram_files=list((SOURCE/'useful_diagram/checkpoints').glob('*.npz'))+list((SOURCE/'matched_noise/diagrams').glob('*.npz'))
    diagram_files += [SOURCE/'matched_noise/closure_diagrams.npz',SOURCE/'ellipse_checks/diagrams.npz']
    for p in diagram_files:
            with np.load(p) as z:
                for key in z.files:
                    if '_M' not in key and not key.endswith('first_connected') and not key.startswith('oracle_') and '_coverage' not in key:continue
                    bars=z[key]
                    if bars.ndim==2 and bars.shape[1]==2 and len(bars):
                        life=bars[:,1]-bars[:,0];positive=life[life>0]
                        if len(positive):minimum=min(minimum,float(positive.min()))
    report=dict(status='passed',original_rank_replays=len(replay),
                maximum_postclosure_grid_e_PH=float(grid.p1.max()),
                maximum_postclosure_grid_increase=float(series.increase.max()),
                maximum_increase_case=dict(seed=int(witness[0]),profile=witness[1],method=witness[2],
                                           first_error=float(series.loc[witness,'min']),last_error=float(series.loc[witness,'max'])),
                paired_seed_calibration=calibration,
                shortest_original_stored_positive_lifetime=minimum,
                calibration_scope='Plug-in convolution of 20 independent seed triples from the 1500-seed empirical distribution; profiles remain paired. Descriptive, not a formal hypothesis test.')
    (OUT/'number_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
