"""Replay paired uniform/trough OR ranks and audit local affine CDF windows."""
from pathlib import Path
import argparse,csv,json
import sys
sys.dont_write_bytecode=True
import numpy as np
import paper1_methods as m
from check_knn_miss_candidate import neighboring_keys,fast_inverse

def audit(root):
    root=Path(root)
    with (root/'results/extra_k/candidate_audit/samples.csv').open() as f:rows=list(csv.DictReader(f))
    result=[]
    for n in [32768,262144,1000000]:
        get=lambda profile:{int(r['seed']):r for r in rows if int(r['n'])==n and r['profile']==profile and r['method']=='OR_kNN'}
        uniform=get('uniform');trough=get('localized_trough')
        disagreement={str(j):[] for j in [0,3]};rank_changes=[]
        counts=dict(critical_windows=0,windows_crossing_density_boundary=0,
                    closing_gaps_in_trough=0,closing_gaps_outside_trough=0,
                    samples_with_all_critical_windows_affine=0,
                    uniform_and_trough_critical_keys_equal=0)
        for t,(seed,a) in enumerate(sorted(uniform.items()),1):
            b=trough[seed];u=np.sort(np.random.default_rng(seed).uniform(size=n))
            phi=float(np.random.default_rng(seed+m.PHASE_SEED_OFFSET).uniform())
            s=fast_inverse(u,'localized_trough',phi)
            ku=neighboring_keys(u)[0]['OR_kNN'];kt=neighboring_keys(s)[0]['OR_kNN']
            for saved,keys in [(a,ku),(b,kt)]:
                top=np.partition(keys,-2)[-2:]
                assert list(map(int,top))==[int(saved['k_conn']),int(saved['k_close'])]
            for j in [0,3]:
                if (int(a['extra_k'])<=j)!=(int(b['extra_k'])<=j):
                    disagreement[str(j)].append(dict(seed=seed,uniform_extra=int(a['extra_k']),trough_extra=int(b['extra_k'])))
            if (a['k_conn'],a['k_close'])!=(b['k_conn'],b['k_close']):rank_changes.append(seed)
            threshold_u=int(np.partition(ku,-2)[-2]);threshold_t=int(np.partition(kt,-2)[-2])
            critical=np.union1d(np.flatnonzero(ku>=threshold_u),np.flatnonzero(kt>=threshold_t))
            gap=np.diff(np.r_[u,u[0]+1])[critical]
            cuts=sorted({v for interval in m.circular_arc_intervals((phi+.5)%1,.25) for v in interval if 0<v<1})
            boundaries=np.array([m.cdf_function(v,'localized_trough',phi) for v in cuts])
            # A window extends one gap behind the first endpoint and one ahead of the second.
            crosses=np.any(np.mod(boundaries[None,:]-(u[critical]-gap)[:,None],1)<=3*gap[:,None],axis=1)
            counts['critical_windows']+=len(critical)
            counts['windows_crossing_density_boundary']+=int(crosses.sum())
            counts['samples_with_all_critical_windows_affine']+=int(not crosses.any())
            counts['uniform_and_trough_critical_keys_equal']+=int(np.array_equal(ku[critical],kt[critical]))
            maximum=np.flatnonzero(kt==kt.max())
            middle=np.mod(s[maximum]+np.diff(np.r_[s,s[0]+1])[maximum]/2,1)
            inside=np.abs(np.mod(middle-(phi+.5)%1+.5,1)-.5)<.125
            counts['closing_gaps_in_trough']+=int(inside.sum());counts['closing_gaps_outside_trough']+=int((~inside).sum())
            if t%50==0:print('rank pairing',n,t,'of',len(uniform),flush=True)
        out=dict(n=n,samples=len(uniform),indicator_disagreements=disagreement,
                 different_top_rank_pairs=rank_changes,**counts)
        result.append(out);print('completed pairing n=',n,json.dumps(out),flush=True)
    return dict(status='passed',scope='Paired uniform/trough OR ranks; full raw-indicator check and seeded replay at all three largest sizes.',
        samples_replayed=sum(r['samples'] for r in result),
        mechanism='On a CDF branch with constant density, inverse-CDF mapping is affine. If the two outer gap-count windows stay on that branch, both neighbor counts, and hence the OR key, are unchanged. Density-boundary windows can change keys.',
        results=result)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--output',type=Path);a=ap.parse_args()
    out=a.output or a.release/'results/rank_pairing_replay.json'
    out.write_text(json.dumps(audit(a.release),indent=2)+'\n')
