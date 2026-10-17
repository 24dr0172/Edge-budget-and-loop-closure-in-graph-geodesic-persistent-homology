"""Verify saved noise and ellipse diagrams, selectors and seeded inputs."""
from pathlib import Path
import json, sys
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import paper1_methods as m
from run_studies import one_bar_distance, nn_crust

def verify_studies(F):
    R=F/'results'
    count=0
    def require(name,condition,detail=None):
        nonlocal count
        if not bool(condition): raise AssertionError(name)
        count+=1
    noise=pd.read_csv(R/'noise/cells.csv');selectors=pd.read_csv(R/'noise/selectors.csv')
    require('noise dimensions',len(selectors)==1440 and len(noise[noise.grid])==17280)
    require('all noise blocks',len(list((R/'noise/checkpoints').glob('*.csv')))==220 and len(list((R/'noise/checkpoints').glob('*.npz')))==220)
    checked=0;replayed={}
    for path in sorted((R/'noise/checkpoints').glob('*.csv')):
        d=pd.read_csv(path);row=d.iloc[0];seed=int(row.seed);model=row.model;delta=float(row.delta)
        if seed not in replayed:replayed[seed]=m.sample_cloud(seed,'uniform')
        latent,clean=replayed[seed];diam=float(m.intrinsic(latent).max())
        if model=='radial':
            normal=clean*(2*np.pi);draw=np.random.default_rng(seed+m.NOISE_SEED_OFFSET).uniform(-1,1,len(latent));displacement=draw[:,None]*normal
        else:
            rng=np.random.default_rng(seed+90000011);radius=np.sqrt(rng.uniform(size=len(latent)));angle=rng.uniform(0,2*np.pi,len(latent));displacement=radius[:,None]*np.c_[np.cos(angle),np.sin(angle)]
        with np.load(path.with_suffix('.npz')) as z:
            require('noise latent '+path.stem,np.array_equal(z['latent'],latent))
            require('noise clean '+path.stem,np.array_equal(z['clean'],clean))
            require('noise draws '+path.stem,np.array_equal(z['displacement'],displacement))
            require('noise observed '+path.stem,np.array_equal(z['observed'],clean+delta*displacement))
            reference=z['reference'];life_ref=float(np.diff(reference,axis=1)[0,0])
            for t in d.itertuples():
                bars=z[f'{t.method}_M{t.M}'];life=bars[:,1]-bars[:,0]
                require('strict positive bars',len(bars)==t.bars and (life>0).all() and np.isfinite(bars).all())
                if t.connected:
                    require('bottleneck '+path.stem,abs(one_bar_distance(bars,reference)/diam-t.e_PH)<1e-12)
                    longest=float(life.max()) if len(life) else 0.
                    require('lifetime',abs(longest-t.longest)<1e-12 and abs(longest/life_ref-t.longest_over_reference)<1e-12)
                    require('stability',t.e_PH<=t.e_inf+5e-8/diam and t.e_2<=t.e_inf+1e-12)
                else:require('disconnected losses',np.isinf(t.e_PH) and np.isinf(t.e_inf) and np.isinf(t.e_2))
                checked+=1
    # Check grid selections against the saved rows, including the stated tie rule.
    for key,g in noise[noise.grid].groupby(['seed','model','delta','method']):
        valid=g[g.connected&(g.bars>0)].sort_values('M')
        for name in ['first_nonempty','longest_lifetime']:
            q=selectors[(selectors.seed==key[0])&(selectors.model==key[1])&(selectors.delta==key[2])&(selectors.method==key[3])&(selectors.selector==name)].iloc[0]
            require('selector presence',bool(q.selected)==bool(len(valid)))
            if len(valid):
                selected=valid.iloc[0] if name=='first_nonempty' else valid[valid.longest>=valid.longest.max()-1e-7].iloc[0]
                require('selector budget',q.M==selected.M and q.before_closure==(q.M<q.M_close) and abs(q.e_PH-selected.e_PH)<1e-12)
    ellipse=pd.read_csv(R/'ellipse/cells.csv');samples=pd.read_csv(R/'ellipse/samples.csv')
    require('ellipse scope',len(samples)==60 and len(ellipse)==3516)
    require('supported ellipse single bar',(ellipse[ellipse.arc_supported].bars==1).all())
    require('local ellipse certificate',ellipse[ellipse.maximum_chord<.18].arc_supported.all() and (ellipse[ellipse.maximum_chord<.18].bars==1).all())
    with np.load(F/'results/ellipse_checks/coordinates.npz') as coordinates,np.load(F/'results/ellipse_checks/diagrams.npz') as original_bars:
        for t in samples.itertuples():
            xy=coordinates[f'{t.seed}_{t.profile}_xy'];permutation=np.random.default_rng(t.seed+110000009).permutation(len(xy))
            got={tuple(sorted((int(permutation[a]),int(permutation[b])))) for a,b in nn_crust(xy[permutation])}
            ring={tuple(sorted((a,(a+1)%len(xy)))) for a in range(len(xy))}
            require('NN-crust replay',got==ring and t.nn_crust_exact_ring)
            with np.load(R/f'ellipse/checkpoints/{t.profile}_{t.seed}.npz') as z:
                for row in ellipse[(ellipse.seed==t.seed)&(ellipse.profile==t.profile)&ellipse.closure].itertuples():
                    key=f'{row.method}_M{row.M}'
                    require('ellipse closure diagram agreement',np.allclose(z[key],original_bars[f'{t.seed}_{t.profile}_'+key],atol=1e-6,rtol=0))
    summary=pd.read_csv(R/'noise/selector_summary.csv')
    for key,g in selectors.groupby(['model','delta','method','selector']):
        q=summary[(summary.model==key[0])&(summary.delta==key[1])&(summary.method==key[2])&(summary.selector==key[3])].iloc[0]
        require('error quartiles',abs(q.q1_e_PH-g.e_PH.quantile(.25,interpolation='linear'))<1e-12 and abs(q.q3_e_PH-g.e_PH.quantile(.75,interpolation='linear'))<1e-12)
    from check_cknn_precision import check as exact_cknn_check
    precision=json.loads((R/'order_checks/precision_check.json').read_text())
    require('exact CkNN order and rigorous cost bounds',exact_cknn_check(F)==precision)
    pairing=json.loads((R/'rank_pairing_replay.json').read_text())
    rank=pd.read_csv(R/'extra_k/candidate_audit/samples.csv')
    require('paired replay size',pairing['status']=='passed' and pairing['samples_replayed']==630)
    for block in pairing['results']:
        raw=rank[(rank.n==block['n'])&(rank.method=='OR_kNN')]
        u=raw[raw.profile=='uniform'].set_index('seed');t=raw[raw.profile=='localized_trough'].set_index('seed')
        for j in [0,3]:
            differences=[dict(seed=int(seed),uniform_extra=int(u.loc[seed].extra_k),trough_extra=int(t.loc[seed].extra_k)) for seed in sorted(u.index) if (u.loc[seed].extra_k<=j)!=(t.loc[seed].extra_k<=j)]
            require('paired raw indicators',differences==block['indicator_disagreements'][str(j)])
        changed=[int(seed) for seed in sorted(u.index) if any(u.loc[seed][k]!=t.loc[seed][k] for k in ['k_conn','k_close'])]
        require('paired raw rank pairs',changed==block['different_top_rank_pairs'])
    timings=json.loads((R/'runtime_benchmark.json').read_text())
    require('timing repetitions',timings['repetitions']==3 and len(timings['measurements'])==12)
    require('timing scale',timings['largest_full_diagram_n_in_studies']==512 and timings['largest_rank_only_n']==1000000)
    for row in timings['measurements']:
        require('timing block',row['n']==512 and row['graph_evaluations']==72 and row['full_H1_diagrams']==len(row['per_diagram_calls']) and row['total_seconds']>=row['apsp_seconds']+row['diagram_seconds'])
    for path in (F/'code').glob('*.py'):compile(path.read_text(),str(path),'exec')
    separation=json.loads((R/'separation_verification.json').read_text());require('separation',separation['orientations_tested']==1024 and separation['cyclic_orientations']==0)
    cutoff=json.loads((R/'cutoff_validation.json').read_text());require('adaptive H1 validation',cutoff['status']=='passed' and cutoff['random_positive_weight_graphs']==150 and cutoff['n512_saved_unrestricted_diagrams']==126)
    return dict(status='passed',checks_passed=count,noise_stored_graphs_checked=checked,ellipse_graphs=len(ellipse),nn_crust_replays=len(samples))
