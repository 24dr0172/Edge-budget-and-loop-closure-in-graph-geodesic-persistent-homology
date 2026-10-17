"""Verify saved study comparisons without rerunning the large calculations."""
from pathlib import Path
import hashlib,json,sys
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from check_knn_miss_candidate import wilson
from paper1_methods import intrinsic,sample_cloud
from phase1_bridgeless import full_expected
from run_studies import one_bar_distance
from finite_diagram import matching_distance


def verify_phase1(F):
    R=F/'results/phase1';count=0
    def require(name,condition):
        nonlocal count
        if not bool(condition):raise AssertionError('Study comparisons: '+name)
        count+=1
    from verify_separation import build_report
    separation=json.loads((F/'results/separation_verification.json').read_text())
    require('exact all-order separation replay',json.loads(json.dumps(build_report()))==separation)
    require('common order impossible',separation['filtration']['common_orders']==[])
    all_study={}
    for study,expected in [('clean',120),('sweep',30000)]:
        raw=pd.read_csv(R/f'bridgeless/{study}_samples.csv')
        all_study[study]=raw
        require(study+' rows',len(raw)==expected and not raw.duplicated(['n','seed','profile','method']).any())
        source=pd.read_csv(F/('results/extra_k/samples.csv' if study=='clean' else 'results/extra_k/candidate_audit/samples.csv'))
        if study=='sweep':source=source[source.n<=32768]
        paired=raw.merge(source,on=['n','seed','profile','method'],suffixes=('','_source'),validate='one_to_one')
        require(study+' input scope',len(paired)==expected)
        require(study+' connectivity and closure replay',(paired.k_conn==paired.k_conn_source).all() and (paired.k_close==paired.k_close_source).all())
        require(study+' ordered bridge event',((raw.k_conn<=raw.k_br)&(raw.k_br<=raw.k_close)).all())
        require(study+' closure indicators',(raw.closed_at_connection==(raw.k_conn==raw.k_close)).all() and (raw.closed_at_bridgeless==(raw.k_br==raw.k_close)).all() and (raw.closed_at_plus3==(raw.k_conn+3>=raw.k_close)).all())
        summary=pd.read_csv(R/f'bridgeless/{study}_summary.csv')
        for key,g in raw.groupby(['n','profile','method']):
            q=summary[(summary.n==key[0])&(summary.profile==key[1])&(summary.method==key[2])].iloc[0]
            require('bridge group count',q.samples==len(g))
            for label,col in [('connection','closed_at_connection'),('bridgeless','closed_at_bridgeless'),('plus3','closed_at_plus3')]:
                successes=int(g[col].sum());lo,hi=wilson(successes,len(g))
                require('bridge fractions and Wilson intervals',q[label+'_closed']==successes and abs(q[label+'_rate']-successes/len(g))<1e-14 and abs(q[label+'_wilson_lo']-lo)<1e-14 and abs(q[label+'_wilson_hi']-hi)<1e-14)
    # Independent full Cartesian adjacency and low-link DFS on one original
    # sample from every profile. The computation validates all 60 originals.
    clean=all_study['clean']
    for profile,g in clean.groupby('profile'):
        seed=int(g.seed.min());s,xy=sample_cloud(seed,profile)
        sub=g[g.seed==seed];expected=full_expected(xy,int(sub.k_close.max()))
        for row in expected:
            stored=sub[sub.method==row['method']].iloc[0]
            require('independent full Cartesian/DFS bridge replay',all(stored[k]==v for k,v in row.items()))
    validation=json.loads((R/'bridgeless/validation.json').read_text())
    require('independent bridge validation scope',validation['status']=='passed' and validation['random_graph_method_checks']==360 and validation['circle_method_checks']>=170)

    baseline=pd.read_csv(R/'euclidean/samples.csv')
    require('baseline scope',len(baseline)==240 and baseline.diagram_file.nunique()==220)
    require('baseline input hashes',len(list((R/'euclidean/checkpoints').glob('*.npz')))==220)
    for i,path in enumerate(sorted((R/'euclidean/checkpoints').glob('*.npz'))):
        source=F/'results/noise/checkpoints'/path.name
        with np.load(path) as z,np.load(source) as original:
            require('baseline input identity',str(z['source_sha256'])==hashlib.sha256(source.read_bytes()).hexdigest())
            bars=z['bars'];reference=original['reference']
            require('baseline finite positive bars',np.isfinite(bars).all() and (bars[:,1]>bars[:,0]).all())
            diameter=float(intrinsic(original['latent']).max())
            error=one_bar_distance(bars,reference)/diameter
            require('baseline error',abs(float(z['e_PH'])-error)<1e-13 and np.array_equal(z['reference'],reference) and float(z['diameter'])==diameter)
            distances=cdist(original['observed'],original['observed'])
            max_error=float(np.abs(distances-intrinsic(original['latent'])).max())/diameter
            require('Euclidean stability',error<=max_error+5e-8/diameter)
            if i%20==0:
                require('independent bottleneck matching',abs(matching_distance(bars,reference)/diameter-error)<1e-12)
            for row in baseline[baseline.diagram_file==str(path.relative_to(F))].itertuples():
                require('baseline summary row',abs(row.e_PH-error)<1e-13 and row.bars==len(bars))
    for key,g in baseline.groupby(['model','delta']):
        q=pd.read_csv(R/'euclidean/summary.csv')
        q=q[(q.model==key[0])&(q.delta==key[1])].iloc[0]
        require('baseline group quantiles',q.samples==20 and all(abs(q[field]-g.e_PH.quantile(prob,interpolation='linear'))<1e-12 for field,prob in [('median_e_PH',.5),('q1_e_PH',.25),('q3_e_PH',.75)]))
    selectors=pd.read_csv(F/'results/noise/selectors.csv')
    paired=pd.read_csv(R/'euclidean/paired_selectors.csv')
    reconstructed=selectors.merge(baseline[['seed','model','delta','e_PH']],on=['seed','model','delta'],suffixes=('_selected','_euclidean'),validate='many_to_one')
    require('baseline paired rule rows',len(paired)==1440 and np.allclose(paired.e_PH_selected,reconstructed.e_PH_selected,rtol=0,atol=1e-13) and np.allclose(paired.e_PH_euclidean,reconstructed.e_PH_euclidean,rtol=0,atol=1e-13))
    for key,g in paired.groupby(['model','delta','method','selector']):
        summaries=pd.read_csv(R/'euclidean/paired_summary.csv')
        q=summaries[(summaries.model==key[0])&(summaries.delta==key[1])&(summaries.method==key[2])&(summaries.selector==key[3])].iloc[0]
        differences=g.e_PH_selected-g.e_PH_euclidean
        require('paired baseline comparisons',q.samples==len(g) and abs(q.median_paired_difference-differences.median())<1e-12 and q.selected_better==(differences<0).sum() and q.euclidean_better==(differences>0).sum() and q.equal==(differences==0).sum())

    noise=pd.read_csv(F/'results/noise/cells.csv')
    oracle=pd.read_csv(R/'grid_oracle_samples.csv')
    groups=noise[noise.grid&noise.connected].groupby(['seed','model','delta','method'])
    require('oracle scope',len(oracle)==720 and len(groups)==720)
    oracle_index=oracle.set_index(['seed','model','delta','method'])
    for key,g in groups:
        expected=g.sort_values(['e_PH','M'],kind='stable').iloc[0]
        q=oracle_index.loc[key]
        require('oracle uses connected prescribed grid only',q.oracle_M==expected.M and abs(q.oracle_e_PH-expected.e_PH)<1e-13)
    for key,g in oracle.groupby(['model','delta','method']):
        summaries=pd.read_csv(R/'grid_oracle_summary.csv')
        q=summaries[(summaries.model==key[0])&(summaries.delta==key[1])&(summaries.method==key[2])].iloc[0]
        require('oracle quantiles',q.samples==20 and all(abs(q[field]-g.oracle_e_PH.quantile(prob,interpolation='linear'))<1e-12 for field,prob in [('median_e_PH',.5),('q1_e_PH',.25),('q3_e_PH',.75)]))
    gaps=pd.read_csv(R/'selected_vs_oracle.csv')
    require('oracle rule pairing',len(gaps)==1440 and (gaps.selected_minus_oracle>=-1e-12).all() and np.allclose(gaps.selected_minus_oracle,gaps.e_PH-gaps.oracle_e_PH,rtol=0,atol=1e-12))
    amps=pd.read_csv(R/'all_amplitudes.csv')
    require('amplitude scope',len(amps)==36 and amps.delta.nunique()==6)
    for row in amps.itertuples():
        g=selectors[(selectors.model==row.model)&(selectors.delta==row.delta)&(selectors.method==row.method)]
        for rule,h in g.groupby('selector'):
            require('amplitude selected error quantiles',all(abs(getattr(row,rule+'_'+field)-h.e_PH.quantile(prob,interpolation='linear'))<1e-12 for field,prob in [('median_e_PH',.5),('q1_e_PH',.25),('q3_e_PH',.75)]))
    counts=json.loads((R/'rank_display_counts.json').read_text())
    require('display and complete record counts',[(a['seed_size_combinations'],a['profile_seed_samples'],a['method_records']) for a in counts['subsets']]==[(5000,15000,30000),(5130,15390,30780)])
    inputs=json.loads((R/'summary_inputs.json').read_text())
    for name,digest in inputs.items():require('summary input hash',hashlib.sha256((F/name).read_bytes()).hexdigest()==digest)
    literature=json.loads((R/'literature_checks.json').read_text())
    require('connectivity wording decision',not literature['connectivity_stopping']['stronger_practice_claim_supported'])
    require('spectral parameter record',literature['spectral_graph']['parameters']==dict(k=[15,100],diffusion_time=[8,64]))
    return dict(status='passed',checks_passed=count,bridge_method_records=30120,
                distinct_euclidean_diagrams=220,grid_oracle_samples=720,
                scope='Exact orientation certificate; saved bridge events and independent adjacency/DFS replays; all Euclidean input hashes, bars and bottleneck scores; grid oracles, paired comparisons and quantiles.')


if __name__=='__main__':
    F=Path(__file__).resolve().parents[1]
    report=verify_phase1(F)
    (F/'results/phase1/verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
