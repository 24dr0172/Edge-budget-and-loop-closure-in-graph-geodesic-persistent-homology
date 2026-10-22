"""Verify saved study comparisons without rerunning the large calculations."""
from pathlib import Path
import hashlib,json,sys
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from rank_sweep import wilson
from graph_methods import intrinsic,sample_cloud
from bridge_free import full_expected
from run_studies import one_bar_distance
from finite_diagram import matching_distance


def verify_comparisons(F):
    R=F/'results/comparisons';count=0
    def require(name,condition):
        nonlocal count
        if not bool(condition):raise AssertionError('Study comparisons: '+name)
        count+=1
    from verify_separation import build_report
    separation=json.loads((F/'results/separation_verification.json').read_text())
    require('exact all-order separation replay',json.loads(json.dumps(build_report()))==separation)
    require('common order impossible',separation['filtration']['common_orders']==[])
    all_study={}
    for study,expected in [('clean',120),('sweep',30780)]:
        raw=pd.read_csv(R/f'bridgeless/{study}_samples.csv')
        all_study[study]=raw
        require(study+' rows',len(raw)==expected and not raw.duplicated(['n','seed','profile','method']).any())
        source=pd.read_csv(F/('results/extra_k/samples.csv' if study=='clean' else 'results/extra_k/rank_sweep/samples.csv'))
        if study=='sweep':source=source[source.n<=1000000]
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
    comparison=pd.read_csv(F/'results/extra_k/rank_sweep/n512_comparison_summary.csv')
    require('six-row n=512 comparison export',len(comparison)==6)
    bridge_summary=pd.read_csv(R/'bridgeless/sweep_summary.csv')
    for row in comparison.itertuples():
        q=bridge_summary[(bridge_summary.n==512)&(bridge_summary.profile==row.profile)&(bridge_summary.method==row.method)].iloc[0]
        require('n=512 comparison denominator',row.n==512 and row.samples==q.samples)
        for event in ('connection','bridgeless','plus3'):
            for field in ('rate','wilson_lo','wilson_hi'):
                require('n=512 comparison export values',abs(getattr(row,event+'_'+field+'_percent')-100*q[event+'_'+field])<1e-11)
    # Independent full Cartesian adjacency and low-link DFS on one clean
    # sample from every profile. The computation validates all 60 clean samples.
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
    clean_lookup={(int(row.seed),row.model,row.method,row.selector):
                  (float(row.e_PH_selected),float(row.e_PH_euclidean))
                  for row in reconstructed[reconstructed.delta==0].itertuples()}
    require('paired clean baseline scope',len(clean_lookup)==240)
    for key,g in paired.groupby(['model','delta','method','selector']):
        summaries=pd.read_csv(R/'euclidean/paired_summary.csv')
        q=summaries[(summaries.model==key[0])&(summaries.delta==key[1])&(summaries.method==key[2])&(summaries.selector==key[3])].iloc[0]
        differences=g.e_PH_selected-g.e_PH_euclidean
        require('paired baseline comparisons',q.samples==len(g) and abs(q.median_paired_difference-differences.median())<1e-12 and q.selected_better==(differences<0).sum() and q.euclidean_better==(differences>0).sum() and q.equal==(differences==0).sum())

        graph_change=np.array([float(row.e_PH_selected)-clean_lookup[(int(row.seed),row.model,row.method,row.selector)][0] for row in g.itertuples()])
        euclidean_change=np.array([float(row.e_PH_euclidean)-clean_lookup[(int(row.seed),row.model,row.method,row.selector)][1] for row in g.itertuples()])
        require('noise increment quantiles',all(abs(float(q[field+'_'+name+'_noise_increment'])-float(np.quantile(values,prob)))<1e-12
                for name,values in [('selected',graph_change),('euclidean',euclidean_change)]
                for field,prob in [('median',.5),('q1',.25),('q3',.75)]))
        changes=graph_change-euclidean_change
        require('paired noise increment comparisons',q.samples==20 and
                q.increment_selected_better==int((changes<0).sum()) and
                q.increment_euclidean_better==int((changes>0).sum()) and
                q.increment_equal==int((changes==0).sum()) and
                abs(q.median_paired_increment_difference-float(np.median(changes)))<1e-12)

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
    # Check the 36-row size-sweep summary independently against raw bridge rows.
    full=all_study['sweep']
    sweep=pd.read_csv(R/'size_sweep_summary.csv')
    require('size sweep coverage',len(sweep)==36 and set(sweep.n)=={128,512,4096,32768,262144,1000000})
    for row in sweep.itertuples():
        g=full[(full.n==row.n)&(full.profile==row.profile)&(full.method==row.method)]
        require('size sweep sample counts',row.samples==len(g))
        for key in ['k_conn','k_br','k_close']:
            require('size sweep event medians',getattr(row,'median_'+key)==float(np.median(g[key])))
        for event,column in [('connection','closed_at_connection'),('bridgeless','closed_at_bridgeless'),('plus3','closed_at_plus3')]:
            success=int(g[column].sum());lo,hi=wilson(success,len(g))
            require('size sweep closure rates and intervals',getattr(row,event+'_closed')==success and
                    abs(getattr(row,event+'_rate')-success/len(g))<1e-14 and
                    abs(getattr(row,event+'_wilson_lo')-lo)<1e-14 and abs(getattr(row,event+'_wilson_hi')-hi)<1e-14)
    # Its four primary sizes agree with the Table 3 counts and the Figure 4 export.
    table=pd.read_csv(F/'results/retention_table/retention_counts.csv')
    layers=pd.read_csv(F/'results/retention_extra_layers/retention.csv')
    require('size sweep primary rows',len(table)==24 and len(sweep[sweep.n<=32768])==24)
    for row in sweep[sweep.n<=32768].itertuples():
        q=table[(table.profile==row.profile)&(table.method==row.method)&(table.n==row.n)].iloc[0]
        require('size sweep equals Table 3 counts',q.samples==row.samples and q.closed_at_k_conn==row.connection_closed and
                q.closed_at_k_br==row.bridgeless_closed and q.closed_at_k_conn_plus3==row.plus3_closed)
        for j,event in [(0,'connection'),(3,'plus3')]:
            r=layers[(layers.profile==row.profile)&(layers.method==row.method)&(layers.n==row.n)&(layers.j==j)].iloc[0]
            require('size sweep equals extra-layer export',r.samples==row.samples and r.closed==getattr(row,event+'_closed'))
    cited={('uniform','OR_kNN'):1124,('localized_trough','mutual_kNN'):618}
    for (profile,method),closed in cited.items():
        r=sweep[(sweep.n==512)&(sweep.profile==profile)&(sweep.method==method)].iloc[0]
        require('size sweep cited bridge-free counts',r.bridgeless_closed==closed and r.samples==1500)
    # Recompute the closure-error medians from the per-sample clean records.
    clean=F/'results/clean_closure'
    cells=pd.read_csv(clean/'cells.csv');knn=pd.read_csv(clean/'knn_event_checks.csv')
    exact_rows=pd.concat([cells[(cells.method=='R2_global')&(cells.budget_label=='completion')],knn[knn.landmark=='own_coverage']],ignore_index=True)
    require('exact closure error inputs',len(exact_rows)==180)
    exact=pd.read_csv(F/'results/exact_closure_errors.csv');grid=pd.read_csv(clean/'first_closed_errors.csv')
    selected=pd.read_csv(clean/'closure_selection.csv')
    require('closure error rows',len(exact)==9 and len(grid)==9)
    for e in exact.itertuples():
        g=exact_rows[(exact_rows.profile==e.profile)&(exact_rows.method==e.method)]
        require('exact closure error medians',len(g)==20 and all(abs(getattr(e,k)-float(np.median(g[k])))<1e-15 for k in ['e_inf','e_2','p1']))
    for q in grid.itertuples():
        g=selected[(selected.profile==q.profile)&(selected.method==q.method)]
        require('first closed grid error medians',len(g)==q.samples==20 and
                all(abs(getattr(q,'median_'+k)-float(np.median(g['first_covered_'+s])))<1e-15 for k,s in [('e_inf','e_inf'),('e_2','e_2'),('e_PH','p1')]))
    calibration=json.loads((R/'small_study_calibration.json').read_text())['methods']
    audit=json.loads((F/'results/number_audit.json').read_text())['paired_seed_calibration']
    for a in audit:
        c=calibration[a['method']]
        require('paired calibration',c['original_successes']==a['observed'] and
                abs(c['mean_successes']-a['expected'])<1e-12 and abs(c['upper_tail']-a['probability_at_least_observed'])<1e-12)
    literature=json.loads((R/'literature_checks.json').read_text())
    require('connectivity claim scope',not literature['connectivity_stopping']['stronger_practice_claim_supported'])
    require('spectral parameter record',literature['spectral_graph']['parameters']==dict(k=[15,100],diffusion_time=[8,64]))
    return dict(status='passed',checks_passed=count,bridge_method_records=30900,
                distinct_euclidean_diagrams=220,grid_oracle_samples=720,
                scope='Exact orientation certificate; saved bridge events and independent adjacency/DFS replays; all Euclidean input hashes, bars and bottleneck scores; grid oracles, paired comparisons and quantiles.')


if __name__=='__main__':
    F=Path(__file__).resolve().parents[1]
    report=verify_comparisons(F)
    (F/'results/comparisons/verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
