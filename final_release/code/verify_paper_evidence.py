"""Saved-data checks for the claims of the paper, used by every verifier mode."""
from pathlib import Path
import json
import csv
from fractions import Fraction
from rank_table import verify_table
from verify_small_examples import build_report


def verify_peanut_bound(root, require):
    """Check both fixed-graph bounds against independently generated metric ratios and winding thresholds."""
    import numpy as np
    root=Path(root); folder=root/'results/peanut_neck'
    raw=list(csv.DictReader((folder/'trials.csv').open()))
    report=json.loads((folder/'closure_error_bound.json').read_text())
    records={(v['alpha'],v['seed'],v['method']):v for v in report['records']}
    winding={(float(r['alpha']),int(r['seed']),r['method']):r
             for r in csv.DictReader((root/'results/winding_check/records.csv').open())}
    require('paper_peanut_bound_inventory',len(records)==len(raw)==300 and
            report['selected_edge_ratio_checks']==300 and report['status']=='passed')
    for row in raw:
        a,seed,method=float(row['alpha']),int(row['seed']),row['method']
        v=records[a,seed,method]
        with np.load(folder/'checkpoints'/f'a{a:g}_seed{seed}.npz') as saved:
            s=saved['s']; gap=np.diff(np.r_[s,s[0]+1]).max()
            separation=np.abs(s[:,None]-s[None,:])
            diameter=float(np.minimum(separation,1-separation).max())
            bstar=float(saved['reference_bars'][0,1])
            lam=float(row['closure_min_metric_ratio'])
            bound=max(gap-float(row['B']),(1-lam)*bstar)/diameter
            eligible=row['closure_supported']=='True' and float(row['B'])<lam*bstar
            require(f'paper_peanut_bound_{a}_{seed}_{method}',
                abs(v['lambda_G']-lam)<5e-12 and abs(v['U_G']-bound)<1e-11 and
                v['bound_applicable']==eligible and
                v['strict_bound_improvement']==(eligible and bound<float(row['euclidean_error'])))
            w=winding[a,seed,method]
            W_G=float(w['W_G']); support=row['closure_supported']=='True'
            wbound=min(max(gap-float(row['B']),bstar-max(float(row['B']),W_G)),(bstar-float(row['B']))/2)/diameter
            require(f'paper_peanut_winding_bound_{a}_{seed}_{method}',
                abs(v['W_G']-W_G)<1e-15 and abs(v['winding_bound']-wbound)<1e-11 and
                v['winding_survival']==(W_G>float(row['B'])) and v['winding_bound_applicable']==support and
                v['winding_bound_improvement']==(support and wbound<float(row['euclidean_error'])) and
                (not support or float(row['closure_error'])<=wbound+5e-8/diameter) and
                (not eligible or W_G>=lam*bstar-5e-8))
            if eligible:
                bars=saved[f'{method}_closure_bars']
                require(f'paper_peanut_death_bound_{a}_{seed}_{method}',len(bars)==1 and
                    lam*bstar-5e-8<=bars[0,1]<=bstar+5e-8 and
                    float(row['closure_error'])<=bound+5e-8/diameter)
    groups=list(csv.DictReader((folder/'closure_error_bound.csv').open()))
    require('paper_peanut_summary_rows',len(groups)==len(report['groups'])==15)
    for q in report['groups']:
        selected=[r for r in raw if float(r['alpha'])==q['alpha'] and r['method']==q['method']]
        paired=[v for v in report['records'] if v['alpha']==q['alpha'] and v['method']==q['method']]
        export=next(r for r in groups if float(r['alpha'])==q['alpha'] and r['method']==q['method'])
        require(f"paper_peanut_summary_{q['alpha']}_{q['method']}",
            q['samples']==len(selected)==20 and
            q['supported']==sum(r['closure_supported']=='True' for r in selected) and
            q['closure_single']==sum(int(r['closure_bars'])==1 for r in selected) and
            q['strict_bound_improvement']==sum(v['strict_bound_improvement'] for v in paired) and
            q['edgewise_survival']==sum(v['survival'] for v in paired) and
            q['winding_survival']==sum(v['winding_survival'] for v in paired) and
            q['winding_bound_improvement']==sum(v['winding_bound_improvement'] for v in paired) and
            all(abs(q[field+'_'+end]-{'min':min,'max':max}[end](float(r[field]) for r in selected))<1e-14
                for field in ('euclidean_error','closure_error') for end in ('min','max')) and
            all(str(value) == export[key] or (value is None and export[key]=='') for key,value in q.items()))
    return dict(status='passed',closure_graphs=300,method_by_width_rows=15,
                source='Frozen diagrams and all-pairs metric ratios; new builder checks selected edges.')


def verify_clean_crosscheck(root, require):
    """Check the 60 saved full-diagram samples; generate no diagrams."""
    folder=Path(root)/'results/clean_closure'
    def read(name):
        with (folder/name).open() as stream:
            return list(csv.DictReader(stream))
    clouds=read('clouds.csv')
    keys={(r['seed'],r['profile']) for r in clouds}
    require('paper_full_diagram_sample_scope',len(clouds)==len(keys)==60 and
            all(int(r['n'])==512 and float(r['g_max'])<1/3 for r in clouds))
    cells=read('cells.csv')
    knn=read('knn_event_checks.csv')
    events=read('event_table.csv')
    selected=[]
    for method in ('R2_global','OR_kNN','mutual_kNN'):
        for stage,knn_stage in [('first_connected','first_connected'),
                                ('pre_completion','pre_coverage_last'),('completion','own_coverage')]:
            source=cells if method=='R2_global' else knn
            column='budget_label' if method=='R2_global' else 'landmark'
            value=stage if method=='R2_global' else knn_stage
            rows=[r for r in source if r['method']==method and r[column]==value]
            require(f'paper_full_diagram_inventory_{method}_{stage}',
                    len(rows)==60 and {(r['seed'],r['profile']) for r in rows}==keys)
            for row in rows:
                require(f"paper_full_diagram_{method}_{stage}_{row['seed']}_{row['profile']}",
                        row['connected']=='True' and int(row['n_components'])==1 and
                        Fraction(row['est_h1_bars']).denominator==1 and
                        (Fraction(row['est_h1_bars'])>0)==(stage=='completion'))
            aggregate=[r for r in events if r['method']==method and r['event']==stage]
            require(f'paper_full_diagram_aggregate_{method}_{stage}',len(aggregate)==1 and
                    int(aggregate[0]['attempted'])==int(aggregate[0]['connected'])==60 and
                    int(aggregate[0]['empty'])==(0 if stage=='completion' else 60) and
                    int(aggregate[0]['nonempty'])==(60 if stage=='completion' else 0))
            selected.extend(rows)
    return dict(status='passed',samples=60,n=512,methods=3,event_records=len(selected),
                expected_diagrams=['empty','empty','nonempty'],
                events=['first connected exact prefix','immediately before closure','closure'],
                source='Saved full-diagram event records; no n=512 persistence recomputation.')


def verify_paper_evidence(root, require, include_additional=False):
    from verify_ellipse_window import verify as verify_window
    window=verify_window(root)
    saved_window=json.loads((Path(root)/'results/certified_ellipse_window/certificate.json').read_text())
    require('certified_ellipse_window',window==saved_window and window['status']=='passed')
    from tolerance_recheck import verify as tolerance_recheck
    require('tolerance_recheck',tolerance_recheck(root)['status']=='passed')
    table = verify_table(root, require)
    path = Path(root)/'results/small_examples_verification.json'
    require('paper_examples_report_present', path.is_file())
    saved = json.loads(path.read_text())
    current = json.loads(json.dumps(build_report()))
    require('paper_examples_replay', saved == current)
    require('paper_examples_status', saved['status'] == 'passed')
    require('paper_live_thresholds', [q['threshold_float'] for q in saved['separation']['proof_thresholds']] == [1.5, 1.9]
            and not saved['separation']['common_orders'])
    natural=Path(root)/'results/natural_instance_verification.json'
    require('paper_natural_ellipse_report_present',natural.is_file())
    require('paper_natural_ellipse_report_replay',json.loads(natural.read_text())==saved['natural_ellipse'])
    require('paper_natural_ellipse_live_noncyclic',saved['natural_ellipse']['status']=='passed' and
            saved['natural_ellipse']['proof_thresholds'][1]['admissible_orders']==[])
    from check_practical_claims import verify_practical
    practical=verify_practical(root,require)
    clean=verify_clean_crosscheck(root,require)
    from verify_noncircular_locality import verify_noncircular
    noncircular=verify_noncircular(root,require)
    from verify_retention_cost_witness import verify_retention_cost_witness
    retention_cost_witness=verify_retention_cost_witness(root,require)
    import math
    from scipy.integrate import quad
    from scipy.special import gammaln
    for n in (4,8,23,128):
        for alpha in (1.25,1.5,2.):
            integral=n*alpha*quad(lambda u:u**(alpha-1)*(1-u)**(n-1),0,1,epsabs=1e-13)[0]
            gamma=math.exp(gammaln(alpha+1)+gammaln(n+1)-gammaln(n+alpha))
            require(f'radius_beta_integral_{n}_{alpha}',abs(integral-gamma)<2e-11)
    from retention_extra_layers import wilson
    records=list(csv.DictReader((Path(root)/'results/extra_k/rank_sweep/samples.csv').open()))
    export=list(csv.DictReader((Path(root)/'results/retention_extra_layers/retention.csv').open()))
    require('focused_extra_layers_inventory',len(export)==144)
    for row in export:
        n,j=int(row['n']),int(row['j'])
        group=[v for v in records if int(v['n'])==n and v['profile']==row['profile'] and v['method']==row['method']]
        closed=sum(int(v['k_close'])<=int(v['k_conn'])+j for v in group)
        lo,hi=wilson(closed,len(group))
        require('focused_extra_layers_'+str((n,row['profile'],row['method'],j)),
                int(row['closed'])==closed and int(row['samples'])==len(group) and
                abs(float(row['rate'])-closed/len(group))<1e-14 and
                abs(float(row['wilson_lo'])-lo)<1e-14 and abs(float(row['wilson_hi'])-hi)<1e-14)
    from verify_peanut_neck import verify_peanut_neck
    peanut=verify_peanut_neck(root,require,replay=include_additional)
    peanut_bound=verify_peanut_bound(root,require)
    from verify_first_nonempty import verify as verify_stopping
    stopping=verify_stopping(root,require,replay=include_additional)
    additional={}
    if include_additional:
        from verify_reported_claims import verify_reported_claims
        from verify_extra_layers_and_peanut import verify_extra_layers_and_peanut
        from build_practical_choices import verify_practical_choices
        additional=dict(reported_claims=verify_reported_claims(root,require),
                        extra_layers_and_peanut=verify_extra_layers_and_peanut(root,require),
                        practical_choices=verify_practical_choices(root,require))
    return dict(status='passed',table=table,full_diagram_crosscheck=clean,
                practical_claims=practical,noncircular_locality_check=noncircular,
                retention_cost_witness=retention_cost_witness,radius_integral_checks=12,
                extra_layer_rows=144,peanut_neck=peanut,peanut_bound=peanut_bound,first_nonempty=stopping,additional=additional,
                examples='Live-threshold separation, natural ellipse, hexagon, shortcut and ties; finite-example appendices.',
                scope='Focused-paper evidence including the saved varying-waist study and edgewise bounds; other additional cohorts only in complete mode.')
