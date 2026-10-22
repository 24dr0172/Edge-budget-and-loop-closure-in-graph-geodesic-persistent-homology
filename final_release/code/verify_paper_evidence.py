"""Saved-data checks for the claims of the paper, used by every verifier mode."""
from pathlib import Path
import json
import csv
from fractions import Fraction
from rank_table import verify_table
from verify_small_examples import build_report


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
    additional={}
    if include_additional:
        from verify_reported_claims import verify_reported_claims
        from verify_extra_layers_and_peanut import verify_extra_layers_and_peanut
        from build_practical_choices import verify_practical_choices
        from verify_peanut_neck import verify_peanut_neck
        additional=dict(reported_claims=verify_reported_claims(root,require),
                        extra_layers_and_peanut=verify_extra_layers_and_peanut(root,require),
                        practical_choices=verify_practical_choices(root,require),
                        peanut_neck=verify_peanut_neck(root,require))
    return dict(status='passed',table=table,full_diagram_crosscheck=clean,
                practical_claims=practical,noncircular_locality_check=noncircular,
                retention_cost_witness=retention_cost_witness,radius_integral_checks=12,
                extra_layer_rows=144,additional=additional,
                examples='Live-threshold separation, natural ellipse, hexagon, shortcut and ties; Appendix C.',
                scope='Focused-paper evidence and shared certificate exports; additional cohorts only in complete mode.')
