"""Check the fixed-k, bridge-free, heuristic-model and diagram-cost numbers of the paper from saved data."""
from pathlib import Path
import csv
import json
import math


def build_report(root=None, require=None):
    root=Path(root) if root is not None else Path(__file__).resolve().parents[1]
    checks=[]
    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)
        if require is not None:
            require(name, condition)
    def rows(name):
        with (root/'results'/name).open() as stream:
            return list(csv.DictReader(stream))
    raw=rows('extra_k/rank_sweep/samples.csv')
    fixed=[]
    for saved in rows('extra_k/rank_sweep/fixed_k_summary.csv'):
        n,k=int(saved['n']),int(saved['k'])
        group=[r for r in raw if r['profile']==saved['profile'] and
               r['method']==saved['method'] and int(r['n'])==n]
        count=sum(int(r['k_close'])<=k for r in group)
        check('practical_fixed_'+saved['method']+f'_{n}_{k}',
              len(group)==len({r['seed'] for r in group})==int(saved['samples']) and
              all(float(r['largest_gap'])<1/3 for r in group) and
              count==int(saved['closed']) and
              abs(100*count/len(group)-float(saved['closure_percent']))<1e-11)
        fixed.append(dict(profile=saved['profile'],method=saved['method'],n=n,k=k,
                          samples=len(group),closed=count,closure_percent=100*count/len(group)))
    check('practical_fixed_headlines',[(r['n'],r['closed'],r['samples']) for r in fixed
          if r['method']=='mutual_kNN' and r['k']==10 and r['n'] in (512,4096)]==
          [(512,593,1500),(4096,1,1500)])

    bridge=rows('comparisons/bridgeless/sweep_samples.csv')
    bridge_summary=rows('comparisons/bridgeless/sweep_summary.csv')
    references={(r['seed'],r['profile'],r['method']):r for r in raw if int(r['n'])==512}
    n512_rows=[]
    for summary in bridge_summary:
        if int(summary['n'])!=512:
            continue
        profile,method=summary['profile'],summary['method']
        group=[r for r in bridge if int(r['n'])==512 and r['profile']==profile and r['method']==method]
        count=sum(int(r['k_br'])==int(r['k_close']) for r in group)
        check('practical_bridge_'+profile+'_'+method,
              len(group)==len({r['seed'] for r in group})==int(summary['samples'])==1500 and
              count==int(summary['bridgeless_closed']) and
              abs(count/1500-float(summary['bridgeless_rate']))<1e-14 and
              all(int(r['k_conn'])<=int(r['k_br'])<=int(r['k_close']) and
                  (r['closed_at_bridgeless']=='True')==(int(r['k_br'])==int(r['k_close'])) and
                  all(r[key]==references[r['seed'],profile,method][key] for key in ('k_conn','k_close'))
                  for r in group))
        n512_rows.append(dict(profile=profile,method=method,n=512,samples=1500,
                             closed=count,closure_percent=100*count/1500))
    check('practical_bridge_headlines',
          {(r['profile'],r['method']):r['closed'] for r in n512_rows}[
              'uniform','OR_kNN']==1124 and
          {(r['profile'],r['method']):r['closed'] for r in n512_rows}[
              'localized_trough','mutual_kNN']==618)

    from rank_sweep import tails, iid_extra_cdf
    distributions=tails()
    predictions=[]
    for row in rows('extra_k/rank_sweep/iid_predictions.csv'):
        if int(row['n']) not in (128,512,4096,32768):
            continue
        value=iid_extra_cdf(int(row['n']),distributions[row['method']],0)
        check('practical_independent_model_'+row['method']+'_'+row['n'],
              abs(value-float(row['cdf_0']))<1e-14)
        predictions.append(dict(n=int(row['n']),method=row['method'],probability=value))
    check('practical_model_scope',len(predictions)==8)

    clouds=rows('clean_closure/clouds.csv')
    cloud_keys={(r['seed'],r['profile']) for r in clouds}
    check('practical_error_sample_scope',len(clouds)==len(cloud_keys)==60 and
          all(int(r['n'])==512 for r in clouds))
    baseline=[float(r['oracle_longest_lifetime'])/2 for r in clouds]
    own={(r['seed'],r['profile'],r['method'],r['budget_label'],r['M']):r
         for r in rows('clean_closure/own_closure_cells.csv')}
    grid=[]
    for row in rows('clean_closure/cells.csv'):
        key=tuple(row[k] for k in ('seed','profile','method','budget_label','M'))
        if row['budget_label'].startswith('a') and row['connected']=='True' and own[key]['covered']=='True':
            check('practical_postclosure_budget_'+':'.join(key),
                  int(row['M'])>=int(own[key]['M_cov']) and float(row['p1'])>=0 and
                  (row['seed'],row['profile']) in cloud_keys)
            grid.append(row)
    maximum=max(float(r['p1']) for r in grid)
    audit=json.loads((root/'results/number_audit.json').read_text())
    check('practical_saved_error_maximum',abs(maximum-audit['maximum_postclosure_grid_e_PH'])<1e-14 and
          f'{maximum:.3f}'=='0.044')
    check('practical_empty_cost_range',f'{min(baseline):.3f}'=='0.244' and
          f'{max(baseline):.3f}'=='0.324')
    # Stored endpoints are float32; use the release's endpoint tolerance.
    prefixes=rows('useful_diagram/saved_prefixes.csv')
    for row in prefixes:
        g,b,diameter=float(row['g']),float(row['true_death']),float(row['diameter'])
        chi=lambda x: math.sin(math.pi*x)/math.pi
        ceiling=max(g-chi(g),b-chi(b))/diameter
        check('practical_ceiling_'+':'.join(row[k] for k in ('seed','profile','method','M')),
              int(row['bars'])==1 and float(row['e_PH'])<=ceiling+1e-6)
    limit=2/3-math.sqrt(3)/math.pi
    ratio=(1/3)/limit
    check('practical_limit_values',f'{limit:.3f}'=='0.115' and f'{ratio:.1f}'=='2.9')
    return dict(status='passed',checks=len(checks),
                fixed_k=fixed,bridge_free_n512=n512_rows,independent_key_model=predictions,
                diagram_error=dict(samples=60,n=512,stored_postclosure_grid_cells=len(grid),
                                   maximum=maximum,empty_cost_minimum=min(baseline),
                                   empty_cost_maximum=max(baseline),
                                   saved_prefix_ceiling_checks=len(prefixes)),
                dense_circle_limit=dict(ceiling=limit,empty_cost=1/3,cost_ratio=ratio),
                scope='Saved-data checks of fixed-k, bridge-free, heuristic-model and cost claims. '
                      'The model assumes independent keys; the cost ratio is a dense-sample limit. '
                      'No measurements or persistence diagrams are generated.')


def verify_practical(root, require):
    report=build_report(root,require)
    saved=Path(root)/'results/practical_claims_verification.json'
    require('check_practical_claims_report_replay',saved.is_file() and json.loads(saved.read_text())==report)
    return report


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    report=build_report(root)
    (root/'results/practical_claims_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f"Practical claims passed: {report['checks']} saved-data checks.")
