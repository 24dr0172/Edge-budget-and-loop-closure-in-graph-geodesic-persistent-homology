"""Recheck saved error comparisons at two tolerances, from frozen diagrams only."""
from pathlib import Path
import csv,json,math
ROOT=Path(__file__).resolve().parents[1]
def verify(root=ROOT):
    root=Path(root)
    read=lambda f:list(csv.DictReader((root/'results'/f).open()))
    peanuts=read('peanut_neck/trials.csv')
    wide=[r for r in peanuts if float(r['alpha'])>=.6]
    assert len(wide)==180
    residual=max(abs(float(r['euclidean_error'])-float(r['P_star'])) for r in wide)
    assert residual==0
    for seed in {r['seed'] for r in wide}:
        assert len({float(r['P_star']) for r in wide if r['seed']==seed})==1
    rows=read('useful_diagram/saved_prefixes.csv')
    differences=[]; ceiling_slack=[]; series={}
    for r in rows:
        assert int(r['bars'])==1
        b,d,tb,td,diam=[float(r[k]) for k in ['birth','death','true_birth','true_death','diameter']]
        exact=min(max(abs(b-tb),abs(d-td)),max((d-b)/2,(td-tb)/2))/diam
        differences.append(abs(exact-float(r['e_PH'])))
        chi=lambda x:math.sin(math.pi*x)/math.pi
        ceiling=max(float(r['g'])-chi(float(r['g'])),td-chi(td))/diam
        ceiling_slack.append(float(r['e_PH'])-ceiling)
        series.setdefault(tuple(r[k] for k in ['seed','profile','method']),[]).append(r)
    steps=[]
    for group in series.values():
        q=sorted(group,key=lambda r:int(r['M']))
        steps.extend(float(b['e_PH'])-float(a['e_PH']) for a,b in zip(q,q[1:]))
    outcomes={str(t):dict(matching_failures=sum(v>t for v in differences),
        ceiling_failures=sum(v>t for v in ceiling_slack),
        monotonicity_failures=sum(v < -t for v in steps)) for t in [5e-8,1e-6]}
    assert all(not any(v.values()) for v in outcomes.values())
    return dict(status='passed',peanut_rows_checked=len(wide),euclidean_minus_empty_cost_max=residual,
        clean_prefixes_checked=len(rows),clean_series=len(series),tolerance_comparisons=outcomes,
        maximum_matching_residual=max(differences),minimum_monotonic_step=min(steps),
        scope='Frozen diagrams and errors only; no regeneration; the generation tolerances in code/config.json apply.')
if __name__=='__main__':
    result=verify();(ROOT/'results/tolerance_recheck.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
