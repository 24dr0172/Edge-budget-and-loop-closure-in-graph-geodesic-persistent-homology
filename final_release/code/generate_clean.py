"""Generate all clean-circle experiments from the stated seeds and formulas."""
from pathlib import Path
import argparse,concurrent.futures,json,time,math
import numpy as np
import pandas as pd
import graph_methods as m
from experiment_core import ROOT,RESULTS,CONFIG,table,dump,grid_values,all_orders,evaluate,rank_containment,interval_row,transitions
def clean_job(task):
    seed,profile=task;start=time.monotonic();folder=RESULTS/'clean_closure/checkpoints';folder.mkdir(parents=True,exist_ok=True)
    record=folder/f'{profile}_{seed}.json'
    if record.exists():return seed,profile,'resumed'
    s,points=m.sample_cloud(seed,profile,n=CONFIG['n']);n=len(s)
    true=m.intrinsic(s);reference=m.diagram(true);diameter=float(true.max())
    assert reference.shape==(1,2)
    chord,iu,ju,length,orders,keys=all_orders(points,seed)
    conn={name:m.first_connected(order,iu,ju,n) for name,order in orders.items()}
    close={name:m.closure(order,iu,ju,n) for name,order in orders.items()}
    consecutive=(ju-iu==1)|((iu==0)&(ju==n-1))
    gaps=np.diff(np.r_[s,s[0]+1]);g=float(gaps.max());second=float(np.sort(gaps)[-2])
    life=float(np.max(reference[:,1]-reference[:,0])/diameter)
    cloud=dict(seed=seed,profile=profile,noise='clean',n=n,g_max=g,g_second=second,
               oracle_diameter=diameter,oracle_longest_lifetime=life,R2_M_cover=close['R2_global'],R2_M_pre=close['R2_global']-1)
    cells=[];events=[];ranks=[];certs=[];selection=[];monotone=[];prefixes=[];steps=[];control={};cache={};arrays=dict(latent=s,observed=points,reference=reference)
    order_arrays={'latent':s,'observed':points,'reference':reference}
    def calc(method,count):
        key=(method,int(count))
        if key not in cache:cache[key]=evaluate(chord,iu,ju,length,orders[method],int(count),true,reference,conn[method])
        return cache[key]
    base_budgets={}
    for method in m.METHODS:
        order=orders[method];cloud[method+'_M_conn']=conn[method]
        budgets=dict(grid_values());budgets['first_connected']=conn[method]
        if method=='R2_global':
            budgets.update(pre_completion=close[method]-1,completion=close[method])
        else:
            k=int(keys[method][order[conn[method]-1]]);mk=int(np.sum(keys[method]<=k))
            cloud[method+'_k_conn']=k;budgets['smallest_connected_k']=mk
            topkeys=np.sort(keys[method][consecutive]);positions=np.zeros((n,n),np.int64)
            positions[iu[order],ju[order]]=np.arange(1,len(order)+1);positions+=positions.T
            cp=np.sort(positions[np.arange(n),np.roll(np.arange(n),-1)])
            bad,checked,anti,nearest=rank_containment(positions,s)
            ranks.append(dict(seed=seed,profile=profile,method=method,n=n,g_max=g,g_second=second,
                oracle_longest_lifetime=life,M_conn=conn[method],M_second_to_last_consecutive=int(cp[-2]),
                M_own_cover=close[method],k_conn=k,k_second_largest_consecutive_key=int(topkeys[-2]),
                k_largest_consecutive_key=int(topkeys[-1]),M_of_k_conn=mk,integer_k_covers=k>=topkeys[-1],
                closure_violations=bad,closure_steps_checked=checked,antipodal_pairs=anti,min_abs_arc_minus_half=nearest,
                n_consecutive_keys_tied_at_top=int(np.sum(topkeys==topkeys[-1]))))
            for label,count in [('first_connected',conn[method]),('pre_coverage_last',close[method]-1),
                                ('smallest_connected_k',mk),('own_coverage',close[method])]:
                score,bars=calc(method,count)
                events.append(dict(seed=seed,profile=profile,method=method,landmark=label,M=count,**score))
        base_budgets[method]=set(budgets.values())
        gridrows=[]
        for label,count in budgets.items():
            score,bars=calc(method,count)
            row=dict(seed=seed,profile=profile,noise='clean',method=method,budget_label=label,M=count,a_M=count/(n*math.log(n)),**score)
            cells.append(row)
            if label.startswith('a'):gridrows.append(row)
        cover,holds,examined,witness=m.certificate(s,chord,order,iu,ju)
        certs.append(dict(seed=seed,profile=profile,method=method,M_cov=cover,holds=holds,nonconsecutive_edges_checked=examined,**witness))
        covered=[x for x in gridrows if x['M']>=close[method]]
        nonempty=[x for x in gridrows if x['connected'] and x['est_h1_bars']>0]
        first=min(covered,key=lambda x:x['M']);firstnon=min(nonempty,key=lambda x:x['M'])
        selection.append(dict(seed=seed,profile=profile,method=method,certificate_holds=holds,M_cov=close[method],
            first_covered_M=first['M'],first_nonempty_M=firstnon['M'],nonempty_matches_covered=first['M']==firstnon['M'],
            first_covered_e_inf=first['e_inf'],first_covered_e_2=first['e_2'],first_covered_p1=first['p1']))
        mono=dict(seed=seed,profile=profile,method=method,M_cov=close[method],covered_grid_budgets=len(covered))
        covered=sorted(covered,key=lambda x:x['M'])
        for field in ['e_inf','e_2','p1']:
            differences=np.diff([x[field] for x in covered])
            mono[field+'_decreases']=int(np.sum(differences < -5e-8))
            mono[field+'_minimum_difference']=float(differences.min()) if len(differences) else None
        monotone.append(mono)
        series=[]
        for count in sorted({c for c in base_budgets[method] if c>=close[method]}|{close[method],len(order)}):
            score,bars=calc(method,count);arrays[f'{method}_M{count}']=bars
            row=interval_row(s,reference,diameter,method,count,close[method],bars)
            row.update(seed=seed,profile=profile,is_saved_budget=count in base_budgets[method])
            assert len(bars)==1
            series.append(row);prefixes.append(row)
        for row in transitions(series):row.update(seed=seed,profile=profile);steps.append(row)
    for method in ['R1_local','R3_permuted']:
        cover,holds,examined,witness=m.certificate(s,chord,orders[method],iu,ju)
        score,bars=calc(method,conn[method]);order_arrays[method+'_first_connected']=bars
        if method=='R1_local':
            certs.append(dict(seed=seed,profile=profile,method=method,M_cov=cover,holds=holds,nonconsecutive_edges_checked=examined,**witness))
        else:
            first_closed=next((count for _,count in grid_values() if count>=cover),None)
            first_nonempty=None
            for _,count in grid_values():
                if count<conn[method]:continue
                q,bar=calc(method,count);order_arrays[f'{method}_M{count}']=bar
                if len(bar):first_nonempty=count;break
            control=dict(seed=seed,profile=profile,method=method,short_hops_first_passes=holds,M_close=cover,
                M_first_connected=conn[method],first_connected_H1_bars=len(bars),first_closed_grid_M=first_closed,
                first_nonempty_grid_M=first_nonempty,first_nonempty_matches_first_closed_grid=first_nonempty==first_closed)
    output=RESULTS/'useful_diagram/checkpoints';output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/f'{profile}_{seed}.npz',**arrays)
    output=RESULTS/'order_checks/checkpoints';output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/f'{profile}_{seed}.npz',**order_arrays)
    dump(record,dict(seed=seed,profile=profile,cloud=cloud,cells=cells,events=events,ranks=ranks,certs=certs,
        selection=selection,monotone=monotone,prefixes=prefixes,steps=steps,control=control,seconds=time.monotonic()-start))
    return seed,profile,round(time.monotonic()-start,2)
def aggregate():
    data=[json.loads(p.read_text()) for p in sorted((RESULTS/'clean_closure/checkpoints').glob('*.json'))]
    fields={'clouds.csv':'cloud','cells.csv':'cells','knn_event_checks.csv':'events','rank_landmarks.csv':'ranks',
            'closure_selection.csv':'selection','covered_monotonicity.csv':'monotone'}
    for name,key in fields.items():table('clean_closure/'+name,[r for d in data for r in (d[key] if isinstance(d[key],list) else [d[key]])])
    cert=table('order_checks/certificates.csv',[r for d in data for r in d['certs']])
    table('order_checks/shuffled_control.csv',[d['control'] for d in data])
    failure=cert[~cert.holds]
    dump(RESULTS/'order_checks/failure_witness.json',dict(status='computed',failures=failure.where(pd.notna(failure),None).to_dict('records')))
    prefixes=table('useful_diagram/saved_prefixes.csv',[r for d in data for r in d['prefixes']])
    trans=table('useful_diagram/saved_transitions.csv',[r for d in data for r in d['steps']])
    rows=pd.read_csv(RESULTS/'clean_closure/cells.csv')
    selection=pd.read_csv(RESULTS/'clean_closure/closure_selection.csv')
    merged=rows.merge(selection[['seed','profile','method','M_cov']],on=['seed','profile','method'],validate='many_to_one')
    merged['covered']=merged.M>=merged.M_cov
    table('clean_closure/own_closure_cells.csv',merged)
    latents={}
    for d in data:
        with np.load(RESULTS/f"useful_diagram/checkpoints/{d['profile']}_{d['seed']}.npz") as z:
            latents[f"{d['seed']}_{d['profile']}"]=z['latent']
    np.savez_compressed(RESULTS/'clean_closure/latent_coordinates.npz',**latents)
    dump(RESULTS/'useful_diagram/execution_summary.json',dict(status='computed_from_seeds',samples=len(data),
         prefixes=len(prefixes),transitions=len(trans),single_bar_failures=int((prefixes.bars!=1).sum()),
         maximum_birth_deviation=float(prefixes.birth_deviation.max()),raw_error_drops=int(trans.raw_drop.sum()),
         capped_error_drops=int(trans.clipped_drop.sum()),largest_full_diagram_n=CONFIG['n']))
    print('Clean aggregates:',len(data),'samples,',len(rows),'cells,',len(prefixes),'closed diagrams.',flush=True)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=6);ap.add_argument('--limit',type=int);ap.add_argument('--aggregate-only',action='store_true');args=ap.parse_args()
    if not args.aggregate_only:
        tasks=[(s,p) for s in CONFIG['clean']['seeds'] for p in CONFIG['clean']['profiles']]
        if args.limit:tasks=tasks[:args.limit]
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
            for i,result in enumerate(ex.map(clean_job,tasks),1):print('Clean',i,'/',len(tasks),result,flush=True)
    aggregate()
if __name__=='__main__':main()
