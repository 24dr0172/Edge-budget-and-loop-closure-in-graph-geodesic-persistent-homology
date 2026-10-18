"""Generate the ellipse coordinates and exact connectivity/closure checks from seeds."""
import argparse,concurrent.futures,json,math,time
import numpy as np
import pandas as pd
import paper1_methods as m
from experiment_core import ROOT,RESULTS,CONFIG,dump,table,ellipse_points,evaluate
def coordinates():
    arrays={}
    for seed in CONFIG['ellipse']['seeds']:
        for profile in CONFIG['ellipse']['profiles']:
            s,_=m.sample_cloud(seed,profile,n=CONFIG['n'])
            arrays[f'{seed}_{profile}_latent']=s
            arrays[f'{seed}_{profile}_xy']=ellipse_points(s)
    folder=RESULTS/'ellipse_checks';folder.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(folder/'coordinates.npz',**arrays)
    raw=12*np.sqrt(15)/25;perimeter_bound=np.pi*np.sqrt(10)
    dump(folder/'analytic_constants.json',dict(chord_threshold=.18,raw_min_normal_chord=raw,
        perimeter_upper=perimeter_bound,normal_chord_lower_after_normalization=raw/perimeter_bound,
        global_chord_arc_lower_factor=1/np.pi,half_perimeter_chord_lower_after_normalization=2/perimeter_bound,
        ellipse_axis_ratio='2:1',normalization='Unit perimeter; analytic perimeter upper bound pi*sqrt(10).'))
    print('Generated all ellipse coordinates from seeds and a 400001-point arclength grid.',flush=True)
def job(task):
    seed,profile=task;start=time.monotonic();folder=RESULTS/'ellipse_checks/checkpoints';folder.mkdir(parents=True,exist_ok=True)
    target=folder/f'{profile}_{seed}.json'
    if target.exists():return seed,profile,'resumed'
    s,_=m.sample_cloud(seed,profile,n=CONFIG['n']);points=ellipse_points(s);n=len(s)
    true=m.intrinsic(s);reference=m.diagram(true)
    chord,iu,ju,length,orders,keys=m.build_orders(points)
    gaps=np.diff(np.r_[s,s[0]+1]);g=float(gaps.max());B=float(chord[np.arange(n),np.roll(np.arange(n),-1)].max())
    mask=chord.copy();np.fill_diagonal(mask,np.inf);ordered=np.sort(mask,axis=1)
    consecutive=(ju-iu==1)|((iu==0)&(ju==n-1))
    rows=[];arrays={};cache={}
    for method in m.METHODS:
        order=orders[method];conn=m.first_connected(order,iu,ju,n);close=m.closure(order,iu,ju,n)
        cover,holds,checked,witness=m.certificate(s,chord,order,iu,ju)
        positions=np.zeros((n,n),np.int64);positions[iu[order],ju[order]]=np.arange(1,len(order)+1);positions+=positions.T
        cp=np.sort(positions[np.arange(n),np.roll(np.arange(n),-1)])
        bound=B if method=='R2_global' else float(ordered[:,int(keys[method][consecutive].max())-1].max())
        spans=True
        for count in sorted({conn,close-1,close}):
            score,bars=evaluate(chord,iu,ju,length,order,count,true,reference,conn)
            arrays[f'{seed}_{profile}_{method}_M{count}']=bars
            if count==close:
                spans=bool(np.any((bars[:,0]<=B+1e-6)&(bars[:,1]>=1/(3*np.pi)-1e-6)))
        rows.append(dict(seed=seed,profile=profile,method=method,n=n,g=g,B=B,locality_radius_bound=bound,
            r_loc=.18,**{'lambda':1/np.pi},certificate_passes=holds,M_close=close,
            max_chord_through_coverage=float(length[order[:close]].max()),guaranteed_birth_upper=B,
            guaranteed_death_lower=1/(3*np.pi),analytic_lifetime_lower=1/(3*np.pi)-B,
            saved_bar_spans_bound=spans,endpoint_tolerance=1e-6,M_conn=conn,M_second_to_last_consecutive=int(cp[-2])))
    np.savez_compressed(target.with_suffix('.npz'),**arrays)
    dump(target,dict(rows=rows,seconds=time.monotonic()-start))
    return seed,profile,round(time.monotonic()-start,2)
def aggregate():
    records=[json.loads(p.read_text()) for p in sorted((RESULTS/'ellipse_checks/checkpoints').glob('*.json'))]
    rows=table('ellipse_checks/checks.csv',[r for x in records for r in x['rows']])
    arrays={}
    for p in sorted((RESULTS/'ellipse_checks/checkpoints').glob('*.npz')):
        with np.load(p) as z:arrays.update({key:z[key] for key in z.files})
    np.savez_compressed(RESULTS/'ellipse_checks/diagrams.npz',**arrays)
    dump(RESULTS/'ellipse_checks/scope.json',dict(scope='Seed-generated ellipse coordinates and freshly computed connectivity, pre-closure and closure diagrams.',
        samples=len(records),method_sample_cases=len(rows),methods=list(m.METHODS),seeds=CONFIG['ellipse']['seeds'],
        endpoint_tolerance=1e-6,distinct_diagrams=len(arrays),diagram_event_checks=3*len(rows),
        coincident_event_count=3*len(rows)-len(arrays)))
    print('Ellipse checks:',len(records),'samples,',len(rows),'method cases,',len(arrays),'distinct diagrams.',flush=True)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=6);ap.add_argument('--limit',type=int);ap.add_argument('--coordinates-only',action='store_true');args=ap.parse_args()
    coordinates()
    if args.coordinates_only:return
    tasks=[(s,p) for s in CONFIG['ellipse']['seeds'] for p in CONFIG['ellipse']['profiles']]
    if args.limit:tasks=tasks[:args.limit]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i,value in enumerate(ex.map(job,tasks),1):print('Ellipse checks',i,'/',len(tasks),value,flush=True)
    aggregate()
if __name__=='__main__':main()
