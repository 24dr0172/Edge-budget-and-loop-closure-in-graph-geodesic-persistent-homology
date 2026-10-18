"""Generate the matched radial-noise diagnostic, coordinates and diagrams from seeds."""
import argparse,concurrent.futures,time,json
import numpy as np
import pandas as pd
import paper1_methods as m
from run_studies import measurements
from experiment_core import RESULTS,CONFIG,table,dump

def job(seed):
    start=time.monotonic();out=RESULTS/'matched_noise/checkpoints';out.mkdir(parents=True,exist_ok=True)
    target=out/f'{seed}.json'
    if target.exists():return seed,'resumed'
    s,clean=m.sample_cloud(seed,'uniform',n=CONFIG['n']);n=len(s)
    true=m.intrinsic(s);reference=m.diagram(true);diameter=float(true.max())
    life=float(np.max(reference[:,1]-reference[:,0]))
    normal=clean*(2*np.pi);draw=np.random.default_rng(seed+m.NOISE_SEED_OFFSET).uniform(-1,1,n)
    arrays={f'{seed}_latent':s,f'{seed}_clean':clean,f'{seed}_normal':normal,f'{seed}_signed_draw':draw,
            f'oracle_{seed}':reference}
    rows=[];closing=[]
    for delta in CONFIG['noise']['amplitudes']:
        points=clean+delta*draw[:,None]*normal
        arrays[f'{seed}_delta{delta:g}_observed']=points
        chord,iu,ju,length,orders,keys=m.build_orders(points)
        for method in ('R2_global','OR_kNN'):
            order=orders[method];conn=m.first_connected(order,iu,ju,n);close=m.closure(order,iu,ju,n)
            events=[('first_connected',conn)]
            if delta==.01:events.extend([('pre_coverage',close-1),('coverage',close)])
            cache={}
            for event,count in events:
                if count not in cache:
                    D=m.graph_metric(order,iu,ju,length,count,n)
                    cache[count]=measurements(D,true,reference,maximum_edge=float(length[order[:count]].max()))
                score,bars=cache[count]
                arrays[f'{seed}_delta{delta:g}_{method}_{event}']=bars
                row=dict(seed=seed,n=n,delta=delta,delta_n=delta*n,method=method,event=event,M=count,M_conn=conn,
                    oracle_diameter=diameter,oracle_longest_lifetime=life,connected=score['connected'],n_components=1,
                    empty_H1=len(bars)==0,H1_bars=len(bars),longest_over_reference=score['longest_over_reference'],
                    e_inf=score['e_inf'],e_2=score['e_2'],e_PH=score['e_PH'],stability_ok=score['e_PH']<=score['e_inf']+5e-8/diameter)
                (rows if event=='first_connected' else closing).append(row if event=='first_connected' else dict(row,M_cov=close))
    np.savez_compressed(target.with_suffix('.npz'),**arrays)
    dump(target,dict(rows=rows,closing=closing,seconds=time.monotonic()-start))
    return seed,round(time.monotonic()-start,2)

def aggregate():
    out=RESULTS/'matched_noise';(out/'diagrams').mkdir(exist_ok=True)
    rows=[];closing=[];closure_arrays={};coordinates={}
    for seed in CONFIG['noise']['seeds']:
        record=json.loads((out/f'checkpoints/{seed}.json').read_text());rows.extend(record['rows']);closing.extend(record['closing'])
        arrays={}
        with np.load(out/f'checkpoints/{seed}.npz') as z:
            for name in z.files:
                if name.startswith('oracle_'):
                    arrays[name]=z[name];closure_arrays[name]=z[name]
                elif name.endswith('_first_connected'):arrays[name]=z[name]
                elif name.endswith(('_pre_coverage','_coverage')):closure_arrays[name]=z[name]
                else:coordinates[name]=z[name]
        np.savez_compressed(out/f'diagrams/diagrams_{seed}.npz',**arrays)
    raw=table('matched_noise/cells.csv',rows);table('matched_noise/closure_comparison.csv',closing)
    summary=[]
    for (delta,method,event),g in raw.groupby(['delta','method','event']):
        summary.append(dict(delta=delta,delta_n=delta*CONFIG['n'],method=method,event=event,attempted=len(g),
            connected=int(g.connected.sum()),empty_H1=int(g.empty_H1.sum()),median_longest_over_reference=float(g.longest_over_reference.median()),
            median_M=float(g.M.median()),long_at_half=int((g.longest_over_reference>=.5).sum())))
    table('matched_noise/summary.csv',summary)
    np.savez_compressed(out/'closure_diagrams.npz',**closure_arrays)
    np.savez_compressed(out/'coordinates.npz',**coordinates)
    dump(out/'config.json',dict(n=CONFIG['n'],seeds=CONFIG['noise']['seeds'],delta=CONFIG['noise']['amplitudes'],
        methods=['R2_global','OR_kNN'],events=['first_connected'],noise_model='Uniform radial displacement; the same signed draw is reused across amplitudes.',
        arithmetic='clean + delta * signed_draw[:,None] * normal; normal=clean*(2*pi).',
        source='Fresh seed-generated coordinates and persistence diagrams.'))
    print('Generated 240 first-connectivity rows and 80 before/at-closure rows from seeds.',flush=True)
def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=6);p.add_argument('--aggregate-only',action='store_true');a=p.parse_args()
    if not a.aggregate_only:
        with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
            for result in pool.map(job,CONFIG['noise']['seeds']):print('Radial diagnostic',result,flush=True)
    aggregate()
if __name__=='__main__':main()
