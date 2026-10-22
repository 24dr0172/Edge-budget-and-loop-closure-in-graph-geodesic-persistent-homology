"""Connected bridgeless whole-k layers; sparse reconstruction, no persistence.

The C++ helper contracts two-edge-connected components as edges enter. Its
results are checked against a separate depth-first low-link traversal on
explicit graphs and against the full Cartesian neighbor-rank kernel.
"""
from pathlib import Path
import argparse, hashlib, json, os, resource, struct, subprocess, sys, tempfile, time
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import graph_methods as m
from rank_sweep import fast_inverse, neighboring_keys, wilson

F=Path(__file__).resolve().parents[1]
OUT=F/'results/comparisons/bridgeless'
METHODS=['OR_kNN','mutual_kNN']


def kernel():
    source=F/'code/bridgeless_circle.cpp'
    digest=hashlib.sha256(source.read_bytes()).hexdigest()[:16]
    executable=Path(tempfile.gettempdir())/f'edge-budget-bridgeless-{digest}'
    if not executable.exists():
        subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra',str(source),'-o',str(executable)],check=True)
    return executable


def invoke(sample,cap,edges=None):
    with tempfile.NamedTemporaryFile(prefix='edge-budget-bridge-',suffix='.bin') as stream:
        stream.write(struct.pack('<ii',len(sample),cap))
        if edges is None:
            stream.write(np.asarray(sample,dtype='<f8').tobytes())
        else:
            stream.write(struct.pack('<i',len(edges)))
            stream.write(np.asarray(edges,dtype='<i4').reshape(-1,4).tobytes())
        stream.flush()
        args=[str(kernel())]+([] if edges is None else ['--generic'])+[stream.name]
        result=subprocess.run(args,check=True,capture_output=True,text=True)
    values=result.stdout.split()
    records=[]
    for j,method in enumerate(METHODS):
        conn,br,close,bridges,edge_count=map(int,values[5*j:5*j+5])
        records.append(dict(method=method,k_conn=conn,k_br=br,k_close=close,
                            bridges_at_connection=bridges,edges_at_k_br=edge_count))
    return records,dict(constructed_edges=int(values[10]),construction_seconds=float(values[11]),
                        bridge_seconds=float(values[12]))


def traversal(n,edges):
    """Independent low-link DFS: return number of components and bridges."""
    adj=[[] for _ in range(n)]
    for index,(a,b) in enumerate(edges):
        adj[a].append((b,index));adj[b].append((a,index))
    entered=[-1]*n;low=[0]*n;clock=0;bridges=0
    def dfs(a,previous):
        nonlocal clock,bridges
        entered[a]=low[a]=clock;clock+=1
        for b,index in adj[a]:
            if index==previous:continue
            if entered[b]>=0:low[a]=min(low[a],entered[b])
            else:
                dfs(b,index);low[a]=min(low[a],low[b])
                bridges+=low[b]>entered[a]
    components=0
    for a in range(n):
        if entered[a]<0:components+=1;dfs(a,-1)
    return components,bridges


def full_expected(xy,cap):
    _,iu,ju,_,_,keys=m.build_orders(xy)
    expected=[]
    for method in METHODS:
        conn=br=-1;bc=-1;eb=0
        ring=(ju==iu+1)|((iu==0)&(ju==len(xy)-1))
        close=int(keys[method][ring].max())
        for k in range(1,cap+1):
            use=keys[method]<=k
            edges=list(zip(iu[use],ju[use]))
            components,bridges=traversal(len(xy),edges)
            if conn<0 and components==1:conn=k;bc=int(bridges)
            if components==1 and bridges==0:br=k;eb=len(edges);break
        expected.append(dict(method=method,k_conn=conn,k_br=br,k_close=close,
                             bridges_at_connection=bc,edges_at_k_br=eb))
    return expected


def validate():
    rng=np.random.default_rng(20350119);generic_checks=0;circle_checks=0
    # Disconnected graphs, trees, cycles, ties, and random edge insertions.
    for rep in range(180):
        n=3+rep%21;cap=1+rep%8
        edges=[(a,b,int(rng.integers(1,cap+1)),int(rng.integers(1,cap+1)))
               for a in range(n) for b in range(a+1,n) if rng.uniform()<(.05+.9*(rep%10)/9)]
        got,_=invoke(np.arange(n),cap,edges)
        for j in range(2):
            conn=br=-1;bc=-1;eb=0
            for k in range(1,cap+1):
                chosen=[(a,b) for a,b,x,y in edges if (x if j==0 else y)<=k]
                c,b=traversal(n,chosen)
                if conn<0 and c==1:conn=k;bc=b
                if c==1 and b==0:br=k;eb=len(chosen);break
            assert (got[j]['k_conn'],got[j]['k_br'],got[j]['bridges_at_connection'],got[j]['edges_at_k_br'])==(conn,br,bc,eb)
            generic_checks+=1
    for rep in range(90):
        n=12+rep%61;profile=m.PROFILES[rep%3]
        s,xy=m.sample_cloud(20350101+rep,profile,n=n)
        if np.max(np.diff(np.r_[s,s[0]+1]))>=.5:continue
        keys,_,_=neighboring_keys(s);cap=max(int(a.max()) for a in keys.values())
        got,_=invoke(s,cap)
        assert got==full_expected(xy,cap)
        circle_checks+=2
    report=dict(status='passed',random_graph_method_checks=generic_checks,
                circle_method_checks=circle_checks,
                independent_method='Explicit Cartesian ranks plus independent DFS low-link bridge traversal.',
                tied_layers='All edges of one integer-k layer are inserted before testing.')
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Bridge validation:',report,flush=True)
    return report


def clean():
    saved=pd.read_csv(F/'results/extra_k/samples.csv')
    rows=[];timings=[]
    for (seed,profile),g in saved.groupby(['seed','profile']):
        s,xy=m.sample_cloud(int(seed),profile)
        cap=int(g.k_close.max());start=time.monotonic()
        got,timing=invoke(s,cap)
        # All 60 clean-study samples also get a fully independent graph check.
        expected=full_expected(xy,cap)
        assert got==expected
        for record in got:
            prior=g[g.method==record['method']].iloc[0]
            assert record['k_conn']==prior.k_conn and record['k_close']==prior.k_close
            rows.append(dict(seed=int(seed),profile=profile,n=len(s),**record,
                             closed_at_connection=record['k_conn']==record['k_close'],
                             closed_at_bridgeless=record['k_br']==record['k_close'],
                             closed_at_plus3=record['k_conn']+3>=record['k_close']))
        timings.append(dict(seed=int(seed),profile=profile,n=len(s),seconds=time.monotonic()-start,**timing))
    pd.DataFrame(rows).to_csv(OUT/'clean_samples.csv',index=False)
    pd.DataFrame(timings).to_csv(OUT/'clean_timings.csv',index=False)
    summarize('clean')
    print('Clean study: 60 samples, 120 method records, full Cartesian/DFS replay passed.',flush=True)


def sweep(max_n=1000000,pilot=False):
    source=pd.read_csv(F/'results/extra_k/rank_sweep/samples.csv')
    source=source[source.n<=max_n]
    grouped=list(source.groupby(['n','seed','profile'],sort=True))
    if pilot:
        grouped=[(key,g) for key,g in grouped if key[0] in [512,4096,32768]]
        picked={}
        for key,g in grouped:picked.setdefault((key[0],key[2]),(key,g))
        grouped=list(picked.values())
    output=OUT/('pilot_samples.csv' if pilot else 'sweep_samples.csv')
    tfile=OUT/('pilot_timings.csv' if pilot else 'sweep_timings.csv')
    rows=pd.read_csv(output).to_dict('records') if output.exists() else []
    timings=pd.read_csv(tfile).to_dict('records') if tfile.exists() else []
    completed={(int(r['n']),int(r['seed']),r['profile']) for r in timings}
    start=time.monotonic()
    for index,(key,g) in enumerate(grouped):
        n,seed,profile=key
        if key in completed:continue
        t=time.monotonic()
        uniforms=np.sort(np.random.default_rng(int(seed)).uniform(size=int(n)))
        phase=float(np.random.default_rng(int(seed)+m.PHASE_SEED_OFFSET).uniform())
        s=fast_inverse(uniforms,profile,phase)
        cap=int(g.k_close.max());got,timing=invoke(s,cap)
        for record in got:
            prior=g[g.method==record['method']].iloc[0]
            assert record['k_conn']==prior.k_conn and record['k_close']==prior.k_close
            assert record['k_conn']<=record['k_br']<=record['k_close']
            rows.append(dict(seed=int(seed),profile=profile,n=int(n),**record,
                             closed_at_connection=record['k_conn']==record['k_close'],
                             closed_at_bridgeless=record['k_br']==record['k_close'],
                             closed_at_plus3=record['k_conn']+3>=record['k_close']))
        timings.append(dict(seed=int(seed),profile=profile,n=int(n),seconds=time.monotonic()-t,
                            peak_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,**timing))
        if pilot or (index+1)%100==0 or index+1==len(grouped):
            pd.DataFrame(rows).to_csv(output,index=False);pd.DataFrame(timings).to_csv(tfile,index=False)
            print(f'Bridge {"pilot" if pilot else "sweep"} {index+1}/{len(grouped)}: n={n}, {profile}, {timings[-1]["seconds"]:.4f}s, elapsed {time.monotonic()-start:.1f}s',flush=True)
    if not pilot:summarize('sweep')


def summarize(study):
    raw=pd.read_csv(OUT/f'{study}_samples.csv');rows=[]
    for key,g in raw.groupby(['n','profile','method']):
        record=dict(zip(['n','profile','method'],key));record['samples']=len(g)
        for label,column in [('connection','closed_at_connection'),('bridgeless','closed_at_bridgeless'),('plus3','closed_at_plus3')]:
            successes=int(g[column].sum());lo,hi=wilson(successes,len(g))
            record.update({label+'_closed':successes,label+'_rate':successes/len(g),
                           label+'_wilson_lo':lo,label+'_wilson_hi':hi})
        record['median_extra_k_br']=float(np.median(g.k_br-g.k_conn))
        record['median_remaining_k']=float(np.median(g.k_close-g.k_br))
        rows.append(record)
    pd.DataFrame(rows).to_csv(OUT/f'{study}_summary.csv',index=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--validate',action='store_true')
    parser.add_argument('--clean',action='store_true')
    parser.add_argument('--pilot',action='store_true')
    parser.add_argument('--sweep',action='store_true')
    parser.add_argument('--max-n',type=int,default=1000000)
    args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if args.validate:validate()
    if args.clean:clean()
    if args.pilot:sweep(args.max_n,True)
    if args.sweep:sweep(args.max_n)
