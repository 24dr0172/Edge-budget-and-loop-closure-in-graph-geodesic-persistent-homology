"""Check fundamental-cycle winding in the radial-noise long-interval graphs.

Replays the saved point sets at delta=0.01 and the smallest connected k.
A winding-one graph cycle need not represent a particular persistence interval.
"""
import csv
from pathlib import Path
import numpy as np
import graph_methods as graph

ROOT=Path(__file__).resolve().parents[1]
METHODS=('OR_kNN','mutual_kNN')

def step(a,b):
    return (b-a+0.5)%1.0-0.5

def winding(s,ii,jj,edges):
    n=len(s);parent=list(range(n));adj=[[] for _ in range(n)];extra=[]
    def find(v):
        while parent[v]!=v:
            parent[v]=parent[parent[v]];v=parent[v]
        return v
    for e in edges:
        a,b=int(ii[e]),int(jj[e]);ra,rb=find(a),find(b)
        if ra==rb:extra.append((a,b))
        else:
            parent[ra]=rb;adj[a].append(b);adj[b].append(a)
    potential=np.full(n,np.nan);potential[0]=0.;stack=[0]
    tree_parent=[-1]*n;depth=[0]*n
    while stack:
        a=stack.pop()
        for b in adj[a]:
            if np.isnan(potential[b]):
                potential[b]=potential[a]+step(s[a],s[b]);tree_parent[b]=a;depth[b]=depth[a]+1;stack.append(b)
    assert np.isfinite(potential).all()
    values=np.array([potential[a]+step(s[a],s[b])-potential[b] for a,b in extra])
    integers=np.rint(values)
    assert np.max(np.abs(values-integers),initial=0.)<1e-9
    lengths=[]
    for (a,b),w in zip(extra,integers):
        if abs(w)!=1:continue
        left,right=a,b;distance=0
        while depth[left]>depth[right]:left=tree_parent[left];distance+=1
        while depth[right]>depth[left]:right=tree_parent[right];distance+=1
        while left!=right:left=tree_parent[left];right=tree_parent[right];distance+=2
        lengths.append(distance+1)
    return len(lengths),min(lengths,default=0)

def main():
    with (ROOT/'results/radial_noise_knn.csv').open(newline='') as f:
        expected={(int(r['seed']),r['method']):r for r in csv.DictReader(f)
                  if float(r['delta'])==0.01 and int(r['at_least_half'])==1}
    assert len(expected)==258
    results=[]
    for seed in sorted({k[0] for k in expected}):
        s,clean=graph.sample_cloud(seed,'uniform',n=512)
        normal=clean*(2*np.pi)
        draws=np.random.default_rng(seed+graph.NOISE_SEED_OFFSET).uniform(-1,1,512)
        points=clean+0.01*draws[:,None]*normal
        _,ii,jj,_,orders,ranks=graph.build_orders(points)
        for method in METHODS:
            if (seed,method) not in expected:continue
            row=expected[seed,method];order=orders[method]
            conn=graph.first_connected(order,ii,jj,512)
            k=int(ranks[method][order[conn-1]])
            size=int(np.searchsorted(ranks[method][order],k,side='right'))
            assert k==int(row['k']) and graph.closure(order,ii,jj,512)>size
            count,shortest=winding(s,ii,jj,order[:size])
            assert count>0 and 0<shortest<512
            results.append(dict(seed=seed,method=method,k=k,edges=size,winding_one_fundamental_cycles=count,shortest_winding_one_cycle_vertices=shortest))
        if len(results)%40<2:print('checked',len(results),flush=True)
    path=ROOT/'results/radial_noise_winding.csv'
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
    for method in METHODS:
        group=[r for r in results if r['method']==method]
        print(method,'graphs',len(group),'with winding-one cycle',sum(r['winding_one_fundamental_cycles']>0 for r in group),flush=True)

if __name__=='__main__':main()
