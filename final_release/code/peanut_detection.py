"""Replay the peanut cycle-closure and first-nonempty-prefix measurements.

The saved samples use the same arc-coordinate draws at each shape parameter.
The prefix search tests connected graphs in order and does not assume that
nonemptiness is monotone. --closure-only skips that expensive search.
"""
from pathlib import Path
from functools import lru_cache
import argparse
import csv
import sys
import numpy as np
from ripser import ripser
import graph_methods as graph

ROOT=Path(__file__).resolve().parents[1]
METHODS=('R2_global','OR_kNN','mutual_kNN')
ALPHAS=(0.,.4,.6,.8,.9)
N=512
GRID=400001


@lru_cache(maxsize=5)
def arc_table(alpha):
    theta=np.linspace(0.,2*np.pi,GRID)
    radius=1+alpha*np.cos(2*theta)
    xy=np.c_[radius*np.cos(theta),radius*np.sin(theta)]
    return theta,np.r_[0.,np.cumsum(np.linalg.norm(np.diff(xy,axis=0),axis=1))]


def sample(seed,alpha):
    s=np.sort(np.random.default_rng(seed).random(N))
    theta,cumulative=arc_table(alpha)
    length=float(cumulative[-1]); angles=np.interp(s*length,cumulative,theta)
    radius=1+alpha*np.cos(2*angles)
    return s,np.c_[radius*np.cos(angles),radius*np.sin(angles)]/length,length


def query(D,edge_ceiling=None):
    cap=np.inf if edge_ceiling is None else float(np.nextafter(np.float32(edge_ceiling),np.float32(np.inf)))
    raw=np.asarray(ripser(D,distance_matrix=True,maxdim=1,thresh=cap)['dgms'][1],float).reshape(-1,2)
    return raw[raw[:,1]>raw[:,0]]


def update_metric(D,a,b,cost):
    da,db=D[:,a].copy(),D[:,b].copy()
    np.minimum(D,da[:,None]+cost+db[None,:],out=D)
    np.minimum(D,db[:,None]+cost+da[None,:],out=D)
    np.fill_diagonal(D,0.)


def select(order,ii,jj,weights,n,conn):
    D=graph.graph_metric(order,ii,jj,weights,conn,n)
    ceiling=float(weights[order[:conn]].max())
    for M in range(conn,len(order)+1):
        if M>conn:
            e=int(order[M-1]);update_metric(D,int(ii[e]),int(jj[e]),float(weights[e]))
            ceiling=max(ceiling,float(weights[e]))
        if (M-conn)%64==0 and M>conn:
            fresh=graph.graph_metric(order,ii,jj,weights,M,n)
            assert np.array_equal(fresh.astype(np.float32),D.astype(np.float32))
            D=fresh
        if not len(query(D,ceiling)):continue
        fresh=graph.graph_metric(order,ii,jj,weights,M,n)
        assert np.array_equal(fresh.astype(np.float32),D.astype(np.float32))
        bars=query(fresh)
        assert len(bars)>0
        if M>conn:
            assert len(query(graph.graph_metric(order,ii,jj,weights,M-1,n)))==0
        return M,bars
    raise AssertionError('No detecting graph was found')


def arc_support(D,chords,ii,jj,chosen,tol=0.):
    n=len(D);ring=chords[np.arange(n),np.roll(np.arange(n),-1)]
    valid=np.zeros((n,n),bool)
    for a in range(n):
        for direction in (1,-1):
            arc=(a+direction*np.arange(n))%n
            distances=np.maximum.accumulate(D[a,arc])
            steps=ring[arc[:-1]] if direction==1 else ring[arc[1:]]
            costs=np.r_[0.,np.maximum.accumulate(steps)]
            valid[a,arc]|=(distances<=chords[a,arc]+tol)&(costs<=chords[a,arc]+tol)
    supported=valid[ii[chosen],jj[chosen]]|valid[jj[chosen],ii[chosen]]
    return bool(supported.all())


def read_rows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))


def replay(seed,alpha,method,closure_only=False):
    close_rows=read_rows(ROOT/'results/peanut_neck/trials.csv')
    first_rows=read_rows(ROOT/'results/first_nonempty/trials.csv')
    pick=lambda rows:next(r for r in rows if int(r['seed'])==seed and float(r['alpha'])==alpha and r['method']==method)
    record,first=pick(close_rows),pick(first_rows)
    s,xy,_=sample(seed,alpha)
    filename=f'a{alpha:g}_seed{seed}.npz'
    with np.load(ROOT/'results/peanut_neck/checkpoints'/filename,allow_pickle=False) as saved:
        assert np.array_equal(s,saved['s']) and np.allclose(xy,saved['xy'],rtol=0,atol=1e-15)
        expected_closure=saved[f'{method}_closure_bars'].copy()
    chords,ii,jj,weights,orders,_=graph.build_orders(xy);order=orders[method]
    conn=graph.first_connected(order,ii,jj,N);close=graph.closure(order,ii,jj,N)
    assert (conn,close)==(int(record['M_conn']),int(record['M_close']))
    assert bool(graph.certificate(s,chords,order,ii,jj)[1])==(record['short_hops_first']=='True')
    D=graph.graph_metric(order,ii,jj,weights,close,N)
    closure_bars=query(D)
    assert closure_bars.shape==expected_closure.shape
    assert np.allclose(closure_bars,expected_closure,rtol=0,atol=5e-8)
    assert arc_support(D,chords,ii,jj,order[:close])==(record['closure_supported']=='True')
    if not closure_only:
        M,bars=select(order,ii,jj,weights,N,conn)
        assert M==int(first['M_first']) and len(bars)==int(first['first_bars'])
        fname=f'a{alpha:g}_seed{seed}_{method}.npz'
        with np.load(ROOT/'results/first_nonempty/checkpoints'/fname,allow_pickle=False) as z:
            assert bars.shape==z['first_bars'].shape
            assert np.allclose(bars,z['first_bars'],rtol=0,atol=5e-8)
            if method=='R2_global':
                m2=int(np.searchsorted(weights[order],2*weights[order[conn-1]],side='right'))
                double=query(graph.graph_metric(order,ii,jj,weights,m2,N))
                assert double.shape==z['double_bars'].shape
                assert m2==int(first['M_double']) and np.allclose(double,z['double_bars'],rtol=0,atol=5e-8)
    print(f'alpha={alpha:g}, seed={seed}, {method}: passed',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=int,default=20361021)
    p.add_argument('--alpha',type=float,choices=ALPHAS,default=.4)
    p.add_argument('--method',choices=METHODS,default='R2_global')
    p.add_argument('--all',action='store_true',help='Replay all 300 saved method/sample cases')
    p.add_argument('--closure-only',action='store_true')
    a=p.parse_args()
    jobs=[(seed,alpha,method) for alpha in ALPHAS for seed in range(20361021,20361041) for method in METHODS] if a.all else [(a.seed,a.alpha,a.method)]
    for job in jobs:replay(*job,closure_only=a.closure_only)
    print(f'Peanut replay passed: {len(jobs)} cases.')

if __name__=='__main__':main()
