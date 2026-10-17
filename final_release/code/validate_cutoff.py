"""Validate the graph-metric H1 cutoff against unrestricted Ripser diagrams."""
from pathlib import Path
import json,sys,time
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import shortest_path
from run_studies import finite_graph_diagram,m,PAPER,SUPPORT

def canonical(bars):
    bars=np.asarray(bars).reshape(-1,2)
    return bars[np.lexsort((bars[:,1],bars[:,0]))] if len(bars) else bars

def same(a,b):
    a,b=canonical(a),canonical(b)
    assert a.shape==b.shape and np.array_equal(a,b),(a,b)

def main():
    start=time.monotonic();rng=np.random.default_rng(2026101501)
    for i in range(150):
        n=int(rng.integers(5,26));w=np.full((n,n),np.inf);np.fill_diagonal(w,0.)
        for a in range(n):w[a,(a+1)%n]=w[(a+1)%n,a]=rng.uniform(.001,2.)
        for a in range(n):
            for b in range(a+1,n):
                if rng.uniform()<.15:w[a,b]=w[b,a]=rng.uniform(.001,2.)
        D=shortest_path(w,directed=False)
        same(finite_graph_diagram(D,float(w[np.isfinite(w)].max())),m.diagram(D))
    checked=0
    for model,delta in [('radial',0.),('radial',.01),('isotropic',.01)]:
        stem=f'{model}_20320101_delta{delta:g}'
        path=SUPPORT/f'results/noise/checkpoints/{stem}.csv';d=pd.read_csv(path)
        with np.load(path.with_suffix('.npz')) as saved:
            points=saved['observed'];chord,iu,ju,length,orders,keys=m.build_orders(points)
            for method in m.METHODS:
                g=d[(d.method==method)&d.connected].sort_values('M')
                # Every exact event, plus alternate grid points, including both extremes.
                chosen=pd.concat([g.iloc[::2],g.iloc[[-1]],g[g.first_connected|g.before_closure|g.closure]]).drop_duplicates('M')
                for row in chosen.itertuples():
                    D=m.graph_metric(orders[method],iu,ju,length,row.M,len(points))
                    maximum=float(length[orders[method][:row.M]].max())
                    same(finite_graph_diagram(D,maximum),saved[f'{method}_M{row.M}'])
                    checked+=1
    report=dict(status='passed',random_positive_weight_graphs=150,n512_saved_unrestricted_diagrams=checked,
                comparison='Exact equality of lexicographically sorted float32 endpoint pairs',elapsed_seconds=time.monotonic()-start,
                justification='Every Rips edge is homologous to its generating shortest graph path. Above the largest graph edge no positive H1 class is born; an H1-empty cutoff therefore completes the barcode. The last cutoff has a universal vertex.')
    (SUPPORT/'results/cutoff_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
