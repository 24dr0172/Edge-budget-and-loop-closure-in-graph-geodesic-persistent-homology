"""Check one seeded uniform circle, without saved measurements."""
from pathlib import Path
import argparse,json,time,sys
sys.dont_write_bytecode=True
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path
import graph_methods as m

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--seed',type=int,default=20310101);ap.add_argument('--n',type=int,default=128)
    args=ap.parse_args()
    if args.n<16:ap.error('--n must be at least 16')
    start=time.monotonic();s,xy=m.sample_cloud(args.seed,'uniform',n=args.n)
    chord,iu,ju,length,orders,keys=m.build_orders(xy)
    assert np.diff(np.r_[s,s[0]+1]).max()<1/3,'Choose another seed or increase n: sample lacks the sufficient gap condition.'
    rows=[]
    for method in m.METHODS:
        order=orders[method];conn=m.first_connected(order,iu,ju,args.n);close=m.closure(order,iu,ju,args.n)
        for label,count in [('first_connected',conn),('before_closure',close-1),('closure',close)]:
            use=order[:count];a=iu[use];b=ju[use];w=length[use]
            graph=coo_matrix((np.r_[w,w],(np.r_[a,b],np.r_[b,a])),shape=(args.n,args.n)).tocsr()
            D=shortest_path(graph,directed=False,method='D');assert np.isfinite(D).all()
            bars=m.diagram(D);assert len(bars)==(1 if label=='closure' else 0),(method,label,bars)
            rows.append(dict(method=method,event=label,M=count,positive_H1_bars=len(bars)))
    print(json.dumps(dict(status='passed',seed=args.seed,profile='uniform',n=args.n,seconds=time.monotonic()-start,checks=rows),indent=2))
if __name__=='__main__':main()
