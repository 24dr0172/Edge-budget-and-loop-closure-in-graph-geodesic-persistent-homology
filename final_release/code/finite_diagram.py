"""Independent F2 reduction and finite bottleneck matching for the retained seven-point witness; thresholds are float64, not exact-real chords."""
import itertools
import numpy as np
POSITIVE_TOL=1e-12
def reduction_h1(d):
    """Exact F2 column algebra at double-precision simplex thresholds."""
    n=len(d);simplices=[]
    for size in (1,2,3):
        for v in itertools.combinations(range(n),size):
            t=max((d[i,j] for i,j in itertools.combinations(v,2)),default=0.)
            simplices.append((float(t),size,v))
    simplices.sort();index={v:k for k,(_,_,v) in enumerate(simplices)}
    pivots={};bars=[];tiny=[]
    for k,(death,size,v) in enumerate(simplices):
        col=0
        if size>1:
            for j in range(size):col^=1<<index[v[:j]+v[j+1:]]
        while col:
            pivot=col.bit_length()-1
            if pivot not in pivots:break
            col^=pivots[pivot]
        if col:
            pivot=col.bit_length()-1;pivots[pivot]=col
            birth,bdim,_=simplices[pivot]
            if bdim==2 and death>birth:
                if death-birth>POSITIVE_TOL:bars.append((birth,death))
                else:tiny.append(death-birth)
    return np.asarray(sorted(bars),float).reshape(-1,2),max(tiny,default=0.)

def matching_distance(a,b):
    """Bottleneck matching by exact threshold candidates and augmenting paths.

    Includes diagonal copies on both sides; supports multiple bars without
    assuming the theorem being tested. 'Exact' refers to finite matching,
    not arbitrary-precision evaluation of sine or shortest-path arithmetic.
    """
    na,nb=len(a),len(b);n=na+nb
    if n==0:return 0.
    cost=np.full((n,n),np.inf)
    if na and nb:cost[:na,:nb]=np.max(np.abs(a[:,None,:]-b[None,:,:]),axis=2)
    for i in range(na):cost[i,nb+i]=(a[i,1]-a[i,0])/2
    for j in range(nb):cost[na+j,j]=(b[j,1]-b[j,0])/2
    cost[na:,nb:]=0.
    def feasible(t):
        adj=[np.flatnonzero(row<=t).tolist() for row in cost];owner=[-1]*n
        def add(i,seen):
            for j in adj[i]:
                if seen[j]:continue
                seen[j]=True
                if owner[j]<0 or add(owner[j],seen):owner[j]=i;return True
            return False
        return all(add(i,[False]*n) for i in range(n))
    values=np.unique(cost[np.isfinite(cost)]);lo,hi=0,len(values)-1
    while lo<hi:
        mid=(lo+hi)//2
        if feasible(values[mid]):hi=mid
        else:lo=mid+1
    return float(values[lo])
