"""Independent F2 persistence reduction for finite metric spaces."""
import itertools
import numpy as np
# Tolerance for this independent witness helper only; study diagrams retain all positive Ripser bars.
POSITIVE_TOL=1e-12
def reduction_h1(d):
    """Reduce over F2; return bars above POSITIVE_TOL and the largest omitted positive lifetime."""
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
