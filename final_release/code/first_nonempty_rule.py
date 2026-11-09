"""Distance-only bisection under Corollary 2's circle hypotheses.

The inputs are points and a prescribed graph family. Curve order, arc
coordinates and the intrinsic reference are not used. The global circle
monotonicity is a mathematical precondition, not a test performed here.
Unordered noise-free circle samples are valid inputs under these hypotheses.
The theoretical bisection bound uses endpoint statuses known from the theory.
This implementation also evaluates the first connected graph and, if empty,
the complete graph to validate the bracket. The returned persistence_queries
includes these endpoint tests, including the initial whole-layer test.
For a general curve, use only a separately certified search bracket.
"""
import numpy as np
from graph_methods import build_orders,first_connected,graph_metric
from first_nonempty_check import query

def first_nonempty_circle(points,method='OR_kNN',whole_layers=False):
    if method not in ('R2_global','OR_kNN','mutual_kNN'):
        raise ValueError('Choose R2_global, OR_kNN or mutual_kNN.')
    points=np.asarray(points,float)
    if points.ndim!=2 or len(points)<3 or not np.isfinite(points).all():
        raise ValueError('Supply a finite n-by-d point array with n >= 3.')
    C,iu,ju,w,orders,keys=build_orders(points);order=orders[method];n=len(points)
    conn=first_connected(order,iu,ju,n);N=len(order)
    if whole_layers and method in keys:
        kk=keys[method][order]
        budgets=np.flatnonzero(np.r_[kk[1:]!=kk[:-1],True])+1
    else:budgets=np.arange(1,N+1)
    first=int(np.searchsorted(budgets,conn));cache={};ceilings=np.maximum.accumulate(w[order])
    def positive(index):
        m=int(budgets[index])
        if m not in cache:
            D=C if m==N else graph_metric(order,iu,ju,w,m,n)
            cache[m]=bool(len(query(D,float(ceilings[m-1]))[0]))
        return cache[m]
    if positive(first):answer=first
    else:
        lo,hi=first,len(budgets)-1
        if not positive(hi):
            raise ValueError('The upper graph has an empty diagram; the circle guarantee does not apply.')
        while hi-lo>1:
            mid=(lo+hi)//2
            if positive(mid):hi=mid
            else:lo=mid
        answer=hi
    m=int(budgets[answer]);k=int(keys[method][order[m-1]])if whole_layers and method in keys else None
    return dict(method=method,whole_layers=bool(whole_layers),M=m,k=k,
                persistence_queries=len(cache),M_conn=conn,
                scope='Correct under the circle hypotheses of Corollary 2; global monotonicity is not certified from these inputs.')
