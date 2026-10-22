"""Compute the seven-point example with float64 F2 reduction and exact distinctness."""
from fractions import Fraction
import numpy as np
import graph_methods as m
from finite_diagram import reduction_h1
from experiment_core import RESULTS,dump,interval_row,transitions

def main():
    coordinates=['0','0.08','0.15','0.161','0.29','0.60','0.66']
    rational=list(map(Fraction,coordinates))
    distances=[min(abs(a-b),1-abs(a-b)) for i,a in enumerate(rational) for b in rational[i+1:]]
    assert len(distances)==len(set(distances))
    s=np.array(coordinates,float);true=m.intrinsic(s);ref,tiny=reduction_h1(true);assert tiny==0
    points=np.c_[np.cos(2*np.pi*s),np.sin(2*np.pi*s)]/(2*np.pi)
    chord,iu,ju,length,orders,keys=m.build_orders(points)
    # Analytic circle chords keep the thresholds in the float64 calculation.
    chord=np.sin(np.pi*true)/np.pi;length=chord[iu,ju]
    diameter=float(true.max());methods={}
    for method,count in [('R2_global',13),('OR_kNN',12),('mutual_kNN',16)]:
        close=m.closure(orders[method],iu,ju,len(s));rows=[];diagrams=[]
        for M in [count,count+1]:
            D=m.graph_metric(orders[method],iu,ju,length,M,len(s));bars,tiny=reduction_h1(D);assert tiny==0
            rows.append(interval_row(s,ref,diameter,method,M,close,bars));diagrams.append(bars.tolist())
        transition=transitions(rows)[0]
        methods[method]=dict(before=rows[0],after=rows[1],transition=transition,diagrams=diagrams,
            before_edges=[[int(iu[e]),int(ju[e])] for e in orders[method][:count]],
            added_edge=[int(iu[orders[method][count]]),int(ju[orders[method][count]])])
    report=dict(coordinates=coordinates,pairwise_distinct_intrinsic_distances=True,reference_diagram=ref.tolist(),
        analytic_deaths=dict(before=methods['R2_global']['before']['death'],after=methods['R2_global']['after']['death']),
        methods=methods,scope='Rational distinctness; F2 reduction and bottleneck matching with float64 chords.')
    dump(RESULTS/'useful_diagram/rational_counterexample.json',report)
    print('Computed the seven-point example and both matching options independently.')
if __name__=='__main__':main()
