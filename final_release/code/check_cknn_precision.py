"""Certify the CkNN witness with exact edge scores and rigorous cost intervals."""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal, localcontext
from math import isqrt
import argparse, hashlib, json
import sys
sys.dont_write_bytecode=True
import numpy as np
import graph_methods as m

def circle_costs(s, digits):
    """Independent high-precision circle evaluation from the saved arc parameters."""
    with localcontext() as ctx:
        ctx.prec=digits+15
        def atan(z):
            total=term=z;k=1
            while abs(term)>Decimal(10)**(-digits-8):
                term*=-z*z
                total+=term/(2*k+1);k+=1
            return total
        pi=16*atan(Decimal(1)/5)-4*atan(Decimal(1)/239)
        def chord(a,b):
            d=abs(Decimal.from_float(float(s[a]))-Decimal.from_float(float(s[b])))
            d=min(d,1-d);x=pi*d;term=total=x;k=1
            while abs(term)>Decimal(10)**(-digits-8):
                term*=-(x*x)/(2*k*(2*k+1));total+=term;k+=1
            return total/pi
        outer=chord(190,201);left=sum((chord(i,i+1) for i in range(190,200)),Decimal(0))
        return dict(decimal_digits=digits,outer=str(+outer),left_chain=str(+left),
                    excess=str(+(left-outer)))

def check(root):
    root=Path(root)
    inputs=root/'results/clean_closure/latent_coordinates.npz'
    with np.load(inputs) as z:s=z['20310105_localized_trough']
    replay,xy=m.sample_cloud(20310105,'localized_trough')
    assert np.array_equal(s,replay)
    n=len(xy);ratios=[float(v).as_integer_ratio() for v in xy.flat]
    denominator=max(b for a,b in ratios)
    points=[(ratios[2*i][0]*(denominator//ratios[2*i][1]),
             ratios[2*i+1][0]*(denominator//ratios[2*i+1][1])) for i in range(n)]
    squared=[[0]*n for _ in range(n)];pairs=[]
    for a in range(n):
        for b in range(a+1,n):
            length2=sum((points[a][j]-points[b][j])**2 for j in [0,1])
            squared[a][b]=squared[b][a]=length2;pairs.append((a,b))
    rho2=[sorted(squared[a][b] for b in range(n) if a!=b)[12] for a in range(n)]
    # Fourth powers remove both square roots from the ordering score:
    # (c(a,b)/sqrt(rho(a)*rho(b)))^4 = c(a,b)^4/(rho(a)^2*rho(b)^2).
    order=sorted(pairs,key=lambda e:(Fraction(squared[e[0]][e[1]]**2,
                                    rho2[e[0]]*rho2[e[1]]),*e))
    pos={pair:i+1 for i,pair in enumerate(order)}
    getpos=lambda a,b:pos[tuple(sorted((a,b)))]
    assert getpos(190,201)==5840
    close=max(getpos(a,(a+1)%n) for a in range(n));assert close==6255
    assert all(getpos(a,a+1)<5840 for a in range(190,201))
    assert getpos(190,200)>5840
    digits=100;scale=10**digits
    def bounds(a,b):
        value=squared[a][b]*scale*scale;lower=isqrt(value)
        return lower,lower+(lower*lower!=value)
    outer=bounds(190,201)
    chain=[bounds(a,a+1) for a in range(190,200)]
    left=(sum(a for a,b in chain),sum(b for a,b in chain))
    excess=(left[0]-outer[1],left[1]-outer[0]);assert excess[0]>0
    with localcontext() as ctx:
        ctx.prec=digits+30
        todecimal=lambda value:str(Decimal(value)/Decimal(scale*denominator))
        interval=lambda x:[todecimal(x[0]),todecimal(x[1])]
        result=dict(status='passed',seed=20310105,profile='localized_trough',n=n,k0=13,
            input_sha256=hashlib.sha256(inputs.read_bytes()).hexdigest(),
            input_coordinates_sha256=hashlib.sha256(xy.tobytes()).hexdigest(),
            input_convention='Float64 Cartesian sample reconstructed from saved arc parameters, then treated as exact dyadic rationals.',
            ordering='Exact rational fourth powers of CkNN scores; ties by endpoint labels.',
            cost_bounds='Integer square roots, scaled by 10^100; outward rational bounds.',
            common_coordinate_denominator=str(denominator),decimal_places=digits,
            edge=[190,201],edge_position=getpos(190,201),closure_position=close,
            interior_vertex=200,left_direct_edge_position=getpos(190,200),
            consecutive_route_available=True,outer_cost_interval=interval(outer),
            left_route_cost_interval=interval(left),excess_interval=interval(excess),
            interval_width_upper_bound=str(Decimal(11)/(scale*denominator)),
            ideal_circle_evaluations=[circle_costs(s,p) for p in [60,100]])
    for row in result['ideal_circle_evaluations']:
        assert Decimal(row['excess'])>Decimal('0.0000058')
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--output',type=Path);args=ap.parse_args();report=check(args.release)
    target=args.output or args.release/'results/order_checks/precision_check.json'
    target.write_text(json.dumps(report,indent=2)+'\n')
    print('Exact CkNN check passed; excess interval:',report['excess_interval'],flush=True)
