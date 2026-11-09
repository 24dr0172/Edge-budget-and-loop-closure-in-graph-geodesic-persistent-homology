"""Certify a window on the existing ellipse using exact rational interval bounds.

No random samples, numerical locality search, graph geodesics or persistence
computations are used. The input fixture fixes the angles of one existing cloud.
"""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,csv
ROOT=Path(__file__).resolve().parents[1]
SCALE=10**40

def sincos(x,cos=False):
    """Enclose sin/cos of a rational x using a Taylor sum and alternating tail."""
    assert 0<=x<7
    term=F(1) if cos else x
    total=term
    for k in range(1,41):
        den=(2*k-1)*(2*k) if cos else (2*k)*(2*k+1)
        term=-term*x*x/den;total+=term
    den=81*82 if cos else 82*83
    remainder=abs(term*x*x/den)
    # All remaining term magnitudes decrease; the alternating tail is no
    # larger than the first omitted term. Rational arithmetic is exact.
    lo=(total-remainder)*SCALE;hi=(total+remainder)*SCALE
    return lo.numerator//lo.denominator,-((-hi.numerator)//hi.denominator)

def sq_interval(lo,hi):
    return (0 if lo<=0<=hi else min(lo*lo,hi*hi),max(lo*lo,hi*hi))

def verify(root=ROOT):
    root=Path(root);folder=root/'results/certified_ellipse_window'
    fixture=folder/'angles.json';data=json.loads(fixture.read_text())
    angles=[F.from_float(float.fromhex(v)) for v in data['angles_hex']]
    n=len(angles);assert n==data['n']==512
    assert all(a<b for a,b in zip(angles,angles[1:]))
    gaps=[b-a for a,b in zip(angles,angles[1:])]+[angles[0]+F(44,7)-angles[-1]]
    assert all(0<g<F(111,1000) for g in gaps) # 2*pi < 44/7
    # L >= 4*sqrt(5), speed <= 2, so g < .111/(2*sqrt(5)) < 1/40.
    assert F(111,1000)**2<F(1,20)**2*5
    coords=[]
    for t in angles:
        a,b=sincos(t,cos=True);c,d=sincos(t);coords.append((2*a,2*b,c,d))
    edges=[];intervals=[]
    for i in range(n):
        for j in range(i+1,n):
            x0,x1,y0,y1=coords[i];u0,u1,v0,v1=coords[j]
            lo,hi=sq_interval(x0-u1,x1-u0);a,b=sq_interval(y0-v1,y1-v0)
            edges.append((i,j));intervals.append((lo+a,hi+b))
    exact_order=sorted(range(len(edges)),key=lambda e:sum(intervals[e]))
    assert all(intervals[e][0]>0 for e in exact_order)
    margins=[intervals[b][0]-intervals[a][1] for a,b in zip(exact_order,exact_order[1:])]
    assert min(margins)>0 # Global distinctness and every comparison are certified.
    # Raw locality boundary .5; in length-one units this is r_loc=1/(2L).
    cutoff=SCALE*SCALE//4
    assert all(hi<cutoff or lo>cutoff for lo,hi in intervals)
    nonlocal_edges={e for e,(lo,hi) in enumerate(intervals) if lo>cutoff}
    ring={e for e,(i,j) in enumerate(edges) if j==i+1 or (i==0 and j==n-1)}
    ranks=[0]*n;keys_or=[0]*len(edges);keys_mutual=[0]*len(edges)
    for e in exact_order:
        i,j=edges[e];ranks[i]+=1;ranks[j]+=1
        keys_or[e]=min(ranks[i],ranks[j]);keys_mutual[e]=max(ranks[i],ranks[j])
    ordinal={e:k for k,e in enumerate(exact_order)}
    expected={'radius':(4236,13562,None,None),'OR':(2117,10960,8,40),'mutual':(2599,10756,12,45)}
    results={}
    for method,keys in [('radius',None),('OR',keys_or),('mutual',keys_mutual)]:
        order=exact_order if keys is None else sorted(exact_order,key=lambda e:(keys[e],ordinal[e]))
        pos={e:k+1 for k,e in enumerate(order)}
        close=max(pos[e] for e in ring);geom=min(pos[e] for e in nonlocal_edges)
        kc=None if keys is None else max(keys[e] for e in ring)
        kg=None if keys is None else min(keys[e] for e in nonlocal_edges)
        assert (close,geom,kc,kg)==expected[method]
        assert close<geom and (kc is None or kc<kg)
        results[method]=dict(M_close=close,M_geom=geom,M_end=geom,k_close=kc,k_geom=kg,k_end=kg)
    with (root/'results/noncircular_locality/locality_check_results.csv').open() as f:
        original=next(row for row in csv.DictReader(f) if row['curve']=='ellipse' and int(row['seed'])==data['seed'])
    assert all(int(original[method+'_M_close'])==v['M_close'] for method,v in results.items())
    # Analytic bound: g<1/40; lambda >=3/pi; b*<=1/3+2g; Delta>=1/2-g.
    # pi<22/7 gives 3/pi>21/22, and (46/57)*(1-21/22)<.037.
    error_upper=F(46,57)*(1-F(21,22));assert error_upper<F(37,1000)
    # B<=g<1/40<1/pi<=lambda*b*, hence survival cannot fail before M_geom.
    assert F(1,40)<F(7,22)
    return dict(status='passed',scope='One exact-angle realization of the first existing n=512 ellipse cloud; no new cohort or diagrams.',seed=data['seed'],n=n,fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),squared_chord_intervals=len(edges),strict_adjacent_order_comparisons=len(margins),minimum_squared_order_separation_lower=str(F(min(margins),SCALE*SCALE)),angle_gap_upper='111/1000 radians',intrinsic_gap_upper='1/40',raw_ellipse_reach='1/2',raw_locality_radius='1/2',edgewise_ratio_lower='3/pi',normalized_error_upper='(46/57)*(1-3/pi) < 0.037',perimeter_lower='4*sqrt(5)',windows=results,numerical_locality_estimate_used=False)

if __name__=='__main__':
    report=verify();dest=ROOT/'results/certified_ellipse_window/certificate.json'
    dest.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
