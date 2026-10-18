"""Recompute arc support and the exact common-order obstruction of Example 1."""
from pathlib import Path
import itertools,json,sys
sys.dont_write_bytecode=True
import numpy as np
from scipy.sparse.csgraph import shortest_path
PAPER=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PAPER/'code'))
from verify_arc_support import supported

def rank_mod(matrix,p):
    a=np.array(matrix,int).copy()%p;rank=0
    for j in range(a.shape[1]):
        choices=np.flatnonzero(a[rank:,j])
        if not len(choices):continue
        k=rank+int(choices[0]);a[[rank,k]]=a[[k,rank]]
        a[rank]=a[rank]*pow(int(a[rank,j]),-1,p)%p
        for k in range(rank+1,len(a)):a[k]=(a[k]-a[k,j]*a[rank])%p
        rank+=1
        if rank==len(a):break
    return rank

def exact_metric():
    """Shortest paths in the Q-basis (1,sqrt(2),sqrt(5),sqrt(10))/5.

    Integer square-root bounds certify every comparison. Equal coefficient
    vectors identify tied thresholds without a floating-point tolerance.
    """
    import math
    from fractions import Fraction
    scale=10**80
    roots=[scale]+[math.isqrt(k*scale*scale) for k in [2,5,10]]
    def interval(v):
        lo=sum(a*b for a,b in zip(v,roots))
        hi=lo+sum(v[1:])
        return lo,hi
    def smaller(a,b):
        if a==b:return False
        al,ah=interval(a);bl,bh=interval(b)
        assert ah<bl or bh<al, 'Increase precision for overlapping distinct radical forms'
        return ah<bl
    zero=(0,0,0,0)
    D=[[zero if a==b else None for b in range(6)] for a in range(6)]
    weights={(0,1):(0,1,0,0),(1,2):(6,0,0,0),
             (2,3):(0,0,2,0),(3,4):(0,0,0,1),
             (4,5):(0,5,0,0),(0,5):(0,5,0,0),
             (0,2):(0,5,0,0),(0,3):(0,0,0,3)}
    for (a,b),v in weights.items():D[a][b]=D[b][a]=v
    for k in range(6):
        for a in range(6):
            for b in range(6):
                if D[a][k] is None or D[k][b] is None:continue
                candidate=tuple(x+y for x,y in zip(D[a][k],D[k][b]))
                if D[a][b] is None or smaller(candidate,D[a][b]):D[a][b]=candidate
    def bounds(v):
        lo,hi=interval(v)
        return dict(lower_numerator=str(lo),upper_numerator=str(hi),
                    denominator=str(5*scale))
    levels=sorted({D[a][b] for a in range(6) for b in range(a+1,6)},
                  key=lambda v:interval(v)[0])
    assert all(interval(a)[1]<interval(b)[0] for a,b in zip(levels,levels[1:]))
    numeric=np.array([[float(Fraction(sum(x*y for x,y in zip(v,roots)),5*scale))
                       for v in row] for row in D])
    return D,numeric,[zero]+levels,bounds


def cyclic_orientation(order,edges):
    """Solve the orientation implications; return a full directed certificate.

    Each possible direction is a literal. It implies both required directions
    through every intervening vertex. A missing edge forces its negation.
    Transitive closures plus exhaustive branching over unassigned edges give
    satisfiability; every returned assignment is checked directly.
    """
    literals=[d for a,b in edges for d in [(a,b),(b,a)]]
    ix={d:i for i,d in enumerate(literals)}
    reach=[1<<i for i in range(len(literals))]
    pos={v:i for i,v in enumerate(order)}
    for i,(a,b) in enumerate(literals):
        for step in range(1,(pos[b]-pos[a])%6):
            w=order[(pos[a]+step)%6]
            for pair in [(a,w),(w,b)]:
                reach[i]|=1<<(ix[pair] if pair in ix else i^1)
    for k in range(len(literals)):
        for i in range(len(literals)):
            if reach[i]>>k&1:reach[i]|=reach[k]
    even=sum(1<<i for i in range(0,len(literals),2));odd=even<<1
    def solve(mask):
        if mask&(((mask&even)<<1)|((mask&odd)>>1)):return None
        for i in range(0,len(literals),2):
            if not mask&(3<<i):
                a=solve(mask|reach[i])
                return a if a is not None else solve(mask|reach[i+1])
        return mask
    mask=solve(0)
    if mask is None:return None
    directed=[pair for i,pair in enumerate(literals) if mask>>i&1]
    assert orientation_valid(order,set(directed))
    return directed


def orientation_valid(order,directed):
    pos={v:i for i,v in enumerate(order)}
    return all(all((a,order[(pos[a]+step)%6]) in directed and
                   (order[(pos[a]+step)%6],b) in directed
                   for step in range(1,(pos[b]-pos[a])%6))
               for a,b in directed)


def filtration_certificate():
    exact,D,levels,bounds=exact_metric()
    orders=[(0,)+p for p in itertools.permutations(range(1,6))]
    common=set(range(len(orders)));reports=[];direct_checks=0
    for level in levels:
        lo=bounds(level)
        def present(v):
            if v==level:return True
            # Interval comparisons are exact integers on the common denominator.
            vb=bounds(v)
            if int(vb['upper_numerator'])<int(lo['lower_numerator']):return True
            assert int(lo['upper_numerator'])<int(vb['lower_numerator'])
            return False
        edges=[(a,b) for a,b in itertools.combinations(range(6),2) if present(exact[a][b])]
        valid=[];certificates=[]
        for i,order in enumerate(orders):
            directed=cyclic_orientation(order,edges)
            if directed is not None:
                valid.append(i);certificates.append(directed)
            # Independently enumerate all orientations at the two proof graphs.
            if level in [(6,0,0,0),(0,0,0,3)]:
                brute=sum(orientation_valid(order,
                          {(b,a) if bit else (a,b) for (a,b),bit in zip(edges,bits)})
                          for bits in itertools.product([0,1],repeat=len(edges)))
                assert bool(brute)==(directed is not None)
                direct_checks+=1
        common &= set(valid)
        reports.append(dict(exact_coefficients=list(level),
                            threshold_float=float(sum(level[i]*[1,np.sqrt(2),np.sqrt(5),np.sqrt(10)][i] for i in range(4))/5),
                            threshold_bounds=lo,edges=edges,
                            admissible_orders=[orders[i] for i in valid],
                            orientation_certificates=certificates,
                            common_orders_through_threshold=[orders[i] for i in sorted(common)]))
    assert not common and all(r['admissible_orders'] for r in reports)
    early=next(r for r in reports if r['exact_coefficients']==[6,0,0,0])
    assert early['admissible_orders']==[(0,1,2,3,4,5),(0,5,4,3,2,1)]
    return dict(status='passed',vertex_names=['a','u','v','b','q','z'],
                cyclic_orders_checked=120,threshold_graphs=len(reports),
                exact_metric_basis=['1/5','sqrt(2)/5','sqrt(5)/5','sqrt(10)/5'],
                metric_exact_coefficients=exact,
                comparison_arithmetic='Integer outward square-root bounds at 80 decimal places; equal vectors are tied.',
                direct_enumeration_order_threshold_checks=direct_checks,
                common_orders=[],each_threshold_individually_cyclic=True,
                thresholds=reports)


def build_report():
    xy=np.array([(1,0),(24/25,7/25),(0,1),(-4/5,3/5),(-1,0),(0,-1)])
    chord=np.linalg.norm(xy[:,None]-xy[None,:],axis=2)
    edges=sorted({tuple(sorted((i,(i+1)%6))) for i in range(6)}|{(0,2),(0,3)})
    weights=np.full((6,6),np.inf);np.fill_diagonal(weights,0)
    for a,b in edges:weights[a,b]=weights[b,a]=chord[a,b]
    D=shortest_path(weights,directed=False)
    assert supported(D,chord,edges)
    threshold=1.9
    H=[(a,b) for a,b in itertools.combinations(range(6),2) if D[a,b]<=threshold]
    cyclic_count=0
    for bits in itertools.product([0,1],repeat=len(H)):
        directed={(b,a) if bit else (a,b) for (a,b),bit in zip(H,bits)}
        valid=True
        for a,b in directed:
            for step in range(1,(b-a)%6):
                v=(a+step)%6
                if (a,v) not in directed or (v,b) not in directed:valid=False;break
            if not valid:break
        cyclic_count+=int(valid)
    assert cyclic_count==0
    triangles=[(a,b,c) for a,b,c in itertools.combinations(range(6),3)
               if all(tuple(sorted(e)) in H for e in [(a,b),(a,c),(b,c)])]
    boundary=np.zeros((len(H),len(triangles)),int)
    for j,(a,b,c) in enumerate(triangles):
        for edge,sign in [((b,c),1),((a,c),-1),((a,b),1)]:boundary[H.index(edge),j]=sign
    ranks={str(p):rank_mod(boundary,p) for p in [2,3,5]}
    assert len(H)==10 and len(triangles)==4 and all(r==4 for r in ranks.values())
    assert D[0,3]<threshold<2<min(D[1,3],D[0,4])
    report=dict(status='passed',vertices=6,graph_edges=len(edges),arc_supported=True,
                threshold=threshold,threshold_edges=len(H),triangles=triangles,
                orientations_tested=2**len(H),cyclic_orientations=cyclic_count,
                boundary_ranks=ranks,betti_1=1,
                distances=dict(ab=float(D[0,3]),ub=float(D[1,3]),aq=float(D[0,4])),
                scope='No common cyclic vertex order throughout the threshold filtration, although each threshold is individually cyclic. All-field independence follows from a triangle-specific edge in each boundary.',
                filtration=filtration_certificate())
    return report

def main():
    report=build_report()
    (PAPER/'results/separation_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Six-point separation:',report['orientations_tested'],'orientations, none cyclic; beta1=1.')

if __name__=='__main__':main()
