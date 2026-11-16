"""Local consecutive-edge rank calculations and tied-rank probabilities."""
from __future__ import annotations
import argparse, csv, json, math, time, sys
sys.dont_write_bytecode=True
from fractions import Fraction
from pathlib import Path
import numpy as np
from graph_methods import (sample_cloud, build_orders, first_connected,
                           invert_cdf_vector, PHASE_SEED_OFFSET, PROFILES)


def neighboring_keys(s):
    """Local count identity for sorted circular coordinates in [0, 1).

    The largest gap must be < 1/2 and distances to consecutive neighbors
    must have no ties. Exact floating-point ties at the search boundaries are
    rejected; numerically near-tied inputs need separate care.
    Searching the sorted three-copy lift uses O(n log n) time and O(n) space.
    """
    s=np.asarray(s,dtype=float)
    if s.ndim!=1 or len(s)<3 or not np.isfinite(s).all():
        raise ValueError('Need a finite one-dimensional array of at least three points')
    if np.any(s<0) or np.any(s>=1) or np.any(np.diff(s)<=0):
        raise ValueError('Coordinates must be distinct and sorted in [0, 1)')
    n=len(s);lift=np.concatenate((s-1,s,s+1))
    a=np.arange(n)+n;b=a+1;g=lift[b]-lift[a]
    if g.max()>=.5: raise ValueError('Circle certificate requires g < 1/2')
    left_boundary=lift[a]-g;right_boundary=lift[b]+g
    left_lo=np.searchsorted(lift,left_boundary,side='left')
    left_hi=np.searchsorted(lift,left_boundary,side='right')
    right_lo=np.searchsorted(lift,right_boundary,side='left')
    right_hi=np.searchsorted(lift,right_boundary,side='right')
    if np.any(left_lo!=left_hi) or np.any(right_lo!=right_hi):
        raise ValueError('Equal-distance ties at consecutive neighbors are not supported')
    left=a-left_hi;right=right_lo-b-1
    return {'OR_kNN':1+np.minimum(left,right),
            'mutual_kNN':1+np.maximum(left,right)},left,right


def fast_inverse(u,profile,phi):
    """Same CDF as the released sampler, vectorized with residual <= 2e-13."""
    if profile=='uniform':return u.copy()
    if profile=='localized_trough':
        from graph_methods import circular_arc_intervals,cdf_function
        intervals=circular_arc_intervals((phi+.5)%1.,.25)
        cuts=np.array(sorted({0.,1.,*[z for pair in intervals for z in pair]}))
        cu=np.array([cdf_function(float(v),profile,phi) for v in cuts])
        return np.interp(u,cu,cuts)
    a=.8/(2*np.pi);offset=np.sin(-2*np.pi*phi)
    lo=np.zeros_like(u);hi=np.ones_like(u);s=u.copy()
    for _ in range(48):
        residual=s+a*(np.sin(2*np.pi*(s-phi))-offset)-u
        if np.max(np.abs(residual))<=2e-13: return s
        lo=np.where(residual<0,s,lo);hi=np.where(residual>0,s,hi)
        step=s-residual/(1+.8*np.cos(2*np.pi*(s-phi)))
        proposal=np.where((step>lo)&(step<hi),step,(lo+hi)/2)
        s=np.where(np.abs(residual)<=2e-13,s,proposal)
    raise RuntimeError('Vectorized cosine inverse did not reach tolerance')


def tails(kmax=120):
    """Exact rational one-gap Poisson tails, converted only after cancellation."""
    tor=[Fraction(1)]
    for k in range(1,kmax+1):
        correction=sum((Fraction(math.comb(k+j-1,j),3**(k+j))
                        for j in range(k)),Fraction())
        tor.append(Fraction(2,2**k)-2*correction)
    tmu=[Fraction(1)]+[Fraction(2,2**k)-tor[k] for k in range(1,kmax+1)]
    return {'OR_kNN':np.array([float(t) for t in tor]),
            'mutual_kNN':np.array([float(t) for t in tmu])}


def iid_extra_cdf(n,tail,j):
    # For Z=min(L,R) or max(L,R), tail[m]=P(Z >= m), p[m]=P(Z=m).
    # The gap exceeds j iff exactly one mark is m and all others are < m-j.
    p=tail[:-1]-tail[1:];m=np.arange(j+1,len(p))
    logs=np.log1p(-tail[m-j]);terms=n*p[m]*np.exp((n-1)*logs)
    return float(1-terms.sum())


def wilson(successes,total):
    z=1.959963984540054;p=successes/total;den=1+z*z/total
    center=(p+z*z/(2*total))/den
    half=z*math.sqrt(p*(1-p)/total+z*z/(4*total*total))/den
    return center-half,center+half


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def audit_local(seed_start,tests):
    checks=0;residual=0.;rank_checks=0;connectivity_checks=0
    for index in range(tests):
        seed=seed_start+index;profile=PROFILES[index%3];n=16+(index%49)
        s,xy=sample_cloud(seed,profile,n=n)
        if np.max(np.diff(np.r_[s,s[0]+1]))>=.5: continue
        phi=float(np.random.default_rng(seed+PHASE_SEED_OFFSET).uniform())
        u=np.sort(np.random.default_rng(seed).uniform(size=n))
        sv=fast_inverse(u,profile,phi)
        err=float(np.max(np.abs(s-sv)));residual=max(residual,err)
        if err>2e-12:raise AssertionError('Vectorized sample differs from released sample')
        keys,left,right=neighboring_keys(s)
        vector_keys,_,_=neighboring_keys(sv)
        if any(not np.array_equal(keys[m],vector_keys[m]) for m in keys):
            raise AssertionError('Vectorized sampler changes a neighboring rank')
        chord,iu,ju,length,orders,ranks=build_orders(xy)
        cyc=(ju==iu+1)|((iu==0)&(ju==n-1))
        edge_indices=np.where(cyc)[0]
        for method in ('OR_kNN','mutual_kNN'):
            full=np.empty(n,int)
            for e in edge_indices:
                i,j=iu[e],ju[e];idx=i if j==i+1 else n-1
                full[idx]=ranks[method][e]
            if not np.array_equal(full,keys[method]):raise AssertionError('Local identity failed')
            rank_checks+=n
            top=np.partition(keys[method],-2)[-2:]
            conn_edge=orders[method][first_connected(orders[method],iu,ju,n)-1]
            if ranks[method][conn_edge]!=min(top):raise AssertionError('Full connectivity key failed')
            connectivity_checks+=1
        checks+=1
    return {'samples':checks,'neighbor_keys_checked':rank_checks,
            'full_connectivity_checks':connectivity_checks,
            'maximum_sampler_coordinate_difference':residual}


def audit_neighboring_input_guards():
    """Check that unsupported inputs fail instead of producing misleading ranks."""
    invalid=[[], [0.,.2], [0.,.2,.2,.5], [0.,.2,.4,.6,.8,1.],
             [0.,.25,.5,.75], [0.,.1,.2,.3], [0.,.2,float('nan')],
             [[0.,.2,.4]], [-.1,.2,.4,.7]]
    for sample in invalid:
        try: neighboring_keys(sample)
        except ValueError: continue
        raise AssertionError(f'Unsupported circular coordinates accepted: {sample!r}')
    return len(invalid)

