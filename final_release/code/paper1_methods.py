"""Sampling, edge orders, shortest-path metrics and order certificates for the paper."""
from __future__ import annotations
from typing import List,Tuple
import math
import numpy as np
from scipy.optimize import brentq
from scipy.spatial.distance import cdist
CDF_INVERSION_TOL=1e-12
PHASE_SEED_OFFSET=50000017
NOISE_SEED_OFFSET=70000003
METHODS=('R2_global','OR_kNN','mutual_kNN')
PROFILES=('uniform','strong_cosine','localized_trough')
def circular_arc_intervals(tau: float, w: float = 0.25) -> List[Tuple[float, float]]:
    left = tau - w / 2.0
    right = tau + w / 2.0
    if left < 0.0:
        return [(left + 1.0, 1.0), (0.0, right)]
    if right > 1.0:
        return [(left, 1.0), (0.0, right - 1.0)]
    return [(left, right)]

def arc_intersection_length(s: float, intervals: List[Tuple[float, float]]) -> float:
    length = 0.0
    for a, b in intervals:
        overlap_a = max(0.0, a)
        overlap_b = min(s, b)
        if overlap_b > overlap_a:
            length += (overlap_b - overlap_a)
    return length

def cdf_function(s: float, profile_name: str, phi: float) -> float:
    if s <= 0.0:
        return 0.0
    if s >= 1.0:
        return 1.0
    if profile_name == "uniform":
        return float(s)
    if profile_name == "strong_cosine":
        return float(s + (0.8 / (2.0 * np.pi)) * (
            np.sin(2.0 * np.pi * (s - phi)) - np.sin(-2.0 * np.pi * phi)
        ))
    if profile_name == "localized_trough":
        tau = (phi + 0.5) % 1.0
        intervals = circular_arc_intervals(tau, 0.25)
        inter_len = arc_intersection_length(s, intervals)
        return float((9.0 / 7.0) * s - (8.0 / 7.0) * inter_len)
    raise ValueError(f"Unknown density profile: {profile_name}")

def invert_cdf_vector(u_arr: np.ndarray, profile_name: str, phi: float) -> Tuple[np.ndarray, float]:
    if profile_name == "uniform":
        return u_arr.copy(), 0.0

    s_arr = np.empty_like(u_arr)
    max_residual = 0.0
    for i, u in enumerate(u_arr):
        if u <= 0.0:
            s = 0.0
        elif u >= 1.0:
            s = 1.0
        else:
            s = brentq(lambda val: cdf_function(val, profile_name, phi) - u, 0.0, 1.0, xtol=1e-14, rtol=1e-14)
        res = abs(cdf_function(s, profile_name, phi) - u)
        if res > max_residual:
            max_residual = res
        s_arr[i] = s

    if max_residual > CDF_INVERSION_TOL:
        raise RuntimeError(f"CDF inversion residual {max_residual:.2e} exceeded tolerance {CDF_INVERSION_TOL:.2e}")
    return s_arr, max_residual

def sample_cloud(seed,profile,n=512,sigma=0.):
    phi=float(np.random.default_rng(seed+PHASE_SEED_OFFSET).uniform(0.,1.))
    u=np.sort(np.random.default_rng(seed).uniform(0.,1.,size=n))
    s,_=invert_cdf_vector(u,profile,phi)
    normal=np.column_stack((np.cos(2*np.pi*s),np.sin(2*np.pi*s)))
    displacement=sigma*np.random.default_rng(seed+NOISE_SEED_OFFSET).uniform(-1.,1.,size=n)
    return s,normal/(2*np.pi)+displacement[:,None]*normal

def intrinsic(s):
    sep=np.abs(s[:,None]-s[None,:]);return np.minimum(sep,1.-sep)

def build_orders(observed,include_cknn=False):
    n=len(observed);chord=cdist(observed,observed);iu,ju=np.triu_indices(n,1);length=chord[iu,ju]
    mask=chord.copy();np.fill_diagonal(mask,np.inf)
    nbr=np.argsort(mask,axis=1,kind='mergesort');rank=np.empty((n,n),np.int64)
    rank[np.arange(n)[:,None],nbr]=np.arange(1,n+1)[None,:]
    r_or=np.minimum(rank[iu,ju],rank[ju,iu]);r_mu=np.maximum(rank[iu,ju],rank[ju,iu])
    orders={'R2_global':np.lexsort((ju,iu,length)),
            'OR_kNN':np.lexsort((ju,iu,length,r_or)),
            'mutual_kNN':np.lexsort((ju,iu,length,r_mu))}
    if include_cknn:
        k0=max(8,int(math.ceil(2*math.log(n))));rho=mask[np.arange(n),nbr[:,k0-1]]
        orders['R1_local']=np.lexsort((ju,iu,length/np.sqrt(rho[iu]*rho[ju])))
    return chord,iu,ju,length,orders,{'OR_kNN':r_or,'mutual_kNN':r_mu}

def first_connected(order,iu,ju,n):
    parent=np.arange(n);components=n
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for pos,e in enumerate(order,1):
        a,b=find(iu[e]),find(ju[e])
        if a!=b:
            parent[a]=b;components-=1
            if components==1:return pos
    raise RuntimeError('Complete pair pool did not connect')

def closure(order,iu,ju,n):
    position=np.zeros((n,n),dtype=np.int64);position[iu[order],ju[order]]=np.arange(1,len(order)+1);position+=position.T
    return int(position[np.arange(n),np.roll(np.arange(n),-1)].max())

def graph_metric(order,iu,ju,length,m,n):
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import shortest_path
    e=order[:m];w=length[e]
    graph=csr_matrix((np.r_[w,w],(np.r_[iu[e],ju[e]],np.r_[ju[e],iu[e]])),shape=(n,n))
    d=shortest_path(graph,directed=False,method='D');d=(d+d.T)/2;np.fill_diagonal(d,0.)
    return d

def diagram(d):
    from ripser import ripser
    return np.asarray(ripser(d,distance_matrix=True,maxdim=1)['dgms'][1],float)

def score(d,reference,reference_diagram):
    from persim import bottleneck
    if not np.isfinite(d).all():return {'connected':False}
    diam=float(reference.max());err=(d-reference)[np.triu_indices(len(d),1)];bars=diagram(d)
    life=(bars[:,1]-bars[:,0]).max() if len(bars) else 0.
    return {'connected':True,'e_inf':float(np.max(np.abs(err)))/diam,'e_2':float(np.sqrt(np.mean(err**2)))/diam,
            'p1':float(bottleneck(bars,reference_diagram))/diam,'est_h1_bars':len(bars),'est_longest_lifetime':float(life)/diam}

def certificate(s, chord, order, iu, ju):
    n=len(s)
    pos=np.zeros((n,n),dtype=np.int64)
    pos[iu[order],ju[order]]=np.arange(1,len(order)+1)
    pos+=pos.T
    cyc_pos=pos[np.arange(n),np.roll(np.arange(n),-1)]
    cover=int(cyc_pos.max())
    cyc_cost=chord[np.arange(n),np.roll(np.arange(n),-1)]
    cumulative=np.r_[0,np.cumsum(np.r_[cyc_cost,cyc_cost])]
    doubled_pos=np.r_[cyc_pos,cyc_pos]
    examined=0
    for index in order[:cover]:
        a,b=int(iu[index]),int(ju[index])
        if s[b]-s[a]>.5:
            a,b=b,a
        end=b if b>a else b+n
        interior=np.arange(a+1,end)%n
        if not len(interior):
            continue
        examined+=1
        edge_pos=int(pos[a,b]); cost=float(chord[a,b])
        if np.any(doubled_pos[a:end]>=edge_pos):
            return cover,False,examined,{'failure':'consecutive_closure','edge_position':edge_pos,'a':a,'b':b}
        left=cumulative[a+1:end]-cumulative[a]
        right=cumulative[end]-cumulative[a+1:end]
        left_ok=(pos[a,interior]<edge_pos)&(chord[a,interior]<=cost)
        right_ok=(pos[interior,b]<edge_pos)&(chord[interior,b]<=cost)
        # Replacement paths use the exact earlier consecutive edges above.
        left_ok|=left<=cost
        right_ok|=right<=cost
        failures=np.where(~left_ok|~right_ok)[0]
        if len(failures):
            k=int(failures[0]);v=int(interior[k])
            return cover,False,examined,{'failure':'subpair_replacement','edge_position':edge_pos,
                                       'a':a,'b':b,'interior_vertex':v,
                                       'outer_chord':cost,'left_chain':float(left[k]),'right_chain':float(right[k])}
    return cover,True,examined,{}
