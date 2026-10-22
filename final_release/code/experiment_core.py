"""Seed-defined sampling, graph evaluation and common output helpers."""
from pathlib import Path
import json,math,os
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import graph_methods as m
ROOT=Path(__file__).resolve().parents[1]
CONFIG=json.loads((ROOT/'code/config.json').read_text())
SCHEMAS=json.loads((ROOT/'code/result_schemas.json').read_text())
RESULTS=ROOT/'results'
assert set(CONFIG['noise']['models'])=={'radial','isotropic'}, 'Unknown noise model'
assert CONFIG['n']==512, 'The full-diagram protocol uses n=512'

def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(value,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
    temporary.replace(path)
def table(name,rows):
    path=RESULTS/name;path.parent.mkdir(parents=True,exist_ok=True)
    df=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
    if name in SCHEMAS:df=df.reindex(columns=SCHEMAS[name])
    df.to_csv(path,index=False)
    return df
def grid_values():
    a=CONFIG['clean']['a_grid']
    values=np.round(np.geomspace(a['minimum'],a['maximum'],a['points']),a['round_decimals'])
    return [(f'a{x:g}',int(round(float(x)*CONFIG['n']*math.log(CONFIG['n'])))) for x in values]
def ellipse_points(s):
    theta=np.linspace(0,2*np.pi,CONFIG['ellipse']['arc_grid_points'])
    speed=np.sqrt(4*np.sin(theta)**2+np.cos(theta)**2)
    arc=np.r_[0,np.cumsum((speed[1:]+speed[:-1])*.5*np.diff(theta))]
    perimeter=float(arc[-1]);arc/=perimeter
    t=np.interp(np.asarray(s)%1,arc,theta)
    return np.c_[2*np.cos(t),np.sin(t)]/perimeter
def all_orders(points,seed):
    chord,iu,ju,length,orders,keys=m.build_orders(points,include_cknn=True)
    n=len(points);mask=chord.copy();np.fill_diagonal(mask,np.inf)
    nbr=np.argsort(mask,axis=1,kind='mergesort');k0=max(8,int(math.ceil(2*math.log(n))))
    rho=mask[np.arange(n),nbr[:,k0-1]]
    clean=CONFIG['clean'];permutation=np.random.default_rng(seed+clean['scale_permutation_offset']+clean['scale_permutation_n_factor']*n).permutation(n)
    permuted=rho[permutation]
    orders['R3_permuted']=np.lexsort((ju,iu,length/np.sqrt(permuted[iu]*permuted[ju])))
    return chord,iu,ju,length,orders,keys
def evaluate(chord,iu,ju,length,order,count,true,reference,conn):
    from run_studies import measurements
    n=len(chord);chosen=order[:count]
    if count<conn:
        graph=csr_matrix((np.ones(2*len(chosen)),(np.r_[iu[chosen],ju[chosen]],np.r_[ju[chosen],iu[chosen]])),shape=(n,n))
        components=int(connected_components(graph,directed=False,return_labels=False))
        return dict(connected=False,n_components=components,e_inf=np.nan,e_2=np.nan,p1=np.nan,est_h1_bars=np.nan,est_longest_lifetime=np.nan,stability_ok=np.nan),np.empty((0,2))
    D=chord.copy() if count==len(order) else m.graph_metric(order,iu,ju,length,count,n)
    score,bars=measurements(D,true,reference,maximum_edge=float(length[chosen].max()))
    diameter=float(true.max())
    result=dict(connected=True,n_components=1,e_inf=score['e_inf'],e_2=score['e_2'],p1=score['e_PH'],
                est_h1_bars=score['bars'],est_longest_lifetime=score['longest']/diameter,
                stability_ok=bool(score['e_PH']<=score['e_inf']+5e-8/diameter))
    return result,bars
def rank_containment(position,s):
    n=len(s);bad=steps=antipodal=0;closest=1.
    for a in range(n):
        for direction in (1,-1):
            idx=(a+direction*np.arange(1,n))%n
            arc=((s[idx]-s[a])*direction)%1
            antipodal+=int(np.sum(np.abs(arc-.5)<1e-12))
            closest=min(closest,float(np.min(np.abs(arc-.5))))
            vals=position[a,idx[arc<.5-1e-12]]
            bad+=int(np.sum(vals[1:]<=vals[:-1]));steps+=max(len(vals)-1,0)
    return bad,steps,antipodal//2,closest
def interval_row(s,reference,diameter,method,count,close,bars):
    from run_studies import one_bar_distance
    g=float(np.diff(np.r_[s,s[0]+1]).max());B=float(np.sin(np.pi*g)/np.pi)
    error=one_bar_distance(bars,reference)/diameter
    baseline=float(np.max(reference[:,1]-reference[:,0])/2/diameter)
    return dict(method=method,n=len(s),M=count,M_close=close,g=g,B=B,P_star=baseline,bars=len(bars),
                birth=float(bars[0,0]) if len(bars)==1 else None,death=float(bars[0,1]) if len(bars)==1 else None,
                e_PH=error,clipped_error=min(error,baseline),birth_deviation=float(abs(bars[0,0]-B)) if len(bars)==1 else None,
                true_birth=float(reference[0,0]),true_death=float(reference[0,1]),diameter=diameter)
def transitions(rows,tolerance=1e-7):
    out=[]
    for a,b in zip(rows,rows[1:]):
        out.append(dict(method=a['method'],M_before=a['M'],M_after=b['M'],g=a['g'],P_star=a['P_star'],
                        before=a['e_PH'],after=b['e_PH'],delta=b['e_PH']-a['e_PH'],
                        clipped_delta=b['clipped_error']-a['clipped_error'],
                        death_delta=b['death']-a['death'] if a['bars']==b['bars']==1 else None,
                        raw_drop=b['e_PH']<a['e_PH']-tolerance,clipped_drop=b['clipped_error']<a['clipped_error']-tolerance,
                        crosses_below_baseline=a['e_PH']>=a['P_star']+tolerance and b['e_PH']<b['P_star']-tolerance))
    return out
