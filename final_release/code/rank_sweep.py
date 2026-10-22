"""Rank sweep and probability checks for kNN ring closure; no persistence diagrams are computed.
Exact circular neighboring ranks are compared with the full released chord-rank
kernel; global rank extremes are compared with an independent Poisson surrogate.
All generated files stay inside this release.
"""
from __future__ import annotations
import argparse, csv, json, math, time, sys
sys.dont_write_bytecode=True
from fractions import Fraction
from pathlib import Path
import numpy as np
from graph_methods import (sample_cloud, build_orders, first_connected,
                           invert_cdf_vector, PHASE_SEED_OFFSET, PROFILES)


def neighboring_keys(s):
    """Exact local count identity for sorted cyclic points with largest gap < 1/2.

    Counts use open intervals: no center, pair endpoint, or equal-distance tie.
    Searching the sorted three-copy lift uses O(n log n) time and O(n) space.
    """
    s=np.asarray(s,dtype=float)
    if len(s)<3 or np.any(np.diff(s)<=0):
        raise ValueError('Need at least three distinct points in sorted order')
    n=len(s);lift=np.concatenate((s-1,s,s+1))
    a=np.arange(n)+n;b=a+1;g=lift[b]-lift[a]
    if g.max()>=.5: raise ValueError('Circle certificate requires g < 1/2')
    left=a-np.searchsorted(lift,lift[a]-g,side='right')
    right=np.searchsorted(lift,lift[b]+g,side='left')-b-1
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


def verify_rank_sweep(root=None):
    """Recheck all stored summaries and model values; replay one seed per block."""
    root=Path(root) if root is not None else Path(__file__).resolve().parents[1]
    out=root/'results/extra_k/rank_sweep';checks=0
    with (out/'samples.csv').open() as f: rows=list(csv.DictReader(f))
    with (out/'summary.csv').open() as f: summaries=list(csv.DictReader(f))
    assert len(rows)==30780
    assert len({(r['profile'],r['n'],r['seed'],r['method']) for r in rows})==len(rows)
    for summary in summaries:
        group=[r for r in rows if all(r[v]==summary[v] for v in ('profile','n','method'))]
        extra=np.array([int(r['extra_k']) for r in group]);total=len(group)
        assert total==int(summary['samples']);checks+=1
        for j in range(6):
            count=int(np.count_nonzero(extra<=j));lo,hi=wilson(count,total)
            for suffix,value in (('',count/total),('_lo',lo),('_hi',hi)):
                assert abs(float(summary[f'cdf_{j}{suffix}'])-value)<1e-14;checks+=1
        for field,value in [('median_extra_k',np.median(extra)),('max_extra_k',extra.max()),
                            ('median_k_conn',np.median([int(r['k_conn']) for r in group])),
                            ('median_k_close',np.median([int(r['k_close']) for r in group]))]:
            assert abs(float(summary[field])-value)<1e-14;checks+=1
        n=int(summary['n']);row=group[0];seed=int(row['seed'])
        u=np.sort(np.random.default_rng(seed).uniform(size=n))
        phi=float(np.random.default_rng(seed+PHASE_SEED_OFFSET).uniform())
        s=fast_inverse(u,summary['profile'],phi);keys,_,_=neighboring_keys(s)
        top=np.sort(np.partition(keys[summary['method']],-2)[-2:])
        assert [int(row['k_conn']),int(row['k_close'])]==top.tolist();checks+=1
    assert max(float(r['largest_gap']) for r in rows)<1/3;checks+=1
    ts=tails()
    with (out/'iid_predictions.csv').open() as f:
        for row in csv.DictReader(f):
            for j in range(6):
                assert abs(float(row[f'cdf_{j}'])-iid_extra_cdf(int(row['n']),ts[row['method']],j))<1e-14
                checks+=1
    from scipy.integrate import quad
    from scipy.stats import poisson
    integral_errors=[]
    for k in range(1,13):
        value,_=quad(lambda x:np.exp(-x)*poisson.sf(k-1,x)**2,0,np.inf,epsabs=1e-15,epsrel=1e-12)
        error=abs(value-ts['OR_kNN'][k]);assert error<1e-13
        integral_errors.append(error);checks+=1
    return {'status':'passed','checks':checks,'method_records':len(rows),
            'profile_seed_samples':len(rows)//2,'summary_rows':len(summaries),
            'seed_replays':len(summaries),'poisson_integral_checks':12,
            'maximum_poisson_integral_error':max(integral_errors)}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify-only',action='store_true')
    args=parser.parse_args()
    if args.verify_only:
        print(json.dumps(verify_rank_sweep(),indent=2));return
    out=Path(__file__).resolve().parents[1]/'results/extra_k/rank_sweep'
    out.mkdir(parents=True,exist_ok=True)
    from experiment_core import CONFIG
    seed_start=CONFIG['ranks']['seed_start'];start=time.monotonic()
    validation=audit_local(seed_start,240)
    print('Local/full-rank verification:',validation,flush=True)
    settings=[tuple(pair) for pair in CONFIG['ranks']['settings']]
    rows=[];summary=[];tail=tails();predictions=[]
    if (out/'samples.csv').exists():
        with (out/'samples.csv').open() as f: rows=list(csv.DictReader(f))
        for row in rows:
            for field in ('n','seed','k_conn','k_close','extra_k'):row[field]=int(row[field])
            row['largest_gap']=float(row['largest_gap'])
        with (out/'summary.csv').open() as f: summary=list(csv.DictReader(f))
        for row in summary:
            for field in ('n','samples','max_extra_k'):row[field]=int(row[field])
            for field in row:
                if field not in ('profile','method','n','samples','max_extra_k'):row[field]=float(row[field])
    for n,reps in settings:
        for method,t in tail.items():
            predictions.append({'n':n,'method':method,**{f'cdf_{j}':iid_extra_cdf(n,t,j) for j in range(6)}})
    for profile_index,profile in enumerate(PROFILES):
        for setting_index,(n,reps) in enumerate(settings):
            if sum(r['profile']==profile and r['n']==n for r in rows)==2*reps:
                print('Resuming past completed block:',profile,n,flush=True);continue
            rows=[r for r in rows if not (r['profile']==profile and r['n']==n)]
            summary=[r for r in summary if not (r['profile']==profile and r['n']==n)]
            block=[]
            for rep in range(reps):
                seed=seed_start+10000+100000*setting_index+rep
                # Common uniforms and phase across profiles at each n/seed.
                u=np.sort(np.random.default_rng(seed).uniform(size=n))
                phi=float(np.random.default_rng(seed+PHASE_SEED_OFFSET).uniform())
                s=fast_inverse(u,profile,phi);keys,left,right=neighboring_keys(s)
                for method,key in keys.items():
                    top=np.sort(np.partition(key,-2)[-2:]);extra=int(top[1]-top[0])
                    row={'profile':profile,'n':n,'seed':seed,'method':method,
                         'k_conn':int(top[0]),'k_close':int(top[1]),'extra_k':extra,
                         'largest_gap':float(np.max(np.diff(np.r_[s,s[0]+1])))}
                    block.append(row);rows.append(row)
            for method in ('OR_kNN','mutual_kNN'):
                group=[r for r in block if r['method']==method]
                extra=np.array([r['extra_k'] for r in group])
                row={'profile':profile,'n':n,'method':method,'samples':reps,
                     'median_k_conn':float(np.median([r['k_conn'] for r in group])),
                     'median_k_close':float(np.median([r['k_close'] for r in group])),
                     'median_extra_k':float(np.median(extra)),'max_extra_k':int(extra.max())}
                for j in range(6):
                    count=int(np.count_nonzero(extra<=j));lo,hi=wilson(count,reps)
                    row.update({f'cdf_{j}':count/reps,f'cdf_{j}_lo':lo,f'cdf_{j}_hi':hi})
                summary.append(row)
                print(profile,n,method,'tie',round(row['cdf_0'],4),'+3',round(row['cdf_3'],4),'max',row['max_extra_k'],flush=True)
            write_csv(out/'samples.csv',rows);write_csv(out/'summary.csv',summary)
    write_csv(out/'iid_predictions.csv',predictions)
    tail_rows=[{'k':k,'or_tail':float(tail['OR_kNN'][k]),'mutual_tail':float(tail['mutual_kNN'][k]),
                'or_ratio':float(tail['OR_kNN'][k+1]/tail['OR_kNN'][k]),
                'mutual_ratio':float(tail['mutual_kNN'][k+1]/tail['mutual_kNN'][k])} for k in range(1,120)]
    write_csv(out/'poisson_tails.csv',tail_rows)
    config={'date':__import__('datetime').date.today().isoformat(),'seed_start':seed_start,'profiles':PROFILES,'settings':settings,
            'validation':validation,'samples':len(rows)//2,'method_records':len(rows),
            'protocol':'Independent seeds within n; same uniforms and phase across profiles. Exact circular local ranks; no new persistence computations. Local formula checked against full released chord ranks and actual graph connectivity.',
            'surrogate':'Independent one-gap Poisson marks. This is a comparison model, not a proof for dependent ranks on one circle.',
            'uncertainty':'Pointwise 95% Wilson binomial intervals, conditional on profile/n. Profiles share seeds.',
            'runtime_seconds_this_invocation':time.monotonic()-start}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    print('Completed',config['samples'],'samples in',round(config['runtime_seconds_this_invocation'],2),'seconds',flush=True)

if __name__=='__main__':main()
