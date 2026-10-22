"""Peanut/trough stress check at exact ring closure, with no lifetime cutoff.

150 independent clouds at n=128,256,512; the three methods share each cloud.
The trough density and random phase are those of the main rank sweep. The
curve and 400001-node arc-length interpolation are those of Appendix A.2.
This experiment checks closure only; it does not follow later prefixes.
Run: python3 code/check_peanut_trough.py [--workers 4]
"""
from pathlib import Path
import argparse,csv,hashlib,json,os
import numpy as np
from check_noncircular_locality import (arclength_table,curve,orderings,
    short_hops_first,components_at,graph_dist,edge_supported,h1_bars)
from rank_sweep import fast_inverse
from graph_methods import PHASE_SEED_OFFSET
ROOT=Path(__file__).resolve().parents[1]
BASE=20350401
SIZES=(128,256,512)
SAMPLES=150
METHODS=('radius','OR','mutual')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def setup():
    global TABLE
    TABLE=arclength_table('peanut')
def one(job):
    n,rep,out,r_loc=job;seed=BASE+10000*SIZES.index(n)+rep
    rng=np.random.default_rng(seed);u=np.sort(rng.random(n))
    phase=float(np.random.default_rng(seed+PHASE_SEED_OFFSET).uniform())
    s=fast_inverse(u,'localized_trough',phase)
    theta=np.interp(s*TABLE[1][-1],TABLE[1],TABLE[0])
    xy=curve('peanut',theta)/TABLE[1][-1]
    C=np.linalg.norm(xy[:,None]-xy[None,:],axis=2)
    iu,ju,orders,keys=orderings(C)
    gap=np.diff(np.r_[s,s[0]+1]);rs=C[np.arange(n),(np.arange(n)+1)%n];B=float(rs.max())
    ring_a=np.arange(n);ring_b=(ring_a+1)%n
    pair_id=np.zeros((n,n),int);pair_id[iu,ju]=np.arange(len(iu));pair_id[ju,iu]=np.arange(len(iu))
    intrinsic=np.minimum(np.abs(s[:,None]-s[None,:]),1-np.abs(s[:,None]-s[None,:]))
    ref=h1_bars(intrinsic)
    payload=dict(s=s,xy=xy,intrinsic_bars=np.array(ref,dtype=float).reshape(-1,2))
    rows=[];witness={}
    for name,order in orders.items():
        pos=np.zeros((n,n),int);pos[iu[order],ju[order]]=np.arange(len(order));pos[ju[order],iu[order]]=np.arange(len(order))
        mclose=int(pos[ring_a,ring_b].max())+1;mconn=components_at(order,iu,ju,n)+1
        ea,eb=iu[order[:mclose]],ju[order[:mclose]];D=graph_dist(C,ea,eb)
        assert np.all(D[ea,eb] <= C[ea,eb]*(1+1e-13)+1e-15), 'A graph edge was lost from shortest paths'
        bad=[(int(a),int(b)) for a,b in zip(ea,eb) if not edge_supported(D,C,rs,int(a),int(b),n)]
        bars=h1_bars(D);payload['bars_'+name]=np.asarray(bars,dtype=float).reshape(-1,2)
        shf=short_hops_first(C,pos,iu,ju,order,mclose-1,s)
        if name=='radius':local=B<r_loc;kclose=-1
        else:
            kclose=int(max(keys[name][pair_id[a,b]] for a,b in zip(ring_a,ring_b)))
            local=bool(np.sort(C,axis=1)[:,kclose].max()<r_loc)
        rows.append(dict(n=n,seed=seed,phase=phase,method=name,g=float(gap.max()),B=B,
            M_conn=mconn,M_close=mclose,k_close=kclose,locality_hypothesis=bool(local),
            short_hops_first=bool(shf),arc_supported=not bad,unsupported_edges=len(bad),
            bars=len(bars),bars_born_at_B=sum(abs(b-B)<=1e-6*max(1.,B) for b,d in bars),
            intrinsic_bars=len(ref),max_lifetime=max((d-b for b,d in bars),default=0.)))
        if bad:witness[name]=dict(unsupported_edges=bad,coordinates='see matching checkpoint xy',distances='recompute with graph_dist at M_close')
    stem=f'n{n}_seed{seed}';np.savez_compressed(out/'checkpoints'/f'{stem}.npz',**payload)
    if witness:(out/'witnesses'/f'{stem}.json').write_text(json.dumps(witness,indent=2)+'\n')
    else:(out/'witnesses'/f'{stem}.json').unlink(missing_ok=True)
    return rows
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);a=ap.parse_args()
    out=ROOT/'results/peanut_trough';(out/'checkpoints').mkdir(parents=True,exist_ok=True);(out/'witnesses').mkdir(exist_ok=True)
    loc=json.loads((ROOT/'results/noncircular_locality/locality_check_summary.json').read_text())['r_loc']['peanut']['r_loc']
    config=dict(curve='peanut r=1+0.4 cos(2 theta), length one',profile='localized_trough',
        density=dict(trough=1/7,elsewhere=9/7,trough_arc_fraction=.25),sizes=list(SIZES),samples_per_size=SAMPLES,
        base_seed=BASE,seed='base + 10000*size_index + replicate',phase_seed_offset=PHASE_SEED_OFFSET,
        r_loc=loc,threshold_source='results/noncircular_locality/locality_check_summary.json',arc_length_grid=400001,
        methods=list(METHODS),event='exact first prefix containing every consecutive pair',lifetime_cutoff=None,
        code_sha256=sha(Path(__file__)),distance_kernel_sha256=sha(ROOT/'code/check_noncircular_locality.py'),paired='Methods share clouds; different seeds and sizes are independent.')
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    jobs=[(n,j,out,loc) for n in SIZES for j in range(SAMPLES)];rows=[]
    from multiprocessing import Pool
    with Pool(max(1,min(a.workers,len(jobs))),initializer=setup) as pool:
        for count,group in enumerate(pool.imap_unordered(one,jobs),1):
            rows.extend(group)
            if count%25==0:print('peanut/trough',count,'/',len(jobs),flush=True)
    rows.sort(key=lambda v:(v['n'],v['seed'],v['method']))
    with (out/'closure.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=[]
    for n in SIZES:
        for method in METHODS:
            group=[v for v in rows if v['n']==n and v['method']==method]
            summary.append(dict(n=n,method=method,samples=len(group),locality_failures=sum(not v['locality_hypothesis'] for v in group),
                short_hops_failures=sum(not v['short_hops_first'] for v in group),arc_support_failures=sum(not v['arc_supported'] for v in group),
                multibar=sum(v['bars']>1 for v in group),empty=sum(v['bars']==0 for v in group),
                single_bar_born_at_B=sum(v['bars']==v['bars_born_at_B']==1 for v in group),
                positive_intrinsic=sum(v['intrinsic_bars']==1 for v in group),max_g=max(v['g'] for v in group)))
    result=dict(status='completed',clouds=len(jobs),method_records=len(rows),groups=summary)
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
