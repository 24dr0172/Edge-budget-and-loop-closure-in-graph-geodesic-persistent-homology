"""Audit the first nonempty connected prefix on the existing paired peanut clouds.

Every connected prefix is tested in increasing order, without assuming
monotonicity and without using the curve order to choose the stopping budget.
The known ring is used only afterwards to assess the selected graph.
Nonemptiness queries stop at the largest generating-edge cost: Lemma 1
implies that every later H1 class is an image of a class already present
there. Finite and infinite positive bars at that cutoff both count.
Full unrestricted diagrams verify each selected prefix and its predecessor.
"""
from pathlib import Path
import argparse,csv,hashlib,json,os,sys,time
sys.dont_write_bytecode=True
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_name]='1'
import numpy as np
from ripser import ripser
from graph_methods import build_orders,first_connected,closure,graph_metric
from finite_diagram import matching_distance
from bridge_free import traversal
from peanut_neck_study import ring_death

ROOT=Path(__file__).resolve().parents[1]
METHODS=('R2_global','OR_kNN','mutual_kNN')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def query(D,edge_ceiling=None,cocycles=False):
    cap=np.inf if edge_ceiling is None else float(np.nextafter(np.float32(edge_ceiling),np.float32(np.inf)))
    answer=ripser(D,distance_matrix=True,maxdim=1,thresh=cap,do_cocycles=cocycles)
    raw=np.asarray(answer['dgms'][1],float).reshape(-1,2)
    keep=raw[:,1]>raw[:,0]
    bars=raw[keep]
    co=[a for a,k in zip(answer['cocycles'][1],keep)if k] if cocycles else []
    return bars,co

def update_metric(D,a,b,cost):
    """Positive costs make a shortest new path use the added edge at most once."""
    da,db=D[:,a].copy(),D[:,b].copy()
    np.minimum(D,da[:,None]+cost+db[None,:],out=D)
    np.minimum(D,db[:,None]+cost+da[None,:],out=D)
    np.fill_diagonal(D,0.)

def select(order,iu,ju,w,n):
    conn=first_connected(order,iu,ju,n)
    D=graph_metric(order,iu,ju,w,conn,n)
    ceiling=float(w[order[:conn]].max())
    checks=0;rebuilds=1;drift=0.;quantized_disagreements=0
    for m in range(conn,len(order)+1):
        if m>conn:
            e=int(order[m-1]);update_metric(D,int(iu[e]),int(ju[e]),float(w[e]));ceiling=max(ceiling,float(w[e]))
        if (m-conn)%64==0 and m>conn:
            fresh=graph_metric(order,iu,ju,w,m,n);rebuilds+=1
            drift=max(drift,float(np.abs(fresh-D).max()))
            if not np.array_equal(fresh.astype(np.float32),D.astype(np.float32)):
                quantized_disagreements+=1
                # Do not make an unverified choice after numerical disagreement.
                raise AssertionError(('Incremental float32 disagreement',m,drift))
            D=fresh
        bars,_=query(D,ceiling);checks+=1
        if not len(bars):continue
        fresh=graph_metric(order,iu,ju,w,m,n);rebuilds+=1
        drift=max(drift,float(np.abs(fresh-D).max()))
        assert np.array_equal(fresh.astype(np.float32),D.astype(np.float32))
        full,co=query(fresh,cocycles=True)
        assert len(full)>0,'Truncated positive predicate disagrees with unrestricted persistence.'
        if m>conn:
            previous=graph_metric(order,iu,ju,w,m-1,n);rebuilds+=1
            before,_=query(previous)
            assert len(before)==0,'Selected prefix is not the first positive connected prefix.'
        return m,fresh,full,co,dict(prefix_queries=checks,metric_rebuilds=rebuilds,
            maximum_incremental_drift=drift,quantized_disagreements=quantized_disagreements)
    return None,None,np.empty((0,2)),[],dict(prefix_queries=checks,metric_rebuilds=rebuilds,
        maximum_incremental_drift=drift,quantized_disagreements=quantized_disagreements)

def bridgeless(order,iu,ju,n,conn):
    def passes(m):
        e=order[:m];c,b=traversal(n,list(zip(map(int,iu[e]),map(int,ju[e]))))
        return c==1 and b==0
    if passes(conn):return conn
    lo=conn;hi=min(len(order),max(conn+1,2*conn))
    while not passes(hi):lo=hi;hi=min(len(order),2*hi)
    while hi-lo>1:
        mid=(lo+hi)//2
        if passes(mid):hi=mid
        else:lo=mid
    return hi

def describe(D,bars,reference,diameter):
    n=len(D);ring_birth=float(D.astype(np.float32)[np.arange(n),np.roll(np.arange(n),-1)].max())
    last=float(bars[:,1].max()) if len(bars) else 0.
    return dict(bars=len(bars),error=matching_distance(bars,reference)/diameter,
        ring_chain_birth=ring_birth,last_death=last,
        all_bars_die_before_ring_chain=bool(last<=ring_birth),
        longest_lifetime=float((bars[:,1]-bars[:,0]).max()) if len(bars) else 0.)

def one(record):
    start=time.monotonic();alpha=float(record['alpha']);seed=int(record['seed']);method=record['method']
    file=ROOT/'results/peanut_neck/checkpoints'/f'a{alpha:g}_seed{seed}.npz'
    with np.load(file,allow_pickle=False)as f:xy=f['xy'];s=f['s'];reference=f['reference_bars']
    C,iu,ju,w,orders,keys=build_orders(xy);n=len(xy);order=orders[method]
    conn=first_connected(order,iu,ju,n);close=closure(order,iu,ju,n)
    assert (conn,close)==(int(record['M_conn']),int(record['M_close']))
    m,D,bars,co,trace=select(order,iu,ju,w,n)
    assert m is not None and conn<=m<=close
    intrinsic=np.minimum(np.abs(s[:,None]-s),1-np.abs(s[:,None]-s));diameter=float(intrinsic.max())
    details=describe(D,bars,reference,diameter)
    first_ring_death=ring_death(bars,co,details['ring_chain_birth'],n)
    row=dict(alpha=alpha,seed=seed,method=method,n=n,M_conn=conn,M_close=close,
        M_first=m,at_closure=m==close,early=m<close,short_hops_first=record['short_hops_first']=='True',
        first_ring_nonzero=first_ring_death>details['ring_chain_birth'],
        closure_error=float(record['closure_error']),**{'first_'+k:v for k,v in details.items()},**trace)
    assert not row['short_hops_first'] or m==close
    mb=bridgeless(order,iu,ju,n,conn)
    if method in keys:
        sorted_keys=keys[method][order];kbr=int(sorted_keys[mb-1]);mb_layer=int(np.searchsorted(sorted_keys,kbr,side='right'))
        row.update(k_bridge_free=kbr,k_close=int(sorted_keys[close-1]))
    else:mb_layer=mb;row.update(k_bridge_free=-1,k_close=-1)
    Db=graph_metric(order,iu,ju,w,mb_layer,n);bb,_=query(Db)
    row.update(M_bridge_free=mb,M_bridge_free_layer=mb_layer,bridge_free_ring_present=mb_layer>=close,
               bridge_free_bars=len(bb),bridge_free_error=matching_distance(bb,reference)/diameter)
    payload={'first_bars':bars,'bridge_free_bars':bb}
    if method=='R2_global':
        rconn=float(w[order[conn-1]]);md=int(np.searchsorted(w[order],2*rconn,side='right'))
        Dd=graph_metric(order,iu,ju,w,md,n);bd,_=query(Dd)
        row.update(radius_conn=rconn,M_double=md,double_ring_present=md>=close,
                   double_bars=len(bd),double_error=matching_distance(bd,reference)/diameter)
        payload['double_bars']=bd
    row['seconds']=time.monotonic()-start
    return row,payload

def summarize(rows):
    groups=[]
    for alpha in sorted({r['alpha']for r in rows}):
        for method in METHODS:
            g=[r for r in rows if r['alpha']==alpha and r['method']==method]
            groups.append(dict(alpha=alpha,method=method,samples=len(g),
                at_closure=sum(r['at_closure']for r in g),early=sum(r['early']for r in g),
                early_without_ring_class=sum(r['early']and not r['first_ring_nonzero']for r in g),
                early_bars_die_before_ring=sum(r['early']and r['first_all_bars_die_before_ring_chain']for r in g),
                bridge_free_ring_present=sum(r['bridge_free_ring_present']for r in g),
                bridge_free_nonempty=sum(r['bridge_free_bars']>0 for r in g),
                first_error_min=min(r['first_error']for r in g),first_error_max=max(r['first_error']for r in g),
                bridge_error_min=min(r['bridge_free_error']for r in g),bridge_error_max=max(r['bridge_free_error']for r in g),
                max_premature_edges=max(r['M_close']-r['M_first']for r in g),
                double_ring_present=sum(r.get('double_ring_present',False)for r in g)if method=='R2_global'else None,
                double_nonempty=sum(r.get('double_bars',0)>0 for r in g)if method=='R2_global'else None))
    return dict(status='completed',geometric_clouds=len(rows)//3,paired_seeds=len({r['seed']for r in rows}),
        method_records=len(rows),prefix_queries=sum(r['prefix_queries']for r in rows),
        maximum_incremental_drift=max(r['maximum_incremental_drift']for r in rows),
        quantized_disagreements=sum(r['quantized_disagreements']for r in rows),groups=groups)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--workers',type=int,default=4)
    a=parser.parse_args();out=ROOT/'results/first_nonempty';out.mkdir(parents=True,exist_ok=True);(out/'checkpoints').mkdir(exist_ok=True)
    source=ROOT/'results/peanut_neck/trials.csv';records=list(csv.DictReader(source.open()))
    config=dict(source='results/peanut_neck/trials.csv',source_sha256=sha(source),n=512,
        widths=[0.,.4,.6,.8,.9],seeds=20,methods=list(METHODS),new_clouds=False,
        search='Every connected exact prefix, in increasing order, until the first positive diagram. No monotonicity assumption.',
        numerical_convention='Ripser float32 filtration, F2, all positive returned bars, no lifetime cutoff.',
        empty_query='Persistence truncated just above the largest generating-edge cost; Lemma 1 rules out later births.',
        metric_updates='Exact one-new-edge shortest-path identity in float64; Dijkstra rebuild and float32 agreement every 64 insertions and at the selected prefix.',
        unrestricted_checks='Full diagrams at each selected prefix and its preceding connected prefix; full bridge-free and doubled-radius baselines.',
        comparison='Bridge-free: first connected bridgeless whole layer (radius layers are single insertions). Radius baseline: twice the smallest connecting radius.',
        script_sha256=sha(__file__),helper_sha256={n:sha(ROOT/'code'/n)for n in ['graph_methods.py','finite_diagram.py','bridge_free.py','peanut_neck_study.py']})
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    from multiprocessing import Pool
    rows=[]
    with Pool(max(1,a.workers))as pool:
        for index,(row,payload)in enumerate(pool.imap_unordered(one,records),1):
            rows.append(row);np.savez_compressed(out/'checkpoints'/f"a{row['alpha']:g}_seed{row['seed']}_{row['method']}.npz",**payload)
            if index%10==0 or index==len(records):print(f"First nonempty: {index}/{len(records)}, early so far {sum(r['early']for r in rows)}",flush=True)
    rows.sort(key=lambda r:(r['alpha'],r['seed'],METHODS.index(r['method'])))
    fields=list(rows[0])+[k for r in rows for k in r if k not in rows[0]];fields=list(dict.fromkeys(fields))
    with(out/'trials.csv').open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    summary=summarize(rows);(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
