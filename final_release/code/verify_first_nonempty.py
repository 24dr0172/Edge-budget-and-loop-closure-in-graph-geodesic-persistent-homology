"""Check the stopping-rule audit, baseline choices and bounded independent replays."""
from pathlib import Path
import argparse,csv,json,hashlib,sys
sys.dont_write_bytecode=True
import numpy as np
from graph_methods import build_orders,first_connected,closure,graph_metric
from finite_diagram import matching_distance,reduction_h1
from first_nonempty_check import summarize,query,bridgeless,update_metric,METHODS
ROOT=Path(__file__).resolve().parents[1]

def read(path):
    raw=list(csv.DictReader(path.open()))
    bools={'at_closure','early','short_hops_first','first_ring_nonzero','first_all_bars_die_before_ring_chain','bridge_free_ring_present','double_ring_present'}
    integers={'seed','n','M_conn','M_close','M_first','first_bars','prefix_queries','metric_rebuilds','quantized_disagreements','k_bridge_free','k_close','M_bridge_free','M_bridge_free_layer','bridge_free_bars','M_double','double_bars'}
    for row in raw:
        for k,v in row.items():
            if not v:row[k]=None
            elif k in bools or v in ('True','False'):row[k]=v=='True'
            elif k in integers:row[k]=int(v)
            elif k!='method':row[k]=float(v)
    return raw

def table_rows(summary):
    out=[]
    for alpha in (0.,.4,.6,.8,.9):
        groups=[next(r for r in summary['groups']if r['alpha']==alpha and r['method']==m)for m in METHODS]
        cells=[f'${alpha:g}$']+[str(r['at_closure'])for r in groups]+[str(r['bridge_free_ring_present'])for r in groups]+[str(groups[0]['double_ring_present'])]
        out.append(' & '.join(cells)+r'\\')
    return out

def verify(root=ROOT,require=None,replay=False):
    root=Path(root);checks=[]
    def must(name,ok,detail=None):
        if require is not None:require('first_nonempty_'+name,ok,detail)
        if not ok:raise AssertionError((name,detail))
        checks.append(name)
    folder=root/'results/first_nonempty';cfg=json.loads((folder/'config.json').read_text())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    must('script_hash',cfg['script_sha256']==sha(root/'code/first_nonempty_check.py'))
    must('source_hash',cfg['source_sha256']==sha(root/cfg['source']))
    for name,digest in cfg['helper_sha256'].items():must('helper_'+name,digest==sha(root/'code'/name))
    rows=read(folder/'trials.csv');old=read(root/'results/peanut_neck/trials.csv')
    prior={(r['alpha'],r['seed'],r['method']):r for r in old}
    must('inventory',len(rows)==len(prior)==300 and len({(r['alpha'],r['seed'],r['method'])for r in rows})==300)
    current=summarize(rows);saved=json.loads((folder/'summary.json').read_text());must('summary',saved==current)
    must('no_rounding_disagreement',saved['quantized_disagreements']==0 and saved['maximum_incremental_drift']<1e-12)
    must('paired_seeds',saved['paired_seeds']==20 and saved['geometric_clouds']==100)
    replay_cases=set()
    for a in (0.,.4,.6,.8,.9):
        for method in METHODS:
            group=[r for r in rows if r['alpha']==a and r['method']==method]
            must(f'group_{a}_{method}',len(group)==20 and len({r['seed']for r in group})==20)
            replay_cases.add((a,group[0]['seed'],method))
    for r in rows:
        key=(r['alpha'],r['seed'],r['method']);before=prior[key]
        must('landmarks_'+str(key),(r['M_conn'],r['M_close'])==(before['M_conn'],before['M_close']))
        must('sequential_range_'+str(key),r['M_conn']<=r['M_first']<=r['M_close'] and
             r['prefix_queries']==r['M_first']-r['M_conn']+1 and r['early']==(r['M_first']<r['M_close']) and r['at_closure']==(r['M_first']==r['M_close']))
        if r['short_hops_first']:must('insertion_prediction_'+str(key),r['at_closure'])
        with np.load(root/'results/peanut_neck/checkpoints'/f"a{r['alpha']:g}_seed{r['seed']}.npz")as f:
            xy,s,reference=f['xy'],f['s'],f['reference_bars']
        sep=np.abs(s[:,None]-s);diameter=float(np.minimum(sep,1-sep).max())
        file=folder/'checkpoints'/f"a{r['alpha']:g}_seed{r['seed']}_{r['method']}.npz"
        with np.load(file,allow_pickle=False)as f:
            bars=f['first_bars'];bridge=f['bridge_free_bars']
            must('positive_'+str(key),len(bars)==r['first_bars']>0 and np.isfinite(bars).all()and np.all(bars[:,1]>bars[:,0]))
            must('error_'+str(key),abs(matching_distance(bars,reference)/diameter-r['first_error'])<1e-12)
            must('bridge_error_'+str(key),len(bridge)==r['bridge_free_bars'] and abs(matching_distance(bridge,reference)/diameter-r['bridge_free_error'])<1e-12)
            must('ring_absence_'+str(key),r['first_all_bars_die_before_ring_chain']==(float(bars[:,1].max())<=r['first_ring_chain_birth']) and
                 (not r['first_all_bars_die_before_ring_chain'] or not r['first_ring_nonzero']))
            if r['method']=='R2_global':
                doubled=f['double_bars'];must('doubling_error_'+str(key),len(doubled)==r['double_bars']and abs(matching_distance(doubled,reference)/diameter-r['double_error'])<1e-12)
        if r['early'] and r['first_all_bars_die_before_ring_chain']:
            must('spurious_empty_cost_'+str(key),abs(r['first_error']-before['P_star'])<1e-12 and r['first_last_death']<.01 and r['first_ring_chain_birth']>.97)
        if r['early']:replay_cases.add(key)
        if replay and key in replay_cases:
            C,iu,ju,w,orders,keys=build_orders(xy);order=orders[r['method']];n=len(xy)
            D=graph_metric(order,iu,ju,w,r['M_first'],n);bars,_=query(D)
            if r['early'] and r['first_ring_nonzero']:
                sep=np.abs(s[:,None]-s);intrinsic=np.minimum(sep,1-sep);mask=~np.eye(n,dtype=bool)
                lam=float((D[mask]/intrinsic[mask]).min());BG=float(D[np.arange(n),np.roll(np.arange(n),-1)].max())
                must('independent_winding_'+str(key),lam>.4428 and BG<.01386 and .01386/.4428<1/3)
            with np.load(file)as f:must('independent_first_'+str(key),np.array_equal(bars,f['first_bars']))
            if r['M_first']>r['M_conn']:
                Dp=graph_metric(order,iu,ju,w,r['M_first']-1,n);bp,_=query(Dp)
                must('independent_predecessor_'+str(key),len(bp)==0)
            mb=bridgeless(order,iu,ju,n,r['M_conn']);must('independent_bridge_'+str(key),mb==r['M_bridge_free'])
            if r['method'] in keys:
                kk=keys[r['method']][order];kb=int(kk[mb-1]);ml=int(np.searchsorted(kk,kb,side='right'))
                must('independent_layer_'+str(key),kb==r['k_bridge_free']and ml==r['M_bridge_free_layer'])
            if r['method']=='R2_global':
                md=int(np.searchsorted(w[order],2*w[order[r['M_conn']-1]],side='right'))
                must('independent_doubling_'+str(key),md==r['M_double'])
    bisection_queries=[]
    if replay:
        from first_nonempty_rule import first_nonempty_circle
        for method in METHODS:
            record=next(r for r in rows if r['alpha']==0. and r['method']==method)
            with np.load(root/'results/peanut_neck/checkpoints'/f"a0_seed{record['seed']}.npz")as f:points=f['xy']
            for layers in ([False] if method=='R2_global' else [False,True]):
                got=first_nonempty_circle(points,method,whole_layers=layers)
                C,iu,ju,w,orders,keys=build_orders(points)
                if layers:
                    expected_k=int(keys[method][orders[method][record['M_close']-1]])
                    expected_M=int(np.searchsorted(keys[method][orders[method]],expected_k,side='right'))
                    must('bisection_layer_'+method,got['k']==expected_k and got['M']==expected_M)
                else:must('bisection_exact_'+method,got['M']==record['M_close'])
                must('bisection_query_count_'+str((method,layers)),got['persistence_queries']<=2+int(np.ceil(np.log2(len(w)))))
                bisection_queries.append(dict(method=method,whole_layers=layers,queries=got['persistence_queries']))
    # Independent finite F2 reductions on small connected Euclidean-weighted graphs.
    from scipy.spatial.distance import cdist
    rng=np.random.default_rng(20361025)
    for rep in range(24):
        n=7+rep%8;xy=rng.normal(size=(n,2));C=cdist(xy,xy);iu,ju=np.triu_indices(n,1);w=C[iu,ju]
        edges={(i,i+1)for i in range(n-1)}|{(0,n-1)}
        edges|={tuple(e)for e in zip(iu,ju)if rng.random()<.2}
        chosen=[k for k,e in enumerate(zip(iu,ju))if e in edges];rest=[k for k in range(len(w))if k not in chosen];order=np.array(chosen+rest)
        D=graph_metric(order,iu,ju,w,len(chosen),n);cap=float(w[chosen].max());small,_=query(D,cap);full,_=query(D);independent,_=reduction_h1(D)
        must('cutoff_independent_'+str(rep),bool(len(small))==bool(len(full))==bool(len(independent)))
        if rest:
            e=rest[0];old_D=D.copy();update_metric(D,int(iu[e]),int(ju[e]),float(w[e]));fresh=graph_metric(order,iu,ju,w,len(chosen)+1,n)
            must('metric_update_'+str(rep),np.max(np.abs(D-fresh))<1e-12)
    return dict(status='passed',checks=len(checks),method_records=300,groups=15,prefix_queries=saved['prefix_queries'],
        small_independent_graphs=24,bisection_queries=bisection_queries,replayed_first_prefixes=len(replay_cases)if replay else 0,
        replay_scope='One case per group and every premature stop; Dijkstra plus unrestricted Ripser.'if replay else 'Saved-data and independent small-graph checks; n=512 replay not requested.')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--replay',action='store_true');a=p.parse_args()
    answer=verify(replay=a.replay);(ROOT/'results/first_nonempty/verification.json').write_text(json.dumps(answer,indent=2)+'\n');print(json.dumps(answer,indent=2))
