"""Check the varying-waist records, including cone and winding certificates."""
from pathlib import Path
import csv, json, hashlib
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import shortest_path
from ripser import ripser
from graph_methods import build_orders, first_connected, closure, intrinsic
from peanut_neck_study import ALPHAS, METHODS, BASE_SEED, N, PROBE, sample

def shortcut_transition(job):
    """Replay the prefix immediately before a recorded cross-neck insertion."""
    root, row = job
    root = Path(root)
    alpha, seed, method = row['alpha'], int(row['seed']), row['method']
    with np.load(root/'results/peanut_neck/checkpoints'/f'a{alpha:g}_seed{seed}.npz', allow_pickle=False) as saved:
        C, iu, ju, weights, orders, _ = build_orders(saved['xy'])
        order = orders[method]
        budget = int(row['M_shortcut'])
        selected = order[:budget-1]
        adjacency = np.full_like(C, np.inf)
        np.fill_diagonal(adjacency, 0.)
        adjacency[iu[selected], ju[selected]] = weights[selected]
        adjacency[ju[selected], iu[selected]] = weights[selected]
        distance = shortest_path(adjacency, directed=False, method='FW')
        before = ripser(distance, distance_matrix=True, maxdim=1)['dgms'][1]
        edge = order[budget-1]
        i, j = int(iu[edge]), int(ju[edge])
        separation = abs(float(saved['s'][i]-saved['s'][j]))
        pair_arc=np.abs(saved['s'][iu[order]]-saved['s'][ju[order]])
        first_long=int(np.flatnonzero(np.minimum(pair_arc,1-pair_arc)>=PROBE)[0])+1
        return dict(alpha=alpha, seed=seed, method=method, budget=budget,
                    first_long_budget=first_long,
                    edge_matches=(i==int(row['shortcut_i']) and j==int(row['shortcut_j'])),
                    arc=min(separation, 1-separation),
                    cross_neck=bool(saved['xy'][i,1]*saved['xy'][j,1]<0),
                    before=before.tolist(), after=saved[f'{method}_shortcut_bars'].tolist())

def verify_peanut_neck(root, require, replay=True):
    root = Path(root)
    folder = root/'results/peanut_neck'
    cfg = json.loads((folder/'config.json').read_text())
    summary = json.loads((folder/'summary.json').read_text())
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    require('peanut_neck_protocol', cfg['n']==N and cfg['samples_per_width']==20 and
            cfg['alphas']==list(ALPHAS) and cfg['probe']==PROBE and cfg['noise'] is False and
            cfg['source_sha256']==sha(root/'code/peanut_neck_study.py') and
            all(sha(root/'code'/k)==v for k,v in cfg['helper_sha256'].items()))
    with (folder/'trials.csv').open() as f:
        raw = list(csv.DictReader(f))
    bools = {'short_hops_first','euclidean_reference_discarded','closure_supported','shortcut_supported'}
    rows = [{k:(v=='True' if k in bools else v if k=='method' else float(v)) for k,v in row.items()} for row in raw]
    require('peanut_neck_inventory', len(rows)==300 and len(list((folder/'checkpoints').glob('*.npz')))==100 and
            len({(v['alpha'],v['seed'],v['method']) for v in rows})==300)
    groups = []
    by = {(v['alpha'],int(v['seed']),v['method']):v for v in rows}
    for alpha in ALPHAS:
        for seed in range(BASE_SEED,BASE_SEED+20):
            saved = np.load(folder/'checkpoints'/f'a{alpha:g}_seed{seed}.npz',allow_pickle=False)
            s, xy, length = sample(seed,alpha)
            require(f'peanut_neck_seed_{alpha}_{seed}', np.array_equal(s,saved['s']) and
                    np.array_equal(xy,saved['xy']))
            C = cdist(xy,xy)
            dstar = intrinsic(s)
            reference = saved['reference_bars']
            diameter = float(dstar.max())
            half_reference = float((reference[0,1]-reference[0,0])/2)
            require(f'peanut_neck_reference_{alpha}_{seed}', len(reference)==1 and
                    abs(reference[0,0]-np.diff(np.r_[s,s[0]+1]).max())<5e-8 and
                    reference[0,0]<PROBE<reference[0,1])
            for method in METHODS:
                v=by[alpha,seed,method]
                require(f'peanut_neck_geometry_{alpha}_{seed}_{method}',
                        v['n']==N and abs(v['waist_width']-2*(1-alpha)/length)<1e-14 and
                        abs(v['P_star']-half_reference/diameter)<1e-14 and
                        v['M_conn']<=v['M_close']<v['M_shortcut'])
                for event in ('connection','closure','shortcut','euclidean'):
                    bars=saved['euclidean_bars' if event=='euclidean' else f'{method}_{event}_bars']
                    lives=(bars[:,1]-bars[:,0])/2 if len(bars) else np.array([])
                    direct=[]
                    for i,bar in enumerate(bars):
                        direct.append(max(float(np.abs(bar-reference[0]).max()),
                                          max(np.delete(lives,i),default=0.)))
                    diagonal=max(half_reference,max(lives,default=0.))
                    error=min(diagonal,min(direct,default=float('inf')))/diameter
                    require(f'peanut_neck_barcode_{alpha}_{seed}_{method}_{event}',
                            np.isfinite(bars).all() and (bars[:,1]>bars[:,0]).all() and
                            v[event+'_bars']==len(bars) and
                            v[event+'_alive']==int(((bars[:,0]<=PROBE)&(bars[:,1]>PROBE)).sum()) and
                            abs(v[event+'_error']-error)<1e-12)
                    if event=='euclidean':
                        require(f'peanut_neck_matching_{alpha}_{seed}_{method}',
                                v['euclidean_reference_discarded']==
                                (diagonal+5e-8*diameter<min(direct,default=float('inf'))))
                center=int(v['euclidean_cone_vertex'])
                require(f'peanut_neck_cone_{alpha}_{seed}_{method}',
                        abs(np.sqrt(np.max(np.sum((xy-xy[center])**2,axis=1)))-v['euclidean_cone_radius'])<1e-13 and
                        (v['euclidean_cone_radius']>=PROBE or v['euclidean_alive']==0))
            saved.close()
    for alpha in ALPHAS:
        for method in METHODS:
            group=[v for v in rows if v['alpha']==alpha and v['method']==method]
            rec=next(v for v in summary['groups'] if v['alpha']==alpha and v['method']==method)
            require(f'peanut_neck_group_scope_{alpha}_{method}',len(group)==rec['samples']==20)
            for event in ('connection','closure','shortcut','euclidean'):
                require(f'peanut_neck_group_{alpha}_{method}_{event}',
                        rec[event+'_empty']==sum(v[event+'_bars']==0 for v in group) and
                        rec[event+'_single']==sum(v[event+'_bars']==1 for v in group) and
                        rec[event+'_multiple']==sum(v[event+'_bars']>1 for v in group) and
                        rec[event+'_alive_at_probe']==sum(v[event+'_alive']>0 for v in group) and
                        abs(rec[event+'_error_min']-min(v[event+'_error'] for v in group))<1e-12 and
                        abs(rec[event+'_error_max']-max(v[event+'_error'] for v in group))<1e-12)
            for key,field in [('short_hops_first','short_hops_first'),('closure_supported','closure_supported'),
                              ('shortcut_supported','shortcut_supported'),
                              ('euclidean_reference_discarded','euclidean_reference_discarded')]:
                require(f'peanut_neck_group_count_{alpha}_{method}_{key}',rec[key]==sum(v[field] for v in group))
            require(f'peanut_neck_group_certificates_{alpha}_{method}',
                    rec['closure_winding_certificate']==sum(v['closure_min_metric_ratio']>3*PROBE for v in group) and
                    rec['euclidean_cone_at_probe']==sum(v['euclidean_cone_radius']<PROBE for v in group))
    replays=0
    if replay:
        for alpha,method in [(0.,'R2_global'),(.4,'OR_kNN'),(.9,'mutual_kNN')]:
            saved=np.load(folder/'checkpoints'/f'a{alpha:g}_seed{BASE_SEED}.npz',allow_pickle=False)
            C,iu,ju,weights,orders,keys=build_orders(saved['xy']);order=orders[method]
            v=by[alpha,BASE_SEED,method];n=len(C);dstar=intrinsic(saved['s'])
            require(f'peanut_neck_event_replay_{alpha}_{method}',
                    first_connected(order,iu,ju,n)==v['M_conn'] and closure(order,iu,ju,n)==v['M_close'] and
                    int(np.flatnonzero(dstar[iu[order],ju[order]]>=.25)[0])+1==v['M_first_long_edge'])
            for event,key in [('connection','M_conn'),('closure','M_close'),('shortcut','M_shortcut')]:
                chosen=order[:int(v[key])]
                adjacency=np.full((n,n),np.inf);np.fill_diagonal(adjacency,0.)
                adjacency[iu[chosen],ju[chosen]]=weights[chosen]
                adjacency[ju[chosen],iu[chosen]]=weights[chosen]
                # Floyd-Warshall on an explicit infinity matrix is independent of the CSR Dijkstra generator.
                D=shortest_path(adjacency,directed=False,method='FW')
                bars=ripser(D,distance_matrix=True,maxdim=1)['dgms'][1]
                other=saved[f'{method}_{event}_bars']
                key_sort=lambda a:a[np.lexsort((a[:,1],a[:,0]))]
                require(f'peanut_neck_diagram_replay_{alpha}_{method}_{event}',
                        bars.shape==other.shape and np.allclose(key_sort(bars),key_sort(other),atol=2e-7,rtol=0.))
                if event=='closure':
                    mask=~np.eye(n,dtype=bool)
                    ratio=float(np.min(D[mask]/dstar[mask]))
                    require(f'peanut_neck_winding_replay_{alpha}_{method}',
                            abs(ratio-v['closure_min_metric_ratio'])<1e-10 and
                            (ratio<=3*PROBE or v['closure_alive']>0))
                replays+=1
            saved.close()
    transitions=0
    if replay:
        jobs=[(str(root),v) for v in rows if v['alpha']>0]
        with ProcessPoolExecutor(max_workers=4) as pool:
            for v in pool.map(shortcut_transition,jobs):
                before, after = v['before'], v['after']
                require(f"peanut_neck_single_edge_{v['alpha']}_{v['seed']}_{v['method']}",
                        v['edge_matches'] and v['cross_neck'] and v['arc']>=PROBE and v['first_long_budget']==v['budget'] and
                        len(before)==1 and before[0][0]<=PROBE<before[0][1] and
                        len(after)==2 and all(b<d<=PROBE for b,d in after),v)
                transitions+=1
                if transitions%60==0:
                    print(f'One-edge peanut checks: {transitions}/240',flush=True)
    require('peanut_neck_summary_totals', summary['geometric_clouds']==100 and summary['method_records']==300 and
            summary['independent_seeds']==20)
    return dict(status='passed',geometric_clouds=100,independent_seeds=20,method_records=300,
                independent_diagram_replays=replays,single_edge_transition_replays=transitions,probe=PROBE)
