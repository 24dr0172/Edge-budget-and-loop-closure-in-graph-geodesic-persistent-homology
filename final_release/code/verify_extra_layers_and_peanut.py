"""Check the extra-layer export and peanut/trough cohort against raw evidence.

The finite distance replays use independent Floyd-Warshall shortest paths,
including every structural-failure cloud and the near-coincident-edge case.
The radius proposition is analytical; its beta integral and chord bounds are
checked numerically here as corroboration, not as a proof.
"""
from pathlib import Path
import csv,hashlib,json,math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.integrate import quad
from scipy.special import gammaln
from ripser import ripser
ROOT=Path(__file__).resolve().parents[1]
def truth(v):
    assert v in ('True','False');return v=='True'
def verify_extra_layers_and_peanut(root,require):
    root=Path(root);base=root/'results/extra_k/rank_sweep'
    raw=list(csv.DictReader((base/'samples.csv').open()))
    export=list(csv.DictReader((root/'results/retention_extra_layers/retention.csv').open()))
    require('extra_layer_rows',len(export)==144)
    for row in export:
        n,j=int(row['n']),int(row['j']);group=[v for v in raw if int(v['n'])==n and v['profile']==row['profile'] and v['method']==row['method']]
        k=sum(int(v['k_close'])<=int(v['k_conn'])+j for v in group);p=k/len(group);z=1.959963984540054
        den=1+z*z/len(group);center=(p+z*z/(2*len(group)))/den;half=z*math.sqrt(p*(1-p)/len(group)+z*z/(4*len(group)**2))/den
        require('extra_layer_'+str((n,row['profile'],row['method'],j)),int(row['closed'])==k and int(row['samples'])==len(group) and
                all(abs(float(row[key])-value)<1e-14 for key,value in [('rate',p),('wilson_lo',center-half),('wilson_hi',center+half)]))
    cfg=json.loads((root/'results/retention_extra_layers/config.json').read_text())
    require('extra_layer_input_identity',cfg['source_sha256']==hashlib.sha256((base/'samples.csv').read_bytes()).hexdigest())
    for ext,folder in [('pdf','results/figure_exports'),('png','manuscript_figures')]:
        require('extra_layer_figure_'+ext,(root/f'{folder}/retention_extra_layers.{ext}').is_file())

    folder=root/'results/peanut_trough';cfg=json.loads((folder/'config.json').read_text())
    rows=list(csv.DictReader((folder/'closure.csv').open()));summary=json.loads((folder/'summary.json').read_text())
    require('peanut_trough_rows',len(rows)==1350 and len({(v['n'],v['seed'],v['method']) for v in rows})==1350)
    require('peanut_trough_code_identity',cfg['code_sha256']==hashlib.sha256((root/'code/check_peanut_trough.py').read_bytes()).hexdigest() and
            cfg['distance_kernel_sha256']==hashlib.sha256((root/'code/check_noncircular_locality.py').read_bytes()).hexdigest())
    min_chord=1.;replays=0;failure_clouds={(int(v['n']),int(v['seed'])) for v in rows if not truth(v['arc_supported'])}
    # Load each saved point cloud and compare all reported bars and gap statistics.
    from check_noncircular_locality import curve,arclength_table,orderings,sample
    from rank_sweep import fast_inverse
    from graph_methods import PHASE_SEED_OFFSET
    table=arclength_table('peanut')
    for n in (128,256,512):
        for rep in range(150):
            seed=20350401+10000*(128,256,512).index(n)+rep
            path=folder/'checkpoints'/f'n{n}_seed{seed}.npz';saved=np.load(path)
            u=np.sort(np.random.default_rng(seed).random(n));phase=float(np.random.default_rng(seed+PHASE_SEED_OFFSET).uniform())
            s=fast_inverse(u,'localized_trough',phase);xy=curve('peanut',np.interp(s*table[1][-1],table[1],table[0]))/table[1][-1]
            require(f'peanut_coordinates_{n}_{seed}',np.array_equal(saved['s'],s) and np.array_equal(saved['xy'],xy))
            C=np.linalg.norm(xy[:,None]-xy[None,:],axis=2);near=float(C[np.triu_indices(n,1)].min());min_chord=min(min_chord,near)
            B=float(C[np.arange(n),(np.arange(n)+1)%n].max());g=float(np.diff(np.r_[s,s[0]+1]).max())
            group=[v for v in rows if int(v['n'])==n and int(v['seed'])==seed]
            require(f'peanut_cloud_{n}_{seed}',len(group)==3 and {v['method'] for v in group}=={'radius','OR','mutual'} and
                    all(abs(float(v['B'])-B)<1e-14 and abs(float(v['g'])-g)<1e-14 for v in group))
            for v in group:
                bars=saved['bars_'+v['method']]
                require(f'peanut_bars_{n}_{seed}_{v["method"]}',int(v['bars'])==len(bars) and
                        int(v['bars_born_at_B'])==sum(abs(b-B)<=1e-6*max(1.,B) for b,d in bars) and
                        int(v['intrinsic_bars'])==len(saved['intrinsic_bars']))
            if rep==0 or (n,seed) in failure_clouds or near<1e-8:
                iu,ju,orders,keys=orderings(C)
                for v in group:
                    order=orders[v['method']];mc=int(v['M_close']);ea,eb=iu[order[:mc]],ju[order[:mc]]
                    W=csr_matrix((np.r_[C[ea,eb],C[ea,eb]],(np.r_[ea,eb],np.r_[eb,ea])),shape=(n,n))
                    D=shortest_path(W,directed=False,method='FW')
                    require(f'peanut_edges_preserved_{n}_{seed}_{v["method"]}',np.all(D[ea,eb]<=C[ea,eb]+1e-14))
                    diag=ripser(D,distance_matrix=True,maxdim=1)['dgms'][1];diag=diag[diag[:,1]>diag[:,0]]
                    require(f'peanut_independent_metric_diagram_{n}_{seed}_{v["method"]}',np.array_equal(diag,saved['bars_'+v['method']]))
                    # Compute the fan condition directly from the definition.
                    unsupported=[]
                    for a,b in zip(ea,eb):
                        forward=np.arange(a,b+1);backward=np.r_[np.arange(b,n),np.arange(a+1)]
                        supported=False
                        for arc in (forward,backward):
                            if max(C[arc[:-1],arc[1:]])<=C[a,b] and min(max(D[a,arc]),max(D[b,arc]))<=C[a,b]+1e-14:
                                supported=True;break
                        if not supported:unsupported.append((int(a),int(b)))
                    require(f'peanut_independent_support_{n}_{seed}_{v["method"]}',(not unsupported)==truth(v['arc_supported']) and len(unsupported)==int(v['unsupported_edges']))
                    replays+=1
    for row in summary['groups']:
        group=[v for v in rows if int(v['n'])==row['n'] and v['method']==row['method']]
        values=dict(samples=len(group),locality_failures=sum(not truth(v['locality_hypothesis']) for v in group),
            short_hops_failures=sum(not truth(v['short_hops_first']) for v in group),arc_support_failures=sum(not truth(v['arc_supported']) for v in group),
            multibar=sum(int(v['bars'])>1 for v in group),empty=sum(int(v['bars'])==0 for v in group),
            single_bar_born_at_B=sum(int(v['bars'])==int(v['bars_born_at_B'])==1 for v in group),positive_intrinsic=sum(int(v['intrinsic_bars'])==1 for v in group))
        require('peanut_summary_'+str((row['n'],row['method'])),all(row[k]==val for k,val in values.items()))
    require('tiny_positive_edge_case_present',min_chord<1e-8)

    # The two-curve locality samples have no chord below 1e-8, so a dense graph
    # conversion that drops near-zero weights cannot affect their saved results.
    locality_min=1.
    for kind,base_seed in [('ellipse',20350101),('peanut',20350201),('circle',20350301)]:
        for j in range(20):
            ss,xx=sample(kind,np.random.default_rng(base_seed+j));cc=np.linalg.norm(xx[:,None]-xx[None,:],axis=2)
            locality_min=min(locality_min,float(cc[np.triu_indices(len(ss),1)].min()))
    require('locality_edges_above_dense_null_threshold',locality_min>1e-8)

    for n in (4,8,23,128):
        for a in (1.25,1.5,2.):
            integral=n*a*quad(lambda u:u**(a-1)*(1-u)**(n-1),0,1,epsabs=1e-13)[0]
            gamma=math.exp(gammaln(a+1)+gammaln(n+1)-gammaln(n+a))
            require(f'radius_beta_integral_{n}_{a}',abs(integral-gamma)<2e-11)
    threshold=2/math.pi-.5
    require('heuristic_crossings',next(n for n in range(3,1000) if math.log(n)/n<threshold)==23 and next(n for n in range(3,1000) if 7*math.log(n)/n<threshold)==291)
    return dict(status='passed',extra_layer_rows=len(export),peanut_clouds=450,peanut_method_records=len(rows),
                independent_metric_diagram_replays=replays,smallest_peanut_trough_chord=min_chord,smallest_locality_check_chord=locality_min,
                radius_integral_checks=12,scope='Saved evidence and bounded independent replays; radius asymptotics are proved analytically in the paper.')
if __name__=='__main__':
    checks=[]
    def require(name,condition,detail=None):
        if not condition:raise AssertionError((name,detail))
        checks.append(name)
    report=verify_extra_layers_and_peanut(ROOT,require);report['checks_passed']=len(checks)
    (ROOT/'results/extra_layers_peanut_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
