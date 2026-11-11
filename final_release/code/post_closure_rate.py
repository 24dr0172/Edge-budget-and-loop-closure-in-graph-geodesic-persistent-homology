"""Diagram error after closure against the squared edge length (rate corollary).

Uniform circle samples, n in {64,128,256,512}, 20 seeds each. The n=512 clouds
are the uniform clouds of the clean-circle study (seeds 20310101-20310120);
the other sizes use seeds 20400000+100n+i. For each cloud the script records
the normalized diagram error e_PH, the largest minor-arc length l_G of the
selected graph and the bound (pi^2/6) l_G^2 b*/Delta of the rate corollary at

* radius order at exact closure,
* OR and mutual order at exact closure,
* the radius graph at twice the smallest connecting radius,

together with the Euclidean error and the cost P* of a missed loop. It also
checks the bound on equally spaced samples (a tiny deterministic jitter breaks
chord ties). Run `python3 code/post_closure_rate.py` to generate,
`--plot-only` to rebuild the figure from the saved records and `--verify` to
check the saved records, replay one cloud per size and confirm the bound.
"""
from pathlib import Path
import argparse,csv,hashlib,json,math,os,sys
sys.dont_write_bytecode=True
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_name]='1'
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import numpy as np
import graph_methods as m
from run_studies import measurements

SIZES=(64,128,256,512)
SEEDS=20
FIELDS=('n','seed','graph','M','g','b_star','delta','ell_G','tau_G','e_PH','bound','ratio',
        'euclid_error','P_star','ring_present')
BOUND=lambda ell,b,delta:math.pi**2/6*ell**2*b/delta

def seeds_for(n):
    if n==512:return [20310100+i for i in range(1,SEEDS+1)]
    return [20400000+100*n+i for i in range(1,SEEDS+1)]

def study(s,points,tag):
    """Return one record per selected graph of a single cloud."""
    n=len(s);true=m.intrinsic(s);reference=m.diagram(true);delta=float(true.max())
    assert reference.shape==(1,2)
    gaps=np.diff(np.r_[s,s[0]+1]);g=float(gaps.max());b_star=float(reference[0,1])
    assert abs(reference[0,0]-g)<1e-6  # Ripser returns float32 endpoints
    chord,iu,ju,length,orders,keys=m.build_orders(points)
    euclid,_=measurements(chord,true,reference)
    p_star=(b_star-g)/(2*delta)
    rows=[]
    def add(graph,order,M):
        D=chord if M==len(order) else m.graph_metric(order,iu,ju,length,M,n)
        edges=order[:M];score,bars=measurements(D,true,reference,maximum_edge=float(length[edges].max()))
        ell=float(true[iu[edges],ju[edges]].max())
        present=M>=m.closure(order,iu,ju,n)
        tau=float(bars[:,1].max()) if len(bars) else float('nan')
        bound=BOUND(ell,b_star,delta)
        if present:
            assert len(bars)==1,'single bar expected after closure'
            assert score['e_PH']<=bound*(1+1e-6)+1e-6/delta,(graph,score['e_PH'],bound)
        rows.append(dict(**tag,graph=graph,M=int(M),g=g,b_star=b_star,delta=delta,ell_G=ell,tau_G=tau,
            e_PH=score['e_PH'],bound=bound,ratio=score['e_PH']/bound,euclid_error=euclid['e_PH'],
            P_star=p_star,ring_present=bool(present)))
    for method,graph in (('R2_global','radius_closure'),('OR_kNN','or_closure'),('mutual_kNN','mutual_closure')):
        order=orders[method];add(graph,order,m.closure(order,iu,ju,n))
    order=orders['R2_global'];conn=m.first_connected(order,iu,ju,n)
    r_conn=float(length[order[conn-1]])
    M2=int(np.searchsorted(length[order],2*r_conn,side='right'))
    add('radius_double',order,M2)
    return rows

def generate(out):
    rows=[]
    for n in SIZES:
        for seed in seeds_for(n):
            s,points=m.sample_cloud(seed,'uniform',n=n)
            rows+=study(s,points,dict(n=n,seed=seed))
            print(n,seed,flush=True)
    # equally spaced samples with a tiny deterministic jitter that breaks chord ties
    even=[]
    for n in (64,128,256,512,1024):
        rng=np.random.default_rng(20410000+n)
        s=np.sort(np.arange(n)/n+rng.uniform(-1,1,n)*1e-7/n)%1.
        s=np.sort(s);points=np.column_stack((np.cos(2*np.pi*s),np.sin(2*np.pi*s)))/(2*np.pi)
        r=[v for v in study(s,points,dict(n=n,seed=0)) if v['graph']=='radius_closure'][0]
        even.append(r)
    out.mkdir(parents=True,exist_ok=True)
    with (out/'records.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(rows)
    with (out/'equally_spaced.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(even)
    return rows,even

def summarize(rows,even):
    summary={}
    for graph in ('radius_closure','or_closure','mutual_closure','radius_double'):
        sel=[r for r in rows if r['graph']==graph and r['ring_present']]
        summary[graph]=dict(records=sum(r['graph']==graph for r in rows),ring_present=len(sel),
            max_ratio=max(r['ratio'] for r in sel),median_ratio=float(np.median([r['ratio'] for r in sel])),
            max_e_PH=max(r['e_PH'] for r in sel))
    summary['euclidean_error_range']=[min(r['euclid_error'] for r in rows),max(r['euclid_error'] for r in rows)]
    summary['P_star_range']=[min(r['P_star'] for r in rows),max(r['P_star'] for r in rows)]
    summary['equally_spaced_ratio_by_n']={r['n']:r['ratio'] for r in even}
    return summary

def plot(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':300})
    fig,axes=plt.subplots(1,2,figsize=(6.9,3.0),layout='constrained')
    ok=lambda graph:[r for r in rows if r['graph']==graph and r['ring_present']]
    euclid=float(np.median([r['euclid_error'] for r in rows]));pstar=float(np.median([r['P_star'] for r in rows]))
    ax=axes[0]
    xs=np.geomspace(.004,.2,50)
    ax.plot(xs,math.pi**2/9*xs**2,'k-',lw=1,label=r'$(\pi^2/9)g^2$')
    ax.plot(xs,math.pi**2/9*(2*xs)**2,'k--',lw=1,label=r'$(\pi^2/9)(2g)^2$')
    ax.axhline(euclid,color='0.45',lw=1,ls=':');ax.axhline(pstar,color='0.45',lw=1,ls='-.')
    ax.text(.0042,euclid*1.12,'Euclidean',fontsize=8,color='0.3');ax.text(.0042,pstar*1.08,'empty diagram',fontsize=8,color='0.3')
    for graph,color,mark,label in (('radius_closure','#0072B2','o',r'$\varepsilon=B$'),('radius_double','#D55E00','s',r'$\varepsilon=2\varepsilon_{\rm conn}$')):
        ax.scatter([r['g'] for r in ok(graph)],[r['e_PH'] for r in ok(graph)],s=14,color=color,marker=mark,label=label,alpha=.8,lw=0)
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel(r'Largest gap $g$');ax.set_ylabel(r'Diagram error $e_{\rm PH}$')
    ax.set_title(r'(a) $\varepsilon$-neighborhood graphs',fontsize=10);ax.legend(frameon=False,fontsize=7.5,loc='lower right')
    ax=axes[1]
    xs=np.geomspace(.004,.2,50);ax.plot(xs,math.pi**2/9*xs**2,'k-',lw=1,label=r'$(\pi^2/9)\ell_G^2$')
    for graph,color,mark,label in (('radius_closure','#0072B2','o',r'$\varepsilon$-neighborhood'),('or_closure','#D55E00','s','symmetric kNN'),('mutual_closure','#009E73','^','mutual kNN')):
        ax.scatter([r['ell_G'] for r in ok(graph)],[r['e_PH'] for r in ok(graph)],s=14,color=color,marker=mark,label=label,alpha=.8,lw=0)
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel(r'Longest edge arc $\ell_G$');ax.set_title('(b) Cycle closure',fontsize=10)
    ax.legend(frameon=False,fontsize=7.5,loc='lower right')
    for ax in axes:ax.grid(alpha=.2,which='major')
    (ROOT/'results/figure_exports').mkdir(parents=True,exist_ok=True);(ROOT/'manuscript_figures').mkdir(exist_ok=True)
    fig.savefig(ROOT/'results/figure_exports/post_closure_rate.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(ROOT/'manuscript_figures/post_closure_rate.png');plt.close(fig)


def verify(root,require):
    """Check the saved records of Figure 3 and replay the first cloud of each size."""
    root=Path(root);out=root/'results/post_closure_rate'
    rows=[{k:(v if k in ('graph','ring_present') else float(v)) for k,v in r.items()} for r in csv.DictReader((out/'records.csv').open())]
    for r in rows:r['ring_present']=r['ring_present']=='True'
    even=[{k:(v if k in ('graph','ring_present') else float(v)) for k,v in r.items()} for r in csv.DictReader((out/'equally_spaced.csv').open())]
    config=json.loads((out/'config.json').read_text())
    require('rate_records_hash',config['records_sha256']==hashlib.sha256((out/'records.csv').read_bytes()).hexdigest())
    graphs=('radius_closure','or_closure','mutual_closure','radius_double')
    require('rate_scope',len(rows)==len(SIZES)*SEEDS*len(graphs) and
            all(sum(r['n']==n and r['graph']==g for r in rows)==SEEDS for n in SIZES for g in graphs) and
            [r['n'] for r in even]==[64,128,256,512,1024])
    require('rate_seeds',all(r['seed'] in seeds_for(int(r['n'])) for r in rows))
    # The bound of the rate corollary holds for every graph that contains the ring.
    present=[r for r in rows if r['ring_present']]
    require('rate_closure_present',all(r['ring_present'] for r in rows if r['graph']!='radius_double') and len(present)>=len(rows)-1)
    slack=[r['bound']*(1+1e-6)+1e-6/r['delta']-r['e_PH'] for r in present]
    require('rate_bound_holds',min(slack)>=0,{'graphs_checked':len(present),'maximum_ratio':max(r['ratio'] for r in present)})
    require('rate_euclidean_and_missed_loop',all(r['e_PH']<r['euclid_error']<r['P_star'] for r in present),
            {'euclidean_error_range':[min(r['euclid_error'] for r in rows),max(r['euclid_error'] for r in rows)]})
    # Equally spaced samples approach equality, up to the float32 endpoints returned by Ripser.
    require('rate_equally_spaced',all(abs(r['ratio']-1)<5e-3 for r in even if r['n']<=512) and abs(even[-1]['ratio']-1)<5e-2,
            {str(int(r['n'])):r['ratio'] for r in even})
    summary=json.loads((out/'summary.json').read_text());fresh=summarize(rows,even)
    require('rate_summary_replay',json.dumps(fresh,sort_keys=True,default=float)==json.dumps(summary,sort_keys=True,default=float) or
            all(abs(fresh[g]['max_ratio']-summary[g]['max_ratio'])<1e-12 for g in graphs))
    replayed=0
    for n in SIZES:
        seed=seeds_for(n)[0];s,points=m.sample_cloud(seed,'uniform',n=n)
        for got in study(s,points,dict(n=n,seed=seed)):
            saved=next(r for r in rows if r['n']==n and r['seed']==seed and r['graph']==got['graph'])
            require(f"rate_replay_{n}_{got['graph']}",saved['M']==got['M'] and
                    all(abs(saved[k]-got[k])<1e-9 for k in ('g','b_star','delta','ell_G','e_PH','bound','euclid_error','P_star')) and
                    saved['ring_present']==got['ring_present'])
            replayed+=1
    require('rate_figure_files',(root/'manuscript_figures/post_closure_rate.png').is_file() and
            (root/'results/figure_exports/post_closure_rate.pdf').is_file())
    return dict(status='passed',records=len(rows),replayed_graphs=replayed,maximum_ratio=max(r['ratio'] for r in present))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plot-only',action='store_true');parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();out=ROOT/'results/post_closure_rate'
    if args.verify:
        checks=[]
        def require(name,condition,detail=None):
            if not condition:raise AssertionError((name,detail))
            checks.append(name)
        print(json.dumps(verify(ROOT,require),indent=2));return
    if args.plot_only:
        rows=[{k:(v if k in ('graph','ring_present') else float(v)) for k,v in r.items()} for r in csv.DictReader((out/'records.csv').open())]
        for r in rows:r['ring_present']=r['ring_present']=='True'
        even=[{k:(v if k in ('graph','ring_present') else float(v)) for k,v in r.items()} for r in csv.DictReader((out/'equally_spaced.csv').open())]
    else:rows,even=generate(out)
    summary=summarize(rows,even);plot(rows)
    config=dict(sizes=list(SIZES),seeds_per_size=SEEDS,n512_seeds='20310101-20310120 (uniform clouds of the clean-circle study)',
        other_seeds='20400000 + 100n + i, i=1..20',graphs=['radius_closure','or_closure','mutual_closure','radius_double'],
        bound='(pi^2/6) l_G^2 b*/Delta',records_sha256=hashlib.sha256((out/'records.csv').read_bytes()).hexdigest(),
        equally_spaced='n in 64..1024, jitter 1e-7/n, radius order at closure')
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
