"""Fresh serial timings for representative full H1-diagram grid evaluations."""
from pathlib import Path
import argparse, csv, json, os, platform, time
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import sys
sys.dont_write_bytecode=True
import numpy as np
import paper1_methods as m
from run_studies import finite_graph_diagram

def benchmark(root,repeats=3):
    root=Path(root)
    with (root/'results/clean_closure/budget_summary.csv').open() as f:
        budgets=sorted({int(row['M']) for row in csv.DictReader(f)})
    cases=[]
    s,clean=m.sample_cloud(20310101,'uniform');cases.append(('clean_uniform',20310101,clean,False))
    for model in ['radial','isotropic']:
        with np.load(root/f'results/noise/checkpoints/{model}_20320101_delta0.01.npz') as z:
            cases.append((model+'_delta0.01',20320101,z['observed'].copy(),True))
    with np.load(root/'results/ellipse_checks/coordinates.npz') as z:
        cases.append(('ellipse_uniform',20310101,z['20310101_uniform_xy'].copy(),False))
    m.diagram(np.array([[0.,1.,1.],[1.,0.,1.],[1.,1.,0.]]))
    measurements=[]
    for name,seed,points,adaptive in cases:
        for repetition in range(1,repeats+1):
            start=time.perf_counter()
            chord,iu,ju,length,orders,keys=m.build_orders(points)
            order_seconds=time.perf_counter()-start
            apsp=ph=0.;connected=0;full=0;calls=[]
            for method in m.METHODS:
                order=orders[method]
                for budget in budgets:
                    t=time.perf_counter();D=m.graph_metric(order,iu,ju,length,budget,len(points))
                    sp_seconds=time.perf_counter()-t;apsp+=sp_seconds
                    if not np.isfinite(D).all():continue
                    connected+=1;t=time.perf_counter()
                    if adaptive:
                        bars=finite_graph_diagram(D,float(length[order[:budget]].max()))
                    else:bars=m.diagram(D)
                    diagram_seconds=time.perf_counter()-t;ph+=diagram_seconds;full+=1
                    calls.append(dict(method=method,M=budget,apsp_seconds=sp_seconds,
                        diagram_seconds=diagram_seconds,positive_bars=len(bars)))
            elapsed=time.perf_counter()-start
            row=dict(case=name,seed=seed,repetition=repetition,n=len(points),grid_budgets=len(budgets),
                methods=3,graph_evaluations=3*len(budgets),connected_graphs=connected,
                full_H1_diagrams=full,adaptive_certified_cutoff=adaptive,
                edge_order_seconds=order_seconds,apsp_seconds=apsp,diagram_seconds=ph,
                total_seconds=elapsed,per_diagram_calls=calls)
            measurements.append(row)
            print(name,repetition,'seconds',round(elapsed,3),'full H1 diagrams',full,flush=True)
    cpu=platform.processor()
    if Path('/proc/cpuinfo').exists():
        cpu=next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),cpu)
    return dict(scope='Fresh representative timings, not total study generation times.',
        cpu=cpu,platform=platform.platform(),python=platform.python_version(),numpy=np.__version__,
        scipy_version=__import__('scipy').__version__,ripser_version=__import__('ripser').__version__,
        process_mode='Serial process; OPENBLAS/OMP/MKL thread limits set to 1 before imports.',
        warmup='One three-vertex Ripser call before timing; no disk writes in timed blocks.',
        included='Construct all pair/rank orders; 24 budgets for each of three graph rules; all-pairs Dijkstra; full finite H1 diagrams of connected graphs.',
        excluded='Sampling, intrinsic references, bottleneck errors, theorem certificates, figure generation and file I/O.',
        largest_full_diagram_n_in_studies=512,largest_rank_only_n=1000000,
        repetitions=repeats,measurements=measurements)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--output',type=Path);ap.add_argument('--repeats',type=int,default=3);a=ap.parse_args()
    result=benchmark(a.release,a.repeats)
    (a.output or a.release/'results/runtime_benchmark.json').write_text(json.dumps(result,indent=2)+'\n')
