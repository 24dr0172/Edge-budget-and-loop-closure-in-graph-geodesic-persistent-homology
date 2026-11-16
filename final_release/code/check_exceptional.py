"""Full diagrams for every saved sample outside the estimated good event."""
from pathlib import Path
import concurrent.futures, csv, hashlib, json, sys
import numpy as np

RELEASE=Path(__file__).resolve().parents[1]
HERE=RELEASE/'results/smooth_curve_checks'
sys.path.insert(0,str(RELEASE/'code'))
import smooth_curve_enlargement as smooth

def check(job):
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import minimum_spanning_tree, shortest_path, connected_components
    from ripser import ripser
    curve,index,cutoff=job
    start,gaps=smooth.cloud(curve,128,index)
    s=(start+np.r_[0.,np.cumsum(gaps)[:-1]])%1
    xy=smooth.points(curve,s)
    chord=np.linalg.norm(xy[:,None,:]-xy[None,:,:],axis=2)
    consecutive=np.sort(chord[np.arange(128),(np.arange(128)+1)%128])
    B,second=map(float,consecutive[-2:][::-1])
    rconn=float(minimum_spanning_tree(chord).data.max())
    arrays=dict(s=s,xy=xy,gaps=gaps)
    rows=[]
    for alpha in smooth.ALPHAS:
        graph=csr_matrix(np.where(chord<=alpha*rconn,chord,0.))
        assert connected_components(graph,directed=False,return_labels=False)==1
        D=shortest_path(graph,directed=False,method='D');D=(D+D.T)/2;np.fill_diagonal(D,0.)
        bars=ripser(D,distance_matrix=True,maxdim=1)['dgms'][1]
        bars=bars[bars[:,1]>bars[:,0]]
        arrays[f'bars_{alpha:g}']=bars
        actual_miss=len(bars)==0
        rows.append(dict(curve=curve,n=128,cloud=index,alpha=alpha,largest_gap=float(gaps.max()),
            estimated_cutoff=cutoff,B=B,second_chord=second,r_conn=rconn,
            r_conn_is_second_chord=bool(rconn==second),table_proxy_miss=bool(B>alpha*second),
            mst_chord_miss=bool(B>alpha*rconn),h1_bars=len(bars),diagram_miss=actual_miss,
            table_proxy_agrees=bool(actual_miss==(B>alpha*second))))
    np.savez_compressed(HERE/'exceptional_diagrams'/f'{curve}_{index}.npz',**arrays)
    return rows

def main():
    cfg=json.loads((RELEASE/'results/smooth_curve_enlargement/curve_constants.json').read_text())
    batches=list(csv.DictReader((RELEASE/'results/smooth_curve_enlargement/batches.csv').open()))
    jobs=[]
    for row in batches:
        if int(row['n'])!=128 or int(row['good_event'])==int(row['samples']):continue
        curve=row['curve'];batch=int(row['batch']);cutoff=cfg[curve]['good_event_threshold']
        start,gaps=smooth.draw(smooth.batch_rng(curve,128,batch),smooth.BATCH,128)
        bad=np.flatnonzero(gaps.max(axis=1)>=cutoff)
        assert len(bad)==int(row['samples'])-int(row['good_event'])
        jobs.extend((curve,batch*smooth.BATCH+int(i),cutoff) for i in bad)
    assert sum(j[0]=='ellipse' for j in jobs)==2
    assert sum(j[0]=='peanut' for j in jobs)==39
    (HERE/'exceptional_diagrams').mkdir(exist_ok=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        groups=list(pool.map(check,jobs))
    rows=[r for group in groups for r in group]
    with (HERE/'exceptional_diagrams.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    script_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report=dict(point_sets=len(jobs),diagrams=len(rows),table_proxy_disagreements=sum(not r['table_proxy_agrees'] for r in rows),
        connectivity_threshold_disagreements=sum(not r['r_conn_is_second_chord'] for r in rows),
        circle_equivalence_not_assumed=True,geometric_thresholds='Grid estimates, as in the original study',
        per_curve={c:dict(point_sets=sum(j[0]==c for j in jobs),
            nonempty=sum(not r['diagram_miss'] for r in rows if r['curve']==c),
            empty=sum(r['diagram_miss'] for r in rows if r['curve']==c)) for c in smooth.CURVES},
        generation_script_sha256=script_hash,verification_script_sha256=script_hash)
    (HERE/'exceptional_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
