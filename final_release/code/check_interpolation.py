"""Compare the released arc maps with refined Simpson/cubic inverse maps.

This is numerical evidence, not an interval-arithmetic error certificate.
The largest eight gaps are the candidates used in the released sampler.
"""
from pathlib import Path
import csv, hashlib, json, sys
import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.interpolate import CubicSpline

RELEASE=Path(__file__).resolve().parents[1]
HERE=RELEASE/'results/smooth_curve_checks'
sys.path.insert(0,str(RELEASE/'code'))
import smooth_curve_enlargement as smooth

def reference(curve, count):
    t=np.linspace(0.,2*np.pi,count)
    if curve=='ellipse':speed=np.sqrt(4*np.sin(t)**2+np.cos(t)**2)
    else:speed=np.sqrt((1+.4*np.cos(2*t))**2+.64*np.sin(2*t)**2)
    arc=cumulative_simpson(speed,x=t,initial=0.)
    L=float(arc[-1]);inverse=CubicSpline(arc/L,t)
    def evaluate(s):
        theta=inverse(np.asarray(s)%1.)
        if curve=='ellipse':return np.stack((2*np.cos(theta),np.sin(theta)),axis=-1)/L
        r=1+.4*np.cos(2*theta)
        return np.stack((r*np.cos(theta),r*np.sin(theta)),axis=-1)/L
    return evaluate,L

def main():
    reports=[];records=[]
    for curve in smooth.CURVES:
        ref,L=reference(curve,1600001)
        fine,Lfine=reference(curve,3200001)
        # Phase-shifted arc grid avoids checking only original interpolation nodes.
        probe=(np.arange(200001)+.37)/200001
        uniform_error=float(np.linalg.norm(smooth.points(curve,probe)-fine(probe),axis=-1).max())
        reference_difference=float(np.linalg.norm(ref(probe)-fine(probe),axis=-1).max())
        chord_errors=[];corrections=[];remainder=[];class_changes=0;relative=[];unresolved=0
        rng=np.random.default_rng(2026111502+smooth.CURVE_STREAMS[curve])
        for block in range(20):
            start,gaps=smooth.draw(rng,50,32768)
            indices=np.argpartition(gaps,32768-8,axis=1)[:,-8:]
            left=np.cumsum(gaps,axis=1)-gaps
            starts=start[:,None]+np.take_along_axis(left,indices,axis=1)
            lengths=np.take_along_axis(gaps,indices,axis=1);ends=starts+lengths
            old_a,old_b=smooth.points(curve,starts),smooth.points(curve,ends)
            new_a,new_b=fine(starts),fine(ends)
            old_c=np.linalg.norm(old_b-old_a,axis=-1);new_c=np.linalg.norm(new_b-new_a,axis=-1)
            err=np.abs(old_c-new_c);correction=lengths-new_c
            chord_errors.extend(err.ravel());corrections.extend(correction.ravel())
            ref_c=np.linalg.norm(ref(ends)-ref(starts),axis=-1)
            remainder.extend(np.abs(ref_c-new_c).ravel())
            old_sort=np.sort(old_c,axis=1);new_sort=np.sort(new_c,axis=1)
            for alpha in smooth.ALPHAS:
                class_changes+=int(np.sum((old_sort[:,-1]>alpha*old_sort[:,-2])!=(new_sort[:,-1]>alpha*new_sort[:,-2])))
            # Report the two chord-determining gaps per cloud separately.
            sel=np.argsort(new_c,axis=1)[:,-2:]
            two_e=np.take_along_axis(err,sel,axis=1);two_d=np.take_along_axis(correction,sel,axis=1)
            # Avoid meaningless ratios when the reference correction is near
            # floating-point/refinement uncertainty; absolute errors remain saved.
            resolved=two_d>1e-13
            unresolved+=int(np.sum(~resolved))
            relative.extend((two_e[resolved]/two_d[resolved]).ravel())
            for j in range(50):
                for k in sel[j]:
                    records.append(dict(curve=curve,cloud=block*50+j,gap=float(lengths[j,k]),
                        reference_chord=float(new_c[j,k]),chord_error=float(err[j,k]),
                        chord_arc_difference=float(correction[j,k])))
        e=np.asarray(chord_errors);d=np.asarray(corrections);q=np.asarray(relative)
        reports.append(dict(curve=curve,clouds=1000,n=32768,candidate_chords=len(e),
            original_arc_grid=400001,reference_grids=[1600001,3200001],reference_method='Exact speed; cumulative Simpson integration; cubic inverse interpolation',
            probe_points=len(probe),max_coordinate_error=uniform_error,reference_coordinate_difference=reference_difference,
            max_chord_error=float(e.max()),reference_chord_difference=float(max(remainder)),
            chord_arc_difference_quantiles=dict(zip(['minimum','q01','median','q99','maximum'],map(float,np.quantile(d,[0,.01,.5,.99,1])))),
            top_two_error_over_correction_quantiles=dict(zip(['minimum','median','q99','maximum'],map(float,np.quantile(q,[0,.5,.99,1])))),
            top_two_error_below_correction=int(np.sum(q<1)),top_two_resolved_pairs=len(q),
            top_two_unresolved_corrections=unresolved,correction_resolution_cutoff=1e-13,
            changed_enlargement_classifications=class_changes,comparisons=3000,
            perimeter=Lfine,reference_perimeter_difference=abs(L-Lfine)))
        print(json.dumps(reports[-1],indent=2),flush=True)
    with (HERE/'interpolation_pairs.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    script_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'interpolation.json').write_text(json.dumps(dict(protocol='1000 new independent uniform point sets per curve; all eight largest candidate gaps; same three enlargement factors',
        seed_rule='2026111502+curve_stream',numerical_not_certified=True,
        generation_script_sha256=script_hash,verification_script_sha256=script_hash,
        results=reports),indent=2)+'\n')

if __name__=='__main__':main()
