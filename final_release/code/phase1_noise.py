"""Compute Euclidean Rips baselines and paired selection-rule comparisons."""
from pathlib import Path
import argparse, hashlib, json, os, resource, sys, time
sys.dont_write_bytecode = True
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from run_studies import one_bar_distance
from paper1_methods import intrinsic, diagram

F = Path(__file__).resolve().parents[1]
OUT = F / 'results/phase1/euclidean'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(limit=None):
    OUT.mkdir(parents=True, exist_ok=True)
    checkpoints = OUT / 'checkpoints'
    checkpoints.mkdir(exist_ok=True)
    inputs = sorted((F / 'results/noise/checkpoints').glob('*.npz'))
    start = time.monotonic()
    measured = []
    for index, path in enumerate(inputs):
        if limit is not None and index >= limit:
            break
        target = checkpoints / path.name
        digest = sha(path)
        if target.exists():
            with np.load(target) as z:
                assert str(z['source_sha256']) == digest
            continue
        with np.load(path) as z:
            xy, latent, reference = z['observed'], z['latent'], z['reference']
            assert len(xy) == 512 and reference.shape == (1, 2)
            t = time.monotonic()
            distances = cdist(xy, xy)
            bars = diagram(distances)
            assert np.isfinite(bars).all() and (bars[:, 1] > bars[:, 0]).all()
            diameter = float(intrinsic(latent).max())
            error = one_bar_distance(bars, reference) / diameter
            elapsed = time.monotonic() - t
        np.savez_compressed(target, bars=bars, reference=reference,
                            diameter=diameter, e_PH=error,
                            elapsed_seconds=elapsed, source_sha256=digest)
        measured.append(elapsed)
        if (index + 1) % 10 == 0 or index == 0:
            print(f'Euclidean {index+1}/{len(inputs)}: {elapsed:.3f}s, error {error:.6g}', flush=True)
    timing = dict(status='complete' if limit is None else 'pilot',
                  distinct_input_clouds=len(inputs),
                  checkpoints_completed=len(list(checkpoints.glob('*.npz'))),
                  seconds_this_invocation=time.monotonic()-start,
                  new_diagram_seconds=measured,
                  peak_process_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  scope='Full unrestricted finite H1. Existing checkpoints are reused; individual times exclude input loading and checkpoint writing.',
                  newly_computed_diagrams=len(measured))
    (OUT / ('timing.json' if limit is None else 'pilot_timing.json')).write_text(json.dumps(timing, indent=2)+'\n')
    if limit is None:
        summarize()


def summarize():
    rows = []
    cfg = json.loads((F / 'code/config.json').read_text())
    for seed in cfg['noise']['seeds']:
        for model in ['radial', 'isotropic']:
            for delta in cfg['noise']['amplitudes']:
                stem = f'{"radial" if delta == 0 else model}_{seed}_delta{delta:g}'
                path = OUT / 'checkpoints' / (stem + '.npz')
                with np.load(path) as z:
                    bars=z['bars']
                    rows.append(dict(seed=seed, model=model, delta=delta,
                                     baseline='Euclidean', e_PH=float(z['e_PH']),
                                     bars=len(bars), diagram_file=str(path.relative_to(F)),
                                     reused_clean=bool(delta == 0 and model == 'isotropic')))
    raw=pd.DataFrame(rows)
    raw.to_csv(OUT / 'samples.csv', index=False)
    group=raw.groupby(['model', 'delta']).e_PH
    summary=group.agg(samples='size', median_e_PH='median').reset_index()
    summary['q1_e_PH']=group.quantile(.25, interpolation='linear').values
    summary['q3_e_PH']=group.quantile(.75, interpolation='linear').values
    summary.to_csv(OUT / 'summary.csv', index=False)
    selectors=pd.read_csv(F / 'results/noise/selectors.csv')
    paired=selectors.merge(raw[['seed','model','delta','e_PH']], on=['seed','model','delta'],
                           suffixes=('_selected', '_euclidean'), validate='many_to_one')
    paired['selected_minus_euclidean']=paired.e_PH_selected-paired.e_PH_euclidean
    paired.to_csv(OUT / 'paired_selectors.csv', index=False)
    summaries=[]
    for key,g in paired.groupby(['model','delta','method','selector']):
        differences=g.selected_minus_euclidean
        summaries.append(dict(zip(['model','delta','method','selector'],key),
                              samples=len(g), median_selected_e_PH=g.e_PH_selected.median(),
                              median_euclidean_e_PH=g.e_PH_euclidean.median(),
                              median_paired_difference=differences.median(),
                              selected_better=int((differences<0).sum()),
                              euclidean_better=int((differences>0).sum()),
                              equal=int((differences==0).sum())))
    pd.DataFrame(summaries).to_csv(OUT/'paired_summary.csv',index=False)
    print('Euclidean baseline saved: 220 distinct diagrams; 240 model/amplitude/seed rows.', flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--pilot',type=int)
    parser.add_argument('--summarize-only',action='store_true')
    args=parser.parse_args()
    summarize() if args.summarize_only else run(args.pilot)
