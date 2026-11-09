"""Independent million-cloud circle checks for the two rare alpha=2 rows.

For iid exponentials, the second-largest CDF value is Beta(n-1, 2), the
largest exceeds it by an independent Exp(1), and the other n-2 values are
independent exponentials truncated below the second-largest value. Conditional
sampling lets us omit values that cannot affect the chord-miss decision.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SIZES = (4096, 32768)
REPETITIONS = 1_000_000
BATCH = 8192
INITIAL = 64
SEED = 20261104


def batch(rng, n, m):
    # F(Y) for the second-largest exponential has the stated beta law.
    y = -np.log1p(-rng.beta(n-1, 2, size=m))
    x = y + rng.exponential(size=m)
    cap = -np.expm1(-y)
    values = -np.log1p(-rng.random((m, INITIAL))*cap[:,None])
    partial = values.sum(axis=1)
    lower = x + y + partial
    # Additional independent truncated values certify g_max < 1/3. A
    # subset of the n-2 values already drawn is enough whenever lower > 3x.
    used = np.full(m, INITIAL, dtype=np.int32)
    while True:
        ix = np.flatnonzero(lower <= 3*x)
        if not len(ix): break
        if np.any(used[ix] >= n-2):
            raise RuntimeError('A full cloud has largest gap at least 1/3')
        t = min(INITIAL, int(np.min(n-2-used[ix])))
        extra = -np.log1p(-rng.random((len(ix),t))*cap[ix,None])
        lower[ix] += extra.sum(axis=1)
        used[ix] += t
    # Concavity makes x <= 2y a certain non-miss. Draw all remaining
    # truncated exponentials only for the candidates x > 2y.
    ix = np.flatnonzero(x > 2*y)
    exact = lower[ix].copy()
    for used_count in np.unique(used[ix]):
        sel = np.flatnonzero(used[ix] == used_count)
        count = n-2-int(used_count)
        for start in range(0, count, 1024):
            k=min(1024,count-start)
            extra=-np.log1p(-rng.random((len(sel),k))*cap[ix[sel],None])
            exact[sel]+=extra.sum(axis=1)
    flags = np.sin(np.pi*x[ix]/exact) > 2*np.sin(np.pi*y[ix]/exact)
    return x,y,lower,ix,exact,flags


def generate():
    folder=ROOT/'results/radius_enlargement_extra';folder.mkdir(parents=True,exist_ok=True)
    rows=[]
    for n in SIZES:
        rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED,n])))
        records=[]; misses=0; candidates=0; min_margin=float('inf')
        for start in range(0,REPETITIONS,BATCH):
            m=min(BATCH,REPETITIONS-start)
            x,y,lower,ix,exact,flags=batch(rng,n,m)
            margin=float(np.min(lower-3*x));min_margin=min(min_margin,margin)
            misses+=int(flags.sum());candidates+=len(ix)
            records.append(dict(start=start,samples=m,spacing_candidates=len(ix),
                                chord_misses=int(flags.sum()),min_good_event_margin=margin))
        rows.append(dict(n=n,alpha=2,samples=REPETITIONS,misses=misses,
                         miss_rate=misses/REPETITIONS,spacing_misses=candidates,
                         min_good_event_margin=min_margin,seed=SEED))
        (folder/f'batches_{n}.json').write_text(json.dumps(records,indent=2)+'\n')
        print(f'n={n}: misses={misses}, spacing candidates={candidates}, samples={REPETITIONS}',flush=True)
    with (folder/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (folder/'config.json').write_text(json.dumps(dict(seed=SEED,sizes=list(SIZES),samples=REPETITIONS,
        batch=BATCH,initial=INITIAL,method='Exact conditional exponential order-statistic sampling; partial sums certify g < 1/3; all remaining gaps drawn for spacing-event candidates.',
        generator='PCG64 with SeedSequence([seed, n])',independent_of_original_cohort=True,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')

def verify(root, require):
    folder=Path(root)/'results/radius_enlargement_extra'
    config=json.loads((folder/'config.json').read_text())
    require('extra_circle_protocol', config['samples']==REPETITIONS and config['sizes']==list(SIZES)
            and config['seed']==SEED and config['batch']==BATCH and config['initial']==INITIAL)
    require('extra_circle_script_hash',config['script_sha256']==hashlib.sha256(
        (Path(root)/'code/radius_enlargement_extra.py').read_bytes()).hexdigest())
    with (folder/'summary.csv').open(newline='') as f: summaries={int(row['n']):row for row in csv.DictReader(f)}
    require('extra_circle_rows',set(summaries)==set(SIZES))
    for n in SIZES:
        records=json.loads((folder/f'batches_{n}.json').read_text())
        require(f'extra_circle_batch_count_{n}',len(records)==(REPETITIONS+BATCH-1)//BATCH)
        rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED,n])))
        total=misses=candidates=0; minimum=float('inf')
        for record in records:
            m=min(BATCH,REPETITIONS-total)
            x,y,lower,ix,exact,flags=batch(rng,n,m)
            margin=float(np.min(lower-3*x)); minimum=min(minimum,margin)
            require(f'extra_circle_batch_{n}_{total}',record['start']==total and record['samples']==m
                and record['spacing_candidates']==len(ix) and record['chord_misses']==int(flags.sum())
                and record['min_good_event_margin']==margin and margin>0
                and np.all(y<=x) and np.all(lower>3*x) and np.all(exact>=lower[ix]))
            total+=m; misses+=int(flags.sum());candidates+=len(ix)
        row=summaries[n]
        require(f'extra_circle_summary_{n}',total==REPETITIONS and int(row['samples'])==total
            and int(row['misses'])==misses and int(row['spacing_misses'])==candidates
            and float(row['miss_rate'])==misses/total and float(row['min_good_event_margin'])==minimum)
    return {'samples_per_size':REPETITIONS,'sizes':list(SIZES),'status':'passed'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if args.verify:
        def require(name, condition, detail=None):
            if not condition: raise AssertionError(name)
        print(json.dumps(verify(ROOT,require)))
    else: generate()
