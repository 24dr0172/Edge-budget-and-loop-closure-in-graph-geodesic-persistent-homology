"""Generate or verify the finite-size radius-enlargement study of the finite-size appendix.

Normalized independent exponential variables are exact uniform circular
spacings. The circle event theorem then decides retention without computing
a persistence diagram. The three enlargements share each simulated cloud.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SIZES = (128, 512, 4096, 32768)
EPSILONS = (.25, .5, 1.)
REPETITIONS = 100000
MASTER_SEED = 20261024
BATCH_SIZE = 128
FIELDS = ('n', 'epsilon', 'samples', 'misses', 'miss_rate', 'wilson_lo',
          'wilson_hi', 'spacing_misses', 'spacing_probability',
          'asymptotic_probability')


def wilson(k, n):
    z = 1.959963984540054
    p = k / n
    den = 1 + z*z/n
    centre = (p + z*z/(2*n)) / den
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return centre-half, centre+half


def predictions(n, epsilon):
    # Exact probability for the ARC-gap ratio, not the chord-distance event.
    exact_spacing = math.exp(math.lgamma(epsilon+2) + math.lgamma(n+1)
                             - math.lgamma(n+1+epsilon))
    asymptotic = math.gamma(epsilon+2) * n**(-epsilon)
    return exact_spacing, asymptotic


def summarize(arrays):
    rows = []
    for epsilon in EPSILONS:
        for n in SIZES:
            largest, second = arrays[f'n{n}_largest'], arrays[f'n{n}_second']
            failed = np.sin(np.pi*largest) > (1+epsilon)*np.sin(np.pi*second)
            k = int(failed.sum())
            lo, hi = wilson(k, len(largest))
            spacing, asymptotic = predictions(n, epsilon)
            rows.append(dict(zip(FIELDS, (n, epsilon, len(largest), k,
                k/len(largest), lo, hi, int((largest > (1+epsilon)*second).sum()),
                spacing, asymptotic))))
    return rows


def generate(root=ROOT):
    folder = Path(root)/'results/radius_enlargement'
    folder.mkdir(parents=True, exist_ok=True)
    arrays = {}
    witnesses = {}
    for n in SIZES:
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([MASTER_SEED, n])))
        largest, second = np.empty(REPETITIONS), np.empty(REPETITIONS)
        for start in range(0, REPETITIONS, BATCH_SIZE):
            count = min(BATCH_SIZE, REPETITIONS-start)
            exp = rng.exponential(size=(count, n))
            total = exp.sum(axis=1)
            original = exp.copy() if n == SIZES[0] else None
            exp.partition(n-2, axis=1)
            x, y = exp[:, -1]/total, exp[:, -2]/total
            largest[start:start+count], second[start:start+count] = x, y
            if original is not None:
                for epsilon in EPSILONS:
                    failed = np.sin(np.pi*x) > (1+epsilon)*np.sin(np.pi*y)
                    for outcome in (False, True):
                        key = f'witness_eps{epsilon:g}_miss{int(outcome)}'
                        candidates = np.flatnonzero(failed == outcome)
                        if key not in witnesses and len(candidates):
                            i = candidates[0]
                            witnesses[key] = original[i]/total[i]
        if not np.all(largest < 1/3):
            raise ValueError('A large-gap cloud needs a full diagram; no cloud has been discarded.')
        arrays[f'n{n}_largest'], arrays[f'n{n}_second'] = largest, second
        print(f'n={n}: {REPETITIONS:,} clouds; largest observed gap {largest.max():.6g}', flush=True)
    assert len(witnesses) == 2*len(EPSILONS)
    np.savez_compressed(folder/'gap_extremes.npz', **arrays, **witnesses)
    rows = summarize(arrays)
    with (folder/'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader(); writer.writerows(rows)
    config = dict(sizes=list(SIZES), epsilons=list(EPSILONS), samples_per_size=REPETITIONS,
        master_seed=MASTER_SEED, seed_sequence='[master_seed, n]', generator='PCG64',
        batch_size=BATCH_SIZE, sampling='n independent Exp(1) variables divided by their sum',
        miss_event='sin(pi*g_largest) > (1+epsilon)*sin(pi*g_second)',
        interval='Pointwise two-sided 95% Wilson; epsilon values are paired within each cloud.',
        all_gaps_below_one_third=True, discarded_clouds=0, new_clouds=len(SIZES)*REPETITIONS,
        persistence_diagrams_generated=0,
        spacing_benchmark='Gamma(epsilon+2)*Gamma(n+1)/Gamma(n+1+epsilon); an arc-gap event, not the chord event',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (folder/'config.json').write_text(json.dumps(config, indent=2)+'\n')
    return rows


def table_row(row):
    return (f"${row['epsilon']:g}$ & ${row['n']:,}$ & ${row['misses']:,}$ & "
            f"${100*row['miss_rate']:.4f}$ & "
            f"$[{100*row['wilson_lo']:.4f}, {100*row['wilson_hi']:.4f}]$ & "
            f"${100*row['spacing_probability']:.4f}$ & "
            f"${100*row['asymptotic_probability']:.4f}$")


def verify(root, require, paper_text=None):
    """Recount saved gaps, replay seeded first batches and check graph witnesses."""
    from scipy.integrate import quad
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components, shortest_path
    from ripser import ripser
    folder = Path(root)/'results/radius_enlargement'
    config = json.loads((folder/'config.json').read_text())
    require('radius_check_protocol', config['sizes']==list(SIZES) and config['epsilons']==list(EPSILONS)
        and config['samples_per_size']==REPETITIONS and config['master_seed']==MASTER_SEED
        and config['new_clouds']==400000 and config['discarded_clouds']==0)
    require('radius_check_script_hash', config.get('verification_script_sha256',config['script_sha256'])==hashlib.sha256(
        (Path(root)/'code/radius_enlargement_check.py').read_bytes()).hexdigest())
    with np.load(folder/'gap_extremes.npz', allow_pickle=False) as data:
        expected = {f'n{n}_{key}' for n in SIZES for key in ('largest','second')}
        expected |= {f'witness_eps{epsilon:g}_miss{outcome}' for epsilon in EPSILONS for outcome in (0,1)}
        require('radius_check_array_inventory', set(data.files)==expected)
        for n in SIZES:
            x,y = data[f'n{n}_largest'], data[f'n{n}_second']
            require(f'radius_check_gaps_{n}', x.shape==y.shape==(REPETITIONS,)
                and np.all((0<y)&(y<=x)&(x<1/3)))
            # Sorting is independent of the generator's two-element partition.
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([MASTER_SEED,n])))
            exp = rng.exponential(size=(BATCH_SIZE,n))
            total = exp.sum(axis=1)
            ordered = np.sort(exp,axis=1)
            require(f'radius_check_seed_replay_{n}',
                np.array_equal(x[:BATCH_SIZE],ordered[:,-1]/total)
                and np.array_equal(y[:BATCH_SIZE],ordered[:,-2]/total))
        rows = summarize(data)
        saved = list(csv.DictReader((folder/'summary.csv').open()))
        require('radius_check_csv_schema', len(saved)==12 and tuple(saved[0])==FIELDS)
        for row, record in zip(rows,saved):
            tag=f"{row['n']}_{row['epsilon']:g}"
            require('radius_check_recount_'+tag, all(float(record[k])==row[k] for k in FIELDS))
            require('radius_check_chord_subset_'+tag, row['misses']<=row['spacing_misses'])
            n, alpha = row['n'], 1+row['epsilon']
            def integrand(x):
                u=math.exp(-x/alpha)
                return 0. if u==1 else math.exp(-x+(n-1)*math.log1p(-u))
            integral=n*quad(integrand,0,math.inf,epsabs=1e-15,epsrel=1e-11)[0]
            require('radius_check_beta_integral_'+tag,
                abs(integral-row['spacing_probability'])<1e-9*row['spacing_probability'])
            if paper_text is not None:
                require('radius_check_display_'+tag, table_row(row).replace(',', '{,}') in paper_text)
        for epsilon in EPSILONS:
            for outcome in (0,1):
                gaps=data[f'witness_eps{epsilon:g}_miss{outcome}']
                require(f'radius_check_witness_gaps_{epsilon}_{outcome}',
                    len(gaps)==128 and abs(gaps.sum()-1)<1e-14 and np.all(gaps>0))
                s=np.r_[0.,np.cumsum(gaps[:-1])]
                xy=np.column_stack((np.cos(2*np.pi*s),np.sin(2*np.pi*s)))/(2*np.pi)
                chords=np.linalg.norm(xy[:,None,:]-xy[None,:,:],axis=2)
                y=np.sort(gaps)[-2]
                radius=(1+epsilon)*math.sin(math.pi*y)/math.pi
                edges=(chords<=radius)&~np.eye(128,dtype=bool)
                graph=csr_matrix(np.where(edges,chords,0.))
                require(f'radius_check_witness_connected_{epsilon}_{outcome}',
                    connected_components(graph,directed=False,return_labels=False)==1)
                distances=shortest_path(graph,directed=False)
                bars=ripser(distances,distance_matrix=True,maxdim=1)['dgms'][1]
                has_bar=any(d>b for b,d in bars)
                require(f'radius_check_witness_diagram_{epsilon}_{outcome}',has_bar==(not outcome))
    return dict(status='passed',new_clouds=400000,summary_rows=12,
        saved_gap_pairs=400000,independent_graph_diagrams=6,
        inference='Finite-size illustration; confidence intervals are pointwise, not a proof of an asymptotic law.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true',help='Read saved outputs; generate no Monte Carlo cohort.')
    args=parser.parse_args()
    if args.verify:
        def require(name, condition, detail=None):
            if not bool(condition): raise AssertionError(name)
        print(json.dumps(verify(ROOT,require),indent=2))
    else:
        print(json.dumps(generate(),indent=2))
