"""Compute noise and ellipse studies with resumable per-sample checkpoints."""
from pathlib import Path
import argparse, concurrent.futures, itertools, json, os, sys, time
sys.dont_write_bytecode = True
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
PAPER = Path(__file__).resolve().parents[1]
SUPPORT = PAPER
RELEASE = PAPER
sys.path.insert(0, str(RELEASE / 'code'))
import numpy as np
import pandas as pd
from ripser import ripser
import paper1_methods as m
from verify_arc_support import ring_death

DELTAS = (0., .0005, .001, .002, .005, .01)
METHODS = m.METHODS
LIFETIME_TOL = 1e-7


def one_bar_distance(bars, reference):
    """Exact matching options against one reference interval, in raw units."""
    assert reference.shape == (1, 2)
    base = float((reference[0, 1] - reference[0, 0]) / 2)
    if not len(bars):
        return base
    half = (bars[:, 1] - bars[:, 0]) / 2
    diagonal = max(base, float(half.max()))
    match = np.max(np.abs(bars - reference[0]), axis=1)
    if len(half) == 1:
        others = np.zeros(1)
    else:
        largest = int(np.argmax(half))
        others = np.full(len(half), half[largest])
        others[largest] = np.max(np.delete(half, largest))
    return min(diagonal, float(np.minimum(np.maximum(match, others), np.inf).min()))


def finite_graph_diagram(D, maximum_edge):
    """Exact H1 with an adaptive cutoff for a positive weighted graph metric.

    Above the largest generating graph edge, every Rips edge is homologous
    to its shortest graph path, through the fan based at its first endpoint.
    Thus H1 maps are surjective there: no later positive H1 interval can
    start. Once the cutoff is above that edge and has no surviving H1,
    the entire H1 barcode is known. A cone radius gives a safe last cutoff.
    """
    cone=float(np.min(np.max(D,axis=1)))
    first=float(np.nextafter(np.float32(maximum_edge),np.float32(np.inf)))
    cutoffs=sorted(set([first,max(first,.75*cone),max(first,float(np.nextafter(np.float32(cone),np.float32(np.inf))))]))
    for cutoff in cutoffs:
        bars=np.asarray(ripser(D,distance_matrix=True,maxdim=1,thresh=cutoff)['dgms'][1],float)
        if not len(bars) or np.isfinite(bars[:,1]).all():
            return bars[bars[:,1]>bars[:,0]]
    raise AssertionError('A universal vertex must kill all H1 at the cone radius')


def measurements(D, true, reference, maximum_edge=None):
    if not np.isfinite(D).all():
        return dict(connected=False, e_inf=np.inf, e_2=np.inf, e_PH=np.inf,
                    bars=0, longest=0., longest_over_reference=0.,
                    minimum_positive_lifetime=np.nan), np.empty((0, 2))
    bars = m.diagram(D) if maximum_edge is None else finite_graph_diagram(D,maximum_edge)
    diameter = float(true.max())
    errors = (D - true)[np.triu_indices(len(D), 1)]
    lifetimes = bars[:, 1] - bars[:, 0]
    assert not len(bars) or np.all(lifetimes > 0)
    reference_life = float(reference[0, 1] - reference[0, 0])
    eph = one_bar_distance(bars, reference) / diameter
    einf = float(np.max(np.abs(errors))) / diameter
    assert eph <= einf + 5e-8 / diameter
    return dict(connected=True, e_inf=einf,
                e_2=float(np.sqrt(np.mean(errors**2))) / diameter,
                e_PH=eph, bars=len(bars),
                longest=float(lifetimes.max()) if len(bars) else 0.,
                longest_over_reference=float(lifetimes.max()/reference_life) if len(bars) else 0.,
                minimum_positive_lifetime=float(lifetimes.min()) if len(bars) else np.nan), bars


def ring_tau(chords):
    """Track C by pairing it with the returned persistent F2 cocycles.

    Each representative is evaluated just before its interval's death.
    Intervals dying after C becomes a boundary pair trivially with C.
    The last nonzero pairing identifies C's death, rather than the longest bar.
    """
    n = len(chords)
    B = float(max(chords[i, (i+1) % n] for i in range(n)))
    result = ripser(chords, distance_matrix=True, maxdim=1, do_cocycles=True)
    ring = {tuple(sorted((i, (i+1) % n))) for i in range(n)}
    deaths = []
    for bar, cocycle in zip(result['dgms'][1], result['cocycles'][1]):
        if float(bar[1]) <= B:
            continue
        pairing = sum(int(value) for a, b, value in cocycle
                      if tuple(sorted((int(a), int(b)))) in ring) % 2
        if pairing:
            deaths.append(float(bar[1]))
    return max([B] + deaths)


def arc_support(D, chords, iu, ju, chosen, tol=1e-10):
    """The stated two-bound test, using cumulative maxima on both arcs."""
    n = len(D)
    ring = chords[np.arange(n), np.roll(np.arange(n), -1)]
    valid = np.zeros((n, n), bool)
    for a in range(n):
        for direction in (1, -1):
            arc = (a + direction*np.arange(n)) % n
            distances = np.maximum.accumulate(D[a, arc])
            steps = ring[arc[:-1]] if direction == 1 else ring[arc[1:]]
            costs = np.r_[0., np.maximum.accumulate(steps)]
            valid[a, arc] |= ((distances <= chords[a, arc] + tol) &
                              (costs <= chords[a, arc] + tol))
    supported = valid[iu[chosen], ju[chosen]] | valid[ju[chosen], iu[chosen]]
    return bool(supported.all()), int((~supported).sum())


def nn_crust(points):
    """Nearest-neighbor union, then half-neighbor completion at degree-one vertices.

    Degree tests use the first-stage graph, not edges added during completion.
    The input contains only Cartesian coordinates, with no cyclic labels.
    """
    chord = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    np.fill_diagonal(chord, np.inf)
    nearest = np.argmin(chord, axis=1)
    edges = {tuple(sorted((i, int(j)))) for i, j in enumerate(nearest)}
    neighbors = [set() for _ in points]
    for a, b in edges:
        neighbors[a].add(b); neighbors[b].add(a)
    extra = set()
    for i, adj in enumerate(neighbors):
        if len(adj) != 1:
            continue
        neighbor = next(iter(adj))
        dot = (points - points[i]) @ (points[neighbor] - points[i])
        admissible = np.flatnonzero(dot < 0.)
        if len(admissible):
            j = int(admissible[np.argmin(chord[i, admissible])])
            extra.add(tuple(sorted((i, j))))
    return edges | extra


def grid():
    from experiment_core import grid_values
    return [count for label,count in grid_values()]


def noise_job(seed, model, delta):
    stem = f'{model}_{seed}_delta{delta:g}'
    out = SUPPORT/'results/noise/checkpoints'
    out.mkdir(parents=True,exist_ok=True)
    csv, npz = out/f'{stem}.csv', out/f'{stem}.npz'
    if csv.exists() and npz.exists():
        return stem, 'resumed'
    s, clean = m.sample_cloud(seed, 'uniform')
    true = m.intrinsic(s); reference = m.diagram(true)
    if model == 'radial':
        normal = clean * (2*np.pi)
        draws = np.random.default_rng(seed+m.NOISE_SEED_OFFSET).uniform(-1., 1., len(s))
        displacement = draws[:, None]*normal
    else:
        rng = np.random.default_rng(seed+90000011)
        radius = np.sqrt(rng.uniform(size=len(s)))
        angle = rng.uniform(0., 2*np.pi, len(s))
        displacement = radius[:, None]*np.c_[np.cos(angle), np.sin(angle)]
    points = clean + delta*displacement
    chord, iu, ju, length, orders, keys = m.build_orders(points)
    budgets = grid(); rows = []
    arrays = dict(latent=s, clean=clean, observed=points,
                  displacement=displacement, reference=reference)
    for method in METHODS:
        order = orders[method]
        conn = m.first_connected(order, iu, ju, len(s))
        close = m.closure(order, iu, ju, len(s))
        for M in sorted(set(budgets + [conn, close-1, close])):
            D = m.graph_metric(order, iu, ju, length, M, len(s))
            scores, bars = measurements(D, true, reference,maximum_edge=float(length[order[:M]].max()))
            arrays[f'{method}_M{M}'] = bars
            rows.append(dict(seed=seed, model=model, delta=delta, method=method,
                             M=M, a_M=M/(len(s)*np.log(len(s))), n=len(s),
                             grid=M in budgets, first_connected=M == conn,
                             before_closure=M == close-1, closure=M == close,
                             M_conn=conn, M_close=close, **scores))
    pd.DataFrame(rows).to_csv(csv, index=False)
    np.savez_compressed(npz, **arrays)
    # The Euclidean baseline uses this independently generated observed cloud.
    import hashlib
    from scipy.spatial.distance import cdist
    baseline = SUPPORT/'results/phase1/euclidean/checkpoints'
    baseline.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    bars=m.diagram(cdist(points,points)); diameter=float(true.max())
    np.savez_compressed(baseline/npz.name,bars=bars,reference=reference,
        diameter=diameter,e_PH=one_bar_distance(bars,reference)/diameter,
        elapsed_seconds=time.monotonic()-started,
        source_sha256=hashlib.sha256(npz.read_bytes()).hexdigest())
    return stem, len(rows)


def ellipse_job(seed, profile):
    stem = f'{profile}_{seed}'
    out = SUPPORT/'results/ellipse/checkpoints'
    out.mkdir(parents=True,exist_ok=True)
    csv, npz, info = out/f'{stem}.csv', out/f'{stem}.npz', out/f'{stem}.json'
    if csv.exists() and npz.exists() and info.exists():
        return stem, 'resumed'
    with np.load(RELEASE/'results/ellipse_checks/coordinates.npz') as saved:
        s = saved[f'{seed}_{profile}_latent']; points = saved[f'{seed}_{profile}_xy']
    true = m.intrinsic(s); reference = m.diagram(true)
    chord, iu, ju, length, orders, keys = m.build_orders(points)
    tau = ring_tau(chord)
    B = float(chord[np.arange(len(s)), np.roll(np.arange(len(s)), -1)].max())
    g = float(np.diff(np.r_[s, s[0]+1]).max())
    bstar = float(reference[0, 1])
    permutation = np.random.default_rng(seed+110000009).permutation(len(s))
    recovered = {tuple(sorted((int(permutation[a]), int(permutation[b]))))
                 for a, b in nn_crust(points[permutation])}
    expected = {tuple(sorted((a, (a+1) % len(s)))) for a in range(len(s))}
    summary = dict(seed=seed, profile=profile, tau_E=tau, B=B, g=g,
                   b_star=bstar, margin_tau=tau-(g+bstar)/2,
                   margin_birth=B-(3*g-bstar)/2,
                   condition_5c=bool(tau > (g+bstar)/2+1e-6 and B > (3*g-bstar)/2+1e-6),
                   nn_crust_exact_ring=recovered == expected,
                   nn_crust_extra_edges=len(recovered-expected),
                   nn_crust_missing_edges=len(expected-recovered))
    arrays = dict(latent=s, observed=points, reference=reference,
                  nn_crust_edges=np.array(sorted(recovered), int))
    budgets = grid(); rows = []
    for method in METHODS:
        order = orders[method]; close = m.closure(order, iu, ju, len(s))
        nonlocal_positions = np.flatnonzero(length[order] >= .18)
        first_nonlocal = int(nonlocal_positions[0]+1) if len(nonlocal_positions) else 0
        selected = sorted(set([M for M in budgets if M >= close] + [close, max(close, first_nonlocal)]))
        for M in selected:
            D = m.graph_metric(order, iu, ju, length, M, len(s))
            scores, bars = measurements(D, true, reference,maximum_edge=float(length[order[:M]].max()))
            ok, failures = arc_support(D, chord, iu, ju, order[:M])
            arrays[f'{method}_M{M}'] = bars
            rows.append(dict(seed=seed, profile=profile, method=method, M=M,
                             a_M=M/(len(s)*np.log(len(s))), grid=M in budgets,
                             closure=M == close, first_nonlocal=M == first_nonlocal,
                             M_close=close, M_first_nonlocal=first_nonlocal,
                             maximum_chord=float(length[order[:M]].max()),
                             arc_supported=ok, failing_edges=failures,
                             condition_5c=summary['condition_5c'], **scores))
    pd.DataFrame(rows).to_csv(csv, index=False)
    np.savez_compressed(npz, **arrays)
    info.write_text(json.dumps(summary, indent=2)+'\n')
    return stem, len(rows)


def validate_helpers():
    from persim import bottleneck
    from verify_arc_support import supported
    rng = np.random.default_rng(20261015)
    for _ in range(100):
        bars = np.sort(rng.uniform(size=(int(rng.integers(0, 10)), 2)), axis=1)
        reference = np.sort(rng.uniform(size=(1, 2)), axis=1)
        assert abs(one_bar_distance(bars, reference)-bottleneck(bars, reference)) < 1e-12
    paired = 0
    for _ in range(80):
        s = np.sort(rng.uniform(size=8))
        xy = np.c_[np.cos(2*np.pi*s), np.sin(2*np.pi*s)]
        xy[:, 0] *= rng.uniform(.5, 3.)
        chords, iu, ju, length, orders, keys = m.build_orders(xy)
        exact_tau = ring_death(chords, max(chords[i, (i+1)%8] for i in range(8)))
        assert abs(ring_tau(chords)-exact_tau) < 1e-6
        for M in [m.closure(orders['R2_global'], iu, ju, 8), 20, 28]:
            order = orders['R2_global']; D = m.graph_metric(order, iu, ju, length, M, 8)
            fast, _ = arc_support(D, chords, iu, ju, order[:M])
            edges = [(int(iu[e]), int(ju[e])) for e in order[:M]]
            assert fast == supported(D, chords, edges)
            paired += 1
    report = dict(status='passed', one_bar_matching_cases=100,
                  ring_cocycle_vs_boundary_span_cases=80,
                  cumulative_vs_direct_support_cases=paired)
    (SUPPORT/'results/helper_validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(report, flush=True)


def aggregate():
    noise = pd.concat([pd.read_csv(p) for p in sorted((SUPPORT/'results/noise/checkpoints').glob('*.csv'))], ignore_index=True)
    clean = noise[(noise.model == 'radial') & (noise.delta == 0)].copy(); clean['model'] = 'isotropic'
    noise = pd.concat([noise, clean], ignore_index=True)
    noise.to_csv(SUPPORT/'results/noise/cells.csv', index=False)
    selectors = []
    for key, group in noise[noise.grid].groupby(['seed','model','delta','method']):
        group = group.sort_values('M'); valid = group[group.connected & (group.bars > 0)]
        for selector in ['first_nonempty', 'longest_lifetime']:
            if len(valid):
                chosen = valid.iloc[0] if selector == 'first_nonempty' else valid[valid.longest >= valid.longest.max()-LIFETIME_TOL].iloc[0]
                selectors.append(dict(zip(['seed','model','delta','method'], key)) | dict(selector=selector, selected=True, M=int(chosen.M), M_close=int(chosen.M_close), before_closure=bool(chosen.M < chosen.M_close), e_PH=float(chosen.e_PH), longest_over_reference=float(chosen.longest_over_reference)))
            else:
                selectors.append(dict(zip(['seed','model','delta','method'], key)) | dict(selector=selector, selected=False, M=np.nan, M_close=int(group.M_close.iloc[0]), before_closure=False, e_PH=np.inf, longest_over_reference=0.))
    pd.DataFrame(selectors).to_csv(SUPPORT/'results/noise/selectors.csv', index=False)
    ellipse = pd.concat([pd.read_csv(p) for p in sorted((SUPPORT/'results/ellipse/checkpoints').glob('*.csv'))], ignore_index=True)
    ellipse.to_csv(SUPPORT/'results/ellipse/cells.csv', index=False)
    samples = [json.loads(p.read_text()) for p in sorted((SUPPORT/'results/ellipse/checkpoints').glob('*.json'))]
    pd.DataFrame(samples).to_csv(SUPPORT/'results/ellipse/samples.csv', index=False)
    print('Aggregated',len(noise),'noise rows and',len(ellipse),'ellipse rows.',flush=True)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--workers',type=int,default=6)
    parser.add_argument('--validate-only',action='store_true');parser.add_argument('--smoke',action='store_true')
    parser.add_argument('--aggregate-only',action='store_true'); args = parser.parse_args()
    if args.aggregate_only:aggregate();return
    validate_helpers()
    if args.validate_only:return
    if args.smoke:
        start=time.monotonic();print(noise_job(20320101,'radial',.01),flush=True)
        print('Noise block seconds',time.monotonic()-start,flush=True)
        start=time.monotonic();print(ellipse_job(20310101,'uniform'),flush=True)
        print('Ellipse sample seconds',time.monotonic()-start,flush=True);return
    jobs = []
    # Interleave the two independent studies without changing either release input.
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for seed in range(20310101,20310121):
            for profile in m.PROFILES:jobs.append(pool.submit(ellipse_job,seed,profile))
        for seed in range(20320101,20320121):
            for model in ['radial','isotropic']:
                for delta in DELTAS:
                    if model=='isotropic' and delta==0:continue
                    jobs.append(pool.submit(noise_job,seed,model,delta))
        for i,future in enumerate(concurrent.futures.as_completed(jobs),1):
            print(i,'/',len(jobs),future.result(),flush=True)
    aggregate()


if __name__=='__main__':main()
