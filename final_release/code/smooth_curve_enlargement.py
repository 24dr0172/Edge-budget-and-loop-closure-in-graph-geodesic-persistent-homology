"""Generate or verify the radius-enlargement miss rates on two smooth curves.

Two fixed length-one curves are used: the 2:1 ellipse of the ellipse study and
the a=0.4 peanut of the varying-waist study, with the same arc-length tables.
Uniform arc-length samples come from normalized exponential spacings and a
uniform random starting point. Each curve has its own independent stream of
draws: batch b at size n of curve c is seeded by SeedSequence([MASTER_SEED,
CURVE_STREAMS[c], n, b]), so the ellipse and peanut rows are independent
estimates. A miss at factor alpha is the consecutive-chord event
B > alpha * r2, where B and r2 are the largest and second-largest chords
between consecutive samples. Within each cloud the three factors are paired.

Only the few largest arc gaps are mapped to the curve. Every chord is at most
its arc gap, so the two largest chords are certified to lie among the mapped
gaps when the second-largest mapped chord is at least the largest unmapped gap;
otherwise every gap of that cloud is mapped. The theorem's good event
g < min(r_loc, lambda/3) is recorded with the grid estimates of r_loc
(Appendix A, two-curve check) and of the curve's chord-arc ratio lambda.

Full persistence spot checks at n=128 and n=512 compute the first connecting
Euclidean radius from a minimum spanning tree, build the radius graph at
alpha times it, take chord-weighted shortest paths and compute the complete
H1 diagram with unrestricted Ripser.
"""
from pathlib import Path
from functools import lru_cache
import argparse, csv, hashlib, json, math, os, sys
sys.dont_write_bytecode = True
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FOLDER = 'results/smooth_curve_enlargement'
CURVES = ('ellipse', 'peanut')
CURVE_STREAMS = {'ellipse': 1, 'peanut': 2}
PEANUT_A = .4
SIZES = (128, 512, 4096, 32768)
ALPHAS = (1.25, 1.5, 2.)
REPETITIONS = 1000000
BATCH = 1000
MASTER_SEED = 20261027
CANDIDATES = 8
CHUNK_VALUES = 2**22
ARC_GRID = 400001
LAMBDA_GRIDS = (4000, 8000)
SPOT_ALL = {128: 2000, 512: 150}
SPOT_EXTRA_MISSES = {128: 0, 512: 25}
SPOT_EXTRA_SCAN = {128: 0, 512: 20000}
BATCH_FIELDS = ['curve', 'n', 'batch', 'samples'] + \
    [f'miss_{a:g}' for a in ALPHAS] + [f'spacing_miss_{a:g}' for a in ALPHAS] + \
    ['good_event', 'fallback_clouds', 'largest_gap', 'largest_chord_over_gap_min']
SUMMARY_FIELDS = ('curve', 'n', 'alpha', 'samples', 'misses', 'miss_rate', 'wilson_lo', 'wilson_hi',
                  'spacing_misses', 'spacing_probability', 'asymptotic_probability',
                  'observed_over_asymptotic', 'good_event', 'good_event_rate')
SPOT_FIELDS = ('curve', 'n', 'cloud', 'alpha', 'largest_gap', 'good_event', 'B', 'r_conn',
               'r_conn_is_second_chord', 'chord_miss', 'connected', 'h1_bars', 'nonempty',
               'agrees')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def arc_table(curve):
    """The arc-length tables of experiment_core.ellipse_points and peanut_neck_study.arc_table."""
    if curve == 'ellipse':
        theta = np.linspace(0, 2*np.pi, ARC_GRID)
        speed = np.sqrt(4*np.sin(theta)**2+np.cos(theta)**2)
        arc = np.r_[0, np.cumsum((speed[1:]+speed[:-1])*.5*np.diff(theta))]
        perimeter = float(arc[-1])
        return theta, arc/perimeter, perimeter
    theta = np.linspace(0., 2*np.pi, ARC_GRID)
    radius = 1+PEANUT_A*np.cos(2*theta)
    xy = np.c_[radius*np.cos(theta), radius*np.sin(theta)]
    cumulative = np.r_[0., np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
    length = float(cumulative[-1])
    return theta, cumulative/length, length


def points(curve, s):
    """Length-one curve points at arc positions s (taken modulo one)."""
    theta, arc, length = arc_table(curve)
    t = np.interp(np.asarray(s) % 1, arc, theta)
    if curve == 'ellipse':
        return np.stack((2*np.cos(t), np.sin(t)), axis=-1)/length
    radius = 1+PEANUT_A*np.cos(2*t)
    return np.stack((radius*np.cos(t), radius*np.sin(t)), axis=-1)/length


def chord_arc_ratio(curve, grid):
    """Grid estimate of lambda = min chord / intrinsic distance over distinct curve points."""
    s = np.arange(grid)/grid
    xy = points(curve, s)
    best = np.inf
    for start in range(0, grid, 500):
        d = np.linalg.norm(xy[start:start+500, None, :]-xy[None, :, :], axis=2)
        sep = np.abs(s[start:start+500, None]-s[None, :])
        arc = np.minimum(sep, 1-sep)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(arc > 0, d/arc, np.inf)
        best = min(best, float(ratio.min()))
    return best


def constants(root=ROOT):
    locality = json.loads((Path(root)/'results/noncircular_locality/locality_check_summary.json').read_text())
    out = {}
    for curve in CURVES:
        lam = {str(g): chord_arc_ratio(curve, g) for g in LAMBDA_GRIDS}
        r_loc = float(locality['r_loc'][curve]['r_loc'])
        lam_value = min(lam.values())
        out[curve] = dict(r_loc=r_loc, r_loc_source='results/noncircular_locality/locality_check_summary.json',
                          chord_arc_ratio_grid=lam, chord_arc_ratio=lam_value,
                          good_event_threshold=min(r_loc, lam_value/3),
                          perimeter_before_normalization=arc_table(curve)[2])
    return out


def batch_rng(curve, n, batch):
    """Independent PCG64 stream for one curve, size and batch."""
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence([MASTER_SEED, CURVE_STREAMS[curve], n, batch])))


def draw(rng, count, n):
    """Starting points and spacings of `count` uniform arc-length clouds."""
    start = rng.random(count)
    spacing = rng.exponential(size=(count, n))
    spacing /= spacing.sum(axis=1)[:, None]
    return start, spacing


def chords_of(curve, start, spacing, gaps_index):
    """Chords of the selected consecutive gaps; positions use the cumulative spacings."""
    rows = np.arange(len(start))[:, None]
    left = np.cumsum(spacing, axis=1)-spacing
    a = start[:, None]+left[rows, gaps_index]
    b = a+spacing[rows, gaps_index]
    return np.linalg.norm(points(curve, b)-points(curve, a), axis=-1)


def top_two(curve, start, spacing):
    """Largest and second-largest consecutive chords, certified by the chord <= arc bound."""
    count, n = spacing.shape
    k = min(CANDIDATES, n-1)
    part = np.argpartition(spacing, n-k-1, axis=1)
    candidates, rest_max = part[:, n-k:], np.take_along_axis(spacing, part[:, n-k-1:n-k], axis=1)[:, 0]
    chords = chords_of(curve, start, spacing, candidates)
    ordered = np.sort(chords, axis=1)
    B, r2 = ordered[:, -1].copy(), ordered[:, -2].copy()
    # A chord never exceeds its arc gap; 1e-9 covers the arc-table rounding.
    fallback = ~(r2 >= rest_max*(1+1e-9))
    for i in np.flatnonzero(fallback):
        all_chords = chords_of(curve, start[i:i+1], spacing[i:i+1], np.arange(n)[None, :])[0]
        top = np.sort(all_chords)
        B[i], r2[i] = top[-1], top[-2]
    gaps = np.take_along_axis(spacing, candidates, axis=1)
    largest = gaps.max(axis=1)
    ratio = np.take_along_axis(chords, gaps.argmax(axis=1)[:, None], axis=1)[:, 0]/largest
    second_gap = np.sort(gaps, axis=1)[:, -2]
    return B, r2, largest, second_gap, ratio, fallback


def run_batch(job):
    n, batch, thresholds = job
    rows = {c: dict(curve=c, n=n, batch=batch, samples=0, good_event=0, fallback_clouds=0,
                    largest_gap=0., largest_chord_over_gap_min=1.,
                    **{f'miss_{a:g}': 0 for a in ALPHAS}, **{f'spacing_miss_{a:g}': 0 for a in ALPHAS})
            for c in CURVES}
    chunk = max(1, CHUNK_VALUES//n)
    for curve in CURVES:
        rng = batch_rng(curve, n, batch)
        done = 0
        while done < BATCH:
            count = min(chunk, BATCH-done)
            start, spacing = draw(rng, count, n)
            B, r2, g, g2, ratio, fallback = top_two(curve, start, spacing)
            row = rows[curve]
            row['samples'] += count
            for a in ALPHAS:
                row[f'miss_{a:g}'] += int((B > a*r2).sum())
                row[f'spacing_miss_{a:g}'] += int((g > a*g2).sum())
            row['good_event'] += int((g < thresholds[curve]).sum())
            row['fallback_clouds'] += int(fallback.sum())
            row['largest_gap'] = max(row['largest_gap'], float(g.max()))
            row['largest_chord_over_gap_min'] = min(row['largest_chord_over_gap_min'], float(ratio.min()))
            done += count
    return [rows[c] for c in CURVES]


def wilson(k, n):
    z = 1.959963984540054
    p = k/n
    den = 1+z*z/n
    centre = (p+z*z/(2*n))/den
    half = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return centre-half, centre+half


def predictions(n, alpha):
    spacing = math.exp(math.lgamma(alpha+1)+math.lgamma(n+1)-math.lgamma(n+alpha))
    return spacing, math.gamma(alpha+1)*n**(1-alpha)


def summarize(batches):
    rows = []
    for curve in CURVES:
        for alpha in ALPHAS:
            for n in SIZES:
                group = [b for b in batches if b['curve'] == curve and int(b['n']) == n]
                samples = sum(int(b['samples']) for b in group)
                k = sum(int(b[f'miss_{alpha:g}']) for b in group)
                good = sum(int(b['good_event']) for b in group)
                lo, hi = wilson(k, samples)
                spacing, asymptotic = predictions(n, alpha)
                rows.append(dict(zip(SUMMARY_FIELDS, (curve, n, alpha, samples, k, k/samples, lo, hi,
                    sum(int(b[f'spacing_miss_{alpha:g}']) for b in group), spacing, asymptotic,
                    (k/samples)/asymptotic, good, good/samples))))
    return rows


def spot_clouds(n, thresholds):
    """Cloud indices for the persistence spot check: the first clouds, then early alpha=2 misses."""
    chosen = list(range(SPOT_ALL[n]))
    if SPOT_EXTRA_MISSES[n]:
        found = {c: [] for c in CURVES}
        for batch in range(SPOT_EXTRA_SCAN[n]//BATCH):
            for curve in CURVES:
                start, spacing = draw(batch_rng(curve, n, batch), BATCH, n)
                B, r2, *_ = top_two(curve, start, spacing)
                for i in np.flatnonzero(B > ALPHAS[-1]*r2):
                    index = batch*BATCH+int(i)
                    if index >= SPOT_ALL[n] and len(found[curve]) < SPOT_EXTRA_MISSES[n]:
                        found[curve].append(index)
        return {c: chosen+found[c] for c in CURVES}
    return {c: chosen for c in CURVES}


def cloud(curve, n, index):
    """Replay one cloud of one curve: starting point and spacings."""
    rng = batch_rng(curve, n, index//BATCH)
    chunk = max(1, CHUNK_VALUES//n)
    offset = index % BATCH
    done = 0
    while True:
        count = min(chunk, BATCH-done)
        start, spacing = draw(rng, count, n)
        if offset < done+count:
            i = offset-done
            return start[i], spacing[i]
        done += count


def spot_job(job):
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import shortest_path, minimum_spanning_tree, connected_components
    from ripser import ripser
    curve, n, index, threshold = job
    start, spacing = cloud(curve, n, index)
    s = (start+np.r_[0., np.cumsum(spacing)[:-1]]) % 1
    xy = points(curve, s)
    C = np.linalg.norm(xy[:, None, :]-xy[None, :, :], axis=2)
    consecutive = np.sort(C[np.arange(n), (np.arange(n)+1) % n])
    B, second = float(consecutive[-1]), float(consecutive[-2])
    r_conn = float(minimum_spanning_tree(C).data.max())
    out = []
    for alpha in ALPHAS:
        A = np.where(C <= alpha*r_conn, C, 0.)
        np.fill_diagonal(A, 0.)
        graph = csr_matrix(A)
        connected = connected_components(graph, directed=False, return_labels=False) == 1
        D = shortest_path(graph, directed=False, method='D')
        D = (D+D.T)/2
        np.fill_diagonal(D, 0.)
        bars = np.asarray(ripser(D, distance_matrix=True, maxdim=1)['dgms'][1], float).reshape(-1, 2)
        bars = bars[bars[:, 1] > bars[:, 0]]
        miss = B > alpha*r_conn
        out.append(dict(zip(SPOT_FIELDS, (curve, n, index, alpha, float(spacing.max()),
            bool(spacing.max() < threshold), B, r_conn, bool(r_conn == second), bool(miss), bool(connected),
            len(bars), bool(len(bars) > 0), bool((len(bars) > 0) == (not miss))))))
    return out


def write_csv(path, fields, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fields))
        writer.writeheader(); writer.writerows(rows)


def generate(root=ROOT, workers=4):
    from multiprocessing import Pool
    folder = Path(root)/FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    geometry = constants(root)
    (folder/'curve_constants.json').write_text(json.dumps(geometry, indent=2)+'\n')
    thresholds = {c: geometry[c]['good_event_threshold'] for c in CURVES}
    jobs = [(n, b, thresholds) for n in SIZES for b in range(REPETITIONS//BATCH)]
    batches = []
    with Pool(workers) as pool:
        for i, rows in enumerate(pool.imap(run_batch, jobs, chunksize=4), 1):
            batches.extend(rows)
            if i % 500 == 0:
                print(f'Smooth-curve enlargement: {i}/{len(jobs)} batches', flush=True)
    batches.sort(key=lambda r: (CURVES.index(r['curve']), r['n'], r['batch']))
    write_csv(folder/'batches.csv', BATCH_FIELDS, batches)
    summary = summarize(batches)
    write_csv(folder/'summary.csv', SUMMARY_FIELDS, summary)
    spot_jobs = []
    for n in SPOT_ALL:
        for curve, indices in spot_clouds(n, thresholds).items():
            spot_jobs += [(curve, n, i, thresholds[curve]) for i in indices]
    spots = []
    with Pool(workers) as pool:
        for rows in pool.imap(spot_job, spot_jobs, chunksize=2):
            spots.extend(rows)
    spots.sort(key=lambda r: (CURVES.index(r['curve']), r['n'], r['cloud'], r['alpha']))
    write_csv(folder/'spot_checks.csv', SPOT_FIELDS, spots)
    config = dict(curves={'ellipse': '(2 cos t, sin t), normalized to length one; trapezoid arc table as in experiment_core.ellipse_points',
                          'peanut': 'r = 1 + 0.4 cos 2t, normalized to length one; polyline arc table as in peanut_neck_study.arc_table'},
                  arc_grid_points=ARC_GRID, sizes=list(SIZES), alphas=list(ALPHAS),
                  samples_per_size=REPETITIONS, batch_size=BATCH, master_seed=MASTER_SEED,
                  curve_streams=CURVE_STREAMS,
                  seed_sequence='[master_seed, curve_stream, n, batch]', generator='PCG64',
                  sampling='uniform starting point, then n independent Exp(1) spacings divided by their sum',
                  pairing='Each curve has its own independent draws (its own curve_stream); within a curve the three factors reuse each cloud.',
                  miss_event='B > alpha * r2, the two largest consecutive chords',
                  spacing_event='g_largest > alpha * g_second, the arc-gap event with exact probability Gamma(alpha+1)Gamma(n+1)/Gamma(n+alpha)',
                  asymptotic='Gamma(alpha+1) * n^(1-alpha)',
                  good_event='largest arc gap < min(r_loc, lambda/3), with grid estimates of r_loc and lambda',
                  mapped_gaps=CANDIDATES,
                  mapping_rule='The two largest chords are certified to lie among the mapped gaps when the second-largest mapped chord is at least (1+1e-9) times the largest unmapped gap; otherwise all gaps are mapped.',
                  interval='Pointwise two-sided 95% Wilson',
                  spot_checks=dict(first_clouds={str(k): v for k, v in SPOT_ALL.items()},
                                   added_alpha2_misses={str(k): v for k, v in SPOT_EXTRA_MISSES.items()},
                                   miss_search_clouds={str(k): v for k, v in SPOT_EXTRA_SCAN.items()},
                                   method='r_conn from the Euclidean minimum spanning tree; closed radius graph at alpha*r_conn; chord-weighted Dijkstra; unrestricted Ripser over F2; every bar with death > birth kept'),
                  new_clouds=len(SIZES)*REPETITIONS, discarded_clouds=0,
                  script_sha256=sha(__file__))
    (folder/'config.json').write_text(json.dumps(config, indent=2)+'\n')
    return summary


def table_rows(summary):
    """LaTeX rows for the smooth-curve miss-rate table (percentages)."""
    lines = []
    for alpha in ALPHAS:
        for n in SIZES:
            cells = [f'${alpha:g}$', f'${n:,}$']
            for curve in CURVES:
                r = next(v for v in summary if v['curve'] == curve and v['n'] == n and v['alpha'] == alpha)
                cells.append(f"${100*r['miss_rate']:.4f}$")
                cells.append(f"$[{100*r['wilson_lo']:.4f}, {100*r['wilson_hi']:.4f}]$")
            spacing, asymptotic = predictions(n, alpha)
            cells += [f'${100*spacing:.4f}$', f'${100*asymptotic:.4f}$']
            lines.append((' & '.join(cells)+r'\\').replace(',', '{,}'))
    return lines


def read_summary(folder):
    rows = list(csv.DictReader((folder/'summary.csv').open()))
    for r in rows:
        for k in r:
            if k in ('n', 'samples', 'misses', 'spacing_misses', 'good_event'):
                r[k] = int(r[k])
            elif k != 'curve':
                r[k] = float(r[k])
    return rows


def verify(root, require, replay_batches=2):
    """Recount saved batches, replay seeded batches and re-read the persistence spot checks."""
    root = Path(root)
    folder = root/FOLDER
    config = json.loads((folder/'config.json').read_text())
    require('smooth_enlargement_protocol', config['sizes'] == list(SIZES) and config['alphas'] == list(ALPHAS)
            and config['samples_per_size'] == REPETITIONS and config['master_seed'] == MASTER_SEED
            and config.get('curve_streams') == CURVE_STREAMS
            and config['new_clouds'] == len(SIZES)*REPETITIONS and config['discarded_clouds'] == 0)
    require('smooth_enlargement_script_hash', config['script_sha256'] == sha(root/'code/smooth_curve_enlargement.py'))
    geometry = json.loads((folder/'curve_constants.json').read_text())
    locality = json.loads((root/'results/noncircular_locality/locality_check_summary.json').read_text())
    for curve in CURVES:
        g = geometry[curve]
        require('smooth_enlargement_constants_'+curve, g['r_loc'] == locality['r_loc'][curve]['r_loc'] and
                abs(g['good_event_threshold']-min(g['r_loc'], g['chord_arc_ratio']/3)) < 1e-15 and
                max(g['chord_arc_ratio_grid'].values())-min(g['chord_arc_ratio_grid'].values()) < 1e-9)
    # The cached tables reproduce the study samplers.
    from experiment_core import ellipse_points
    from peanut_neck_study import sample as peanut_sample
    s = np.sort(np.random.default_rng(7).random(64))
    require('smooth_enlargement_ellipse_sampler', np.array_equal(points('ellipse', s), ellipse_points(s)))
    _, xy, _ = peanut_sample(20361021, PEANUT_A)
    s0 = np.sort(np.random.default_rng(20361021).random(512))
    require('smooth_enlargement_peanut_sampler', np.max(np.abs(points('peanut', s0)-xy)) < 1e-15)
    batches = list(csv.DictReader((folder/'batches.csv').open()))
    require('smooth_enlargement_batch_schema', tuple(batches[0]) == tuple(BATCH_FIELDS) and
            len(batches) == len(CURVES)*len(SIZES)*(REPETITIONS//BATCH))
    for curve in CURVES:
        for n in SIZES:
            got = sorted(int(b['batch']) for b in batches if b['curve'] == curve and int(b['n']) == n)
            require(f'smooth_enlargement_batches_{curve}_{n}', got == list(range(REPETITIONS//BATCH)) and
                    all(int(b['samples']) == BATCH for b in batches if b['curve'] == curve and int(b['n']) == n))
    saved = read_summary(folder)
    rows = summarize(batches)
    require('smooth_enlargement_summary_rows', len(saved) == len(rows) == len(CURVES)*len(SIZES)*len(ALPHAS))
    for row, record in zip(rows, saved):
        tag = f"{row['curve']}_{row['n']}_{row['alpha']:g}"
        require('smooth_enlargement_recount_'+tag, all(record[k] == row[k] for k in SUMMARY_FIELDS))
        require('smooth_enlargement_good_event_'+tag, row['good_event'] <= row['samples'])
    thresholds = {c: geometry[c]['good_event_threshold'] for c in CURVES}
    for n in SIZES:
        for batch in sorted({0, REPETITIONS//BATCH-1})[:replay_batches]:
            replay = run_batch((n, batch, thresholds))
            for row in replay:
                record = next(b for b in batches if b['curve'] == row['curve'] and int(b['n']) == n and int(b['batch']) == batch)
                require(f"smooth_enlargement_replay_{row['curve']}_{n}_{batch}",
                        all(str(row[k]) == record[k] if isinstance(row[k], (int, str)) else float(record[k]) == row[k]
                            for k in BATCH_FIELDS))
    spots = list(csv.DictReader((folder/'spot_checks.csv').open()))
    require('smooth_enlargement_spot_schema', tuple(spots[0]) == SPOT_FIELDS)
    counts = {}
    for curve in CURVES:
        for n in SPOT_ALL:
            q = [r for r in spots if r['curve'] == curve and int(r['n']) == n]
            clouds = sorted({int(r['cloud']) for r in q})
            require(f'smooth_enlargement_spot_scope_{curve}_{n}', clouds[:SPOT_ALL[n]] == list(range(SPOT_ALL[n])) and
                    len(clouds) == SPOT_ALL[n]+SPOT_EXTRA_MISSES[n] and len(q) == len(clouds)*len(ALPHAS))
            require(f'smooth_enlargement_spot_agreement_{curve}_{n}', all(r['agrees'] == 'True' and r['connected'] == 'True' for r in q))
            counts[f'{curve}_{n}'] = dict(diagrams=len(q), misses=sum(r['chord_miss'] == 'True' for r in q),
                                         r_conn_second_chord=sum(r['r_conn_is_second_chord'] == 'True' for r in q)//len(ALPHAS),
                                         good_event_clouds=sum(r['good_event'] == 'True' for r in q)//len(ALPHAS))
    # Replay one spot cloud per curve and size, including its persistence diagrams.
    for curve in CURVES:
        for n in SPOT_ALL:
            index = max(int(r['cloud']) for r in spots if r['curve'] == curve and int(r['n']) == n)
            again = spot_job((curve, n, index, thresholds[curve]))
            for row in again:
                record = next(r for r in spots if r['curve'] == curve and int(r['n']) == n and
                              int(r['cloud']) == index and float(r['alpha']) == row['alpha'])
                require(f"smooth_enlargement_spot_replay_{curve}_{n}_{index}_{row['alpha']:g}",
                        all(str(row[k]) == record[k] for k in ('chord_miss', 'nonempty', 'h1_bars', 'agrees')))
    return dict(status='passed', clouds=len(SIZES)*REPETITIONS, curves=list(CURVES), summary_rows=len(rows),
                replayed_batches=len(SIZES)*min(replay_batches, 2), spot_checks=counts,
                inference='Finite-size check of the smooth-curve miss law; intervals are pointwise.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', help='Read saved outputs and replay seeded batches; generate no cohort.')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.verify:
        def require(name, condition, detail=None):
            if not bool(condition): raise AssertionError(name)
        print(json.dumps(verify(ROOT, require), indent=2))
    else:
        rows = generate(ROOT, args.workers)
        for line in table_rows(rows): print(line)


if __name__ == '__main__':
    main()
