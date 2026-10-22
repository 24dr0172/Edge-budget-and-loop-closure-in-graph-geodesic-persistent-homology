"""The 24 cells of Table 3, rebuilt from the saved rank records only."""
from pathlib import Path
import collections
import csv
import math

PROFILES = ('uniform', 'strong_cosine', 'localized_trough')
METHODS = ('OR_kNN', 'mutual_kNN')
SIZES = (128, 512, 4096, 32768)
COLUMNS = ('profile', 'method', 'n', 'samples', 'closed_at_connectivity',
           'closure_percent', 'wilson_lo_percent', 'wilson_hi_percent')


def wilson(successes, total):
    z = 1.959963984540054
    p = successes / total
    den = 1 + z*z/total
    center = (p + z*z/(2*total)) / den
    half = z*math.sqrt(p*(1-p)/total + z*z/(4*total*total)) / den
    return center-half, center+half


def build_table(root):
    """Check the 24 groups against the raw counts and the saved summaries."""
    folder = Path(root)/'results/extra_k/rank_sweep'
    with (folder/'samples.csv').open(newline='') as stream:
        raw = list(csv.DictReader(stream))
    with (folder/'summary.csv').open(newline='') as stream:
        summaries = list(csv.DictReader(stream))
    groups = collections.defaultdict(list)
    for row in raw:
        if int(row['n']) in SIZES:
            groups[row['profile'], row['method'], int(row['n'])].append(row)
    expected = {(p, m, n) for p in PROFILES for m in METHODS for n in SIZES}
    if set(groups) != expected:
        raise ValueError('Missing or unexpected rank groups')
    summary = {}
    for row in summaries:
        if int(row['n']) in SIZES:
            key = row['profile'], row['method'], int(row['n'])
            if key in summary:
                raise ValueError('Duplicate summary group')
            summary[key] = row
    if set(summary) != expected:
        raise ValueError('Missing or unexpected summary groups')
    result = []
    paired_seeds = {}
    for profile in PROFILES:
        for method in METHODS:
            for n in SIZES:
                key = profile, method, n
                group = groups[key]
                seeds = {int(r['seed']) for r in group}
                count = 500 if n == 32768 else 1500
                if len(group) != count or len(seeds) != count:
                    raise ValueError(f'Wrong sample count or duplicate seeds: {key}')
                if n in paired_seeds and seeds != paired_seeds[n]:
                    raise ValueError(f'Broken profile/method seed pairing: {key}')
                paired_seeds[n] = seeds
                if not all(0 < float(r['largest_gap']) < 1/3 and
                           1 <= int(r['k_conn']) <= int(r['k_close']) < n
                           for r in group):
                    raise ValueError(f'Invalid rank or reference-gap condition: {key}')
                successes = sum(int(r['k_conn']) == int(r['k_close']) for r in group)
                lo, hi = wilson(successes, count)
                q = summary[key]
                values = {'cdf_0': successes/count, 'cdf_0_lo': lo, 'cdf_0_hi': hi}
                if int(q['samples']) != count or any(
                        abs(float(q[k])-v) > 1e-14 for k, v in values.items()):
                    raise ValueError(f'Raw-count/summary disagreement: {key}')
                result.append(dict(profile=profile, method=method, n=n, samples=count,
                    closed_at_connectivity=successes, closure_percent=100*successes/count,
                    wilson_lo_percent=100*lo, wilson_hi_percent=100*hi))
    return result


def write_table(root):
    rows = build_table(root)
    path = Path(root)/'results/rank_table.csv'
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def latex_rows(rows):
    names = dict(uniform='Uniform', strong_cosine='Cosine', localized_trough='Trough',
                 OR_kNN='OR-$k$NN', mutual_kNN='Mutual-$k$NN')
    return '\n'.join(' & '.join([names[p], names[m]] + [
        f"{next(r['closure_percent'] for r in rows if (r['profile'],r['method'],r['n']) == (p,m,n)):.1f}"
        for n in SIZES]) + r' \\' for p in PROFILES for m in METHODS)


def verify_table(root, require):
    expected = build_table(root)
    with (Path(root)/'results/rank_table.csv').open(newline='') as stream:
        reader = csv.DictReader(stream)
        require('retention_table_export_columns', reader.fieldnames == list(COLUMNS))
        actual = list(reader)
    require('retention_table_export_scope', len(actual) == len(expected) == 24)
    for i, (got, want) in enumerate(zip(actual, expected)):
        require(f'retention_table_export_cell_{i}', all(
            got[k] == str(v) if k in ('profile', 'method') else
            abs(float(got[k])-v) < 1e-11 for k, v in want.items()))
    return dict(cells=24, presentation_rows=6, size_columns=4,
                method_records=sum(r['samples'] for r in expected),
                source='results/extra_k/rank_sweep/samples.csv',
                reference_condition='Every largest_gap in the table is below 1/3.')
