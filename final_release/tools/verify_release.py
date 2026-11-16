"""Verify the code and saved data without requiring the manuscript files."""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def check(condition, description):
    if not condition:
        raise AssertionError(description)


def main():
    manifest = json.loads((ROOT / 'manifest.json').read_text())['sha256']
    files = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file()
             and p != ROOT / 'manifest.json' and '__pycache__' not in p.parts and p.suffix != '.pyc'}
    check(set(manifest) == files, 'release manifest inventory')
    for name, digest in manifest.items():
        check(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest,
              f'release checksum: {name}')

    sys.path.insert(0, str(ROOT / 'code'))
    import rank_table_sweep as sweep
    rank_dir = ROOT / 'results/detection_rank_sweep'
    summary = {(int(r['n']), r['profile'], r['method']): r for r in rows(rank_dir / 'summary.csv')}
    check(len(summary) == 24, 'rank summary inventory')
    for n, reps in sweep.SETTINGS:
        with np.load(rank_dir / f'samples_{n}.npz', allow_pickle=False) as data:
            ranks, gaps = data['ranks'], data['largest_gap']
        check(ranks.shape == (reps, 3, 2, 2) and gaps.shape == (reps, 3), f'rank arrays {n}')
        check(np.max(gaps) < 1/3 and np.all(ranks[..., 0] <= ranks[..., 1]), f'rank validity {n}')
        for j in (0, reps // 2, reps - 1):
            expected_ranks, expected_gaps = sweep.sample(n, j)
            check(np.array_equal(ranks[j], expected_ranks) and np.array_equal(gaps[j], expected_gaps),
                  f'rank seed replay {n}/{j}')
        for pi, profile in enumerate(sweep.PROFILES):
            for mi, method in enumerate(sweep.METHODS):
                r = summary[n, profile, method]
                extra = ranks[:, pi, mi, 1] - ranks[:, pi, mi, 0]
                count = int(np.count_nonzero(extra == 0))
                lo, hi = sweep.wilson(count, reps)
                check(int(r['samples']) == reps and int(r['tied']) == count and
                      abs(float(r['tie_percent']) - 100*count/reps) < 1e-11 and
                      abs(float(r['wilson_lo_percent']) - 100*lo) < 1e-11 and
                      abs(float(r['wilson_hi_percent']) - 100*hi) < 1e-11 and
                      all(int(r[f'plus{j}']) == np.count_nonzero(extra <= j) for j in (1, 2, 3)),
                      f'rank summary {n}/{profile}/{method}')
    audit = sweep.audit_local(29100000, 24)
    check(audit['samples'] == 24 and audit['full_connectivity_checks'] == 48, 'local/full rank checks')

    noise = rows(ROOT / 'results/radial_noise_knn.csv')
    check(len(noise) == 3200, 'radial-noise inventory')
    for r in noise:
        ratio = float(r['longest_ratio'])
        check(int(r['empty']) == int(ratio == 0) and
              int(r['at_least_quarter']) == int(ratio >= .25) and
              int(r['at_least_half']) == int(ratio >= .5), 'radial-noise flags')
    selected = {(r['seed'], r['method']): r for r in noise
                if float(r['delta']) == .01 and int(r['at_least_half']) == 1}
    winding = rows(ROOT / 'results/radial_noise_winding.csv')
    check(len(selected) == len(winding) == 258 and
          len({(r['seed'], r['method']) for r in winding}) == 258, 'winding inventory')
    for r in winding:
        original = selected[r['seed'], r['method']]
        check(r['k'] == original['k'] and original['all_consecutive'] == '0' and
              int(r['winding_one_fundamental_cycles']) > 0 and
              3 <= int(r['shortest_winding_one_cycle_vertices']) < 512, 'winding witness')

    first = rows(ROOT / 'results/first_nonempty/trials.csv')
    closure = rows(ROOT / 'results/peanut_neck/trials.csv')
    closed = {(r['alpha'], r['seed'], r['method']): r for r in closure}
    check(len(first) == len(closed) == 300, 'peanut inventory')
    for r in first:
        near = closed[r['alpha'], r['seed'], r['method']]
        conn, close, found = map(int, (r['M_conn'], r['M_close'], r['M_first']))
        check((conn, close) == (int(near['M_conn']), int(near['M_close'])) and
              conn <= found <= close and int(r['prefix_queries']) == found-conn+1 and
              (r['at_closure'] == 'True') == (found == close) and
              (r['early'] == 'True') == (found < close), 'peanut event flags')
        prefix = f"a{float(r['alpha']):g}_seed{r['seed']}"
        with np.load(ROOT / 'results/first_nonempty/checkpoints' / f"{prefix}_{r['method']}.npz", allow_pickle=False) as z:
            bars = z['first_bars']
            check(len(bars) == int(r['first_bars']) and len(bars) > 0 and
                  np.all(bars[:, 1] > bars[:, 0]) and
                  abs(float(np.max(bars[:, 1]-bars[:, 0]))-float(r['first_longest_lifetime'])) < 1e-14,
                  'peanut first diagram')
            if r['method'] == 'R2_global':
                check(len(z['double_bars']) == int(r['double_bars']), 'peanut doubled-radius diagram')
        with np.load(ROOT / 'results/peanut_neck/checkpoints' / f'{prefix}.npz', allow_pickle=False) as z:
            for event in ('connection', 'closure'):
                bars = z[f"{r['method']}_{event}_bars"]
                check(len(bars) == int(near[event+'_bars']) and np.all(bars[:, 1] > bars[:, 0]),
                      'peanut closure diagram')

    for result, script in (('exceptional_summary.json', 'check_exceptional.py'),
                           ('interpolation.json', 'check_interpolation.py')):
        meta = json.loads((ROOT / 'results/smooth_curve_checks' / result).read_text())
        check(meta['verification_script_sha256'] == hashlib.sha256((ROOT / 'code' / script).read_bytes()).hexdigest(),
              f'{result} script provenance')
    print('PASS: release checksums, rank arrays and seed replays, noise and winding records,')
    print('      peanut diagrams and event records, and numerical-check provenance.')
    print('This checks saved evidence; it does not check the manuscript or prove its theorems.')


if __name__ == '__main__':
    main()
