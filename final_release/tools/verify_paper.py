"""Check the detection-focused manuscript against its saved evidence.

This checks numerical claims, references, figure inventory, selected seed
replays, and the compiled PDF's basic identity. It is not a proof checker.
"""
from pathlib import Path
import argparse
import csv
import importlib.util
import json
import math
import re
from statistics import median
import subprocess
import sys
import hashlib
import zipfile
sys.dont_write_bytecode = True
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results' / 'detection_rank_sweep'


def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def section(text, start, end):
    start_at = text.index(start)
    return text[start_at:text.index(end, start_at)]


def check(condition, description):
    if not condition:
        raise AssertionError(description)


def main(paper):
    text = (paper / 'combined.tex').read_text()
    labels = re.findall(r'\\label\{([^}]+)\}', text)
    refs = re.findall(r'\\(?:ref|eqref)\{([^}]+)\}', text)
    citations = {c.strip() for block in re.findall(
        r'\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}', text) for c in block.split(',')}
    bibliography = re.findall(r'\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}', text)
    check(len(labels) == len(set(labels)), 'duplicate labels')
    check(not (set(refs) - set(labels)), 'unresolved cross-references')
    check(len(bibliography) == len(set(bibliography)), 'duplicate bibliography keys')
    check(citations == set(bibliography), 'unresolved or unused bibliography entries')
    citation_number = {key: i for i, key in enumerate(bibliography, start=1)}
    for block in re.findall(r'\\citep\*?(?:\[[^\]]*\])*\{([^}]+)\}', text):
        numbers = [citation_number[key.strip()] for key in block.split(',')]
        check(numbers == sorted(numbers), f'numeric citation order: {block}')
    check(all((ROOT / x).exists() for x in re.findall(r'\\path\{([^}]+)\}', text)),
          'manuscript release path missing')
    check(r'\title{Loop detection in graph-geodesic persistent homology}' in text, 'title')
    check(re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', text) ==
          ['figures/arc_support_fan.png'], 'external figures')
    check((paper / 'figures/arc_support_fan.png').read_bytes() ==
          (ROOT / 'manuscript_figures/arc_support_fan.png').read_bytes(), 'figure asset')
    check(re.findall(r'\\label\{(tab:[^}]+)\}', text) ==
          ['tab:rank-scaling', 'tab:radial-noise', 'tab:smooth-enlargement'], 'table labels')
    appendix = text[text.index(r'\appendix'):]
    check(re.findall(r'\\section\{([^}]+)\}', appendix) == [
        'Smooth-curve locality', 'Chord ties and the arc-length interval',
        'Proof of the tied-rank theorem', 'Comparisons with existing graph frameworks',
        'Finite-size check of radius enlargement', 'Protocols and numerical tolerances'],
        'appendix order')

    spec = importlib.util.spec_from_file_location('detection_sweep', ROOT / 'code/rank_table_sweep.py')
    sweep = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sweep)
    import radial_noise_knn as noise_study

    config = json.loads((ROOT / 'code/config.json').read_text())
    summary = {(int(r['n']),r['profile'],r['method']):r for r in rows(DATA/'summary.csv')}
    check(len(summary) == 24, 'rank summary inventory')
    table = section(text,r'\label{tab:rank-scaling}',r'\end{table}')
    method_names = {'OR_kNN':r'Symmetric $k$NN','mutual_kNN':r'Mutual $k$NN'}
    displayed = {}
    for line in table.splitlines():
        cells = [v.strip() for v in line.rstrip('\\ ').split('&')]
        if len(cells)==5 and cells[0] in method_names.values():
            displayed[cells[0]] = cells[1:]
    for n, reps in sweep.SETTINGS:
        with np.load(DATA / f'samples_{n}.npz',allow_pickle=False) as data:
            ranks, gaps = data['ranks'], data['largest_gap']
        check(ranks.shape==(reps,3,2,2) and gaps.shape==(reps,3), f'array shape {n}')
        check(np.max(gaps)<1/3 and np.all(ranks[...,0]<=ranks[...,1]),f'gap/rank validity {n}')
        for j in (0,reps//2,reps-1):
            expected_ranks, expected_gaps = sweep.sample(n,j)
            check(np.array_equal(ranks[j],expected_ranks) and np.array_equal(gaps[j],expected_gaps),
                  f'seed replay {n}/{j}')
        for p, profile in enumerate(sweep.PROFILES):
            for m, method in enumerate(sweep.METHODS):
                rec=summary[n,profile,method]
                extra=ranks[:,p,m,1]-ranks[:,p,m,0]
                count=int(np.count_nonzero(extra==0))
                check(int(rec['samples'])==reps and int(rec['tied'])==count, f'counts {n}/{profile}/{method}')
                check(abs(float(rec['tie_percent'])-100*count/reps)<1e-11, 'rank summary percentage')
                lo,hi=sweep.wilson(count,reps)
                check(abs(float(rec['wilson_lo_percent'])-100*lo)<1e-11 and
                      abs(float(rec['wilson_hi_percent'])-100*hi)<1e-11, 'Wilson bounds')
                check(100*(hi-lo)/2 < (.7 if reps==20000 else 1.4), 'reported interval widths')
                check(all(int(rec[f'plus{j}'])==int(np.count_nonzero(extra<=j)) for j in (1,2,3)),
                      'extra-rank counts')
                if profile != "uniform":
                    continue
                cell=displayed[method_names[method]][[128,512,4096,32768].index(n)]
                check(cell==f'{100*count/reps:.1f}', f'Table 1 entry {n}/{profile}/{method}')
    check(len(displayed)==2, 'Table 1 rows')
    event_records=rows(ROOT/'results/clean_closure/events.csv')
    check(len(event_records)==660,'clean-circle event records')
    for method in ('R2_global','OR_kNN','mutual_kNN'):
        for event in ('first_connected','pre_completion','completion'):
            group=[r for r in event_records if r['method']==method and r['event']==event]
            check(len(group)==60 and len({(r['seed'],r['profile']) for r in group})==60,
                  'clean-circle sample inventory')
            check(all(r['connected']=='True' and int(float(r['h1_bars']))==(1 if event=='completion' else 0) for r in group),
                  'clean-circle event outcomes')
    for method, expected in (('OR_kNN', 22), ('mutual_kNN', 20)):
        group = [r for r in event_records
                 if r['method'] == method and r['event'] == 'smallest_connected_k']
        check(len(group) == 60 and all(r['connected'] == 'True' for r in group)
              and sum(float(r['h1_bars']) > 0 for r in group) == expected,
              f'Appendix F first connected integer k: {method}')
    clouds=rows(ROOT/'results/clean_closure/clouds.csv')
    check(len(clouds)==60,'clean-circle cloud inventory')
    for cloud in clouds:
        check(0<float(cloud['g_second'])<=float(cloud['g_max'])<1/3,'clean-circle gaps')
        sample=[r for r in event_records if r['seed']==cloud['seed'] and r['profile']==cloud['profile']]
        check(len(sample)==11,'clean-circle events per sample')
        for method in ('R2_global','OR_kNN','mutual_kNN'):
            budget={r['event']:int(r['M']) for r in sample if r['method']==method}
            check(budget['first_connected']<=budget['pre_completion'] and
                  budget['pre_completion']+1==budget['completion'],'clean-circle event budgets')
    large={(int(r['n']),r['method']):r for r in rows(DATA/'uniform_100k/rank_sweep.csv')}
    for n,reported in ((128,1.74),(512,.94),(4096,.27),(32768,.015)):
        record=large[n,'mutual_kNN']
        difference=100*(float(record['tie_rate'])-float(record['formula']))
        tolerance=.0005 if n==32768 else .005
        check(int(record['samples'])==100000 and abs(difference-reported)<tolerance,
              f'100,000-sample mutual discrepancy at {n}')
        check(abs(difference-float(record['difference_percentage_points']))<1e-11,
              '100,000-sample discrepancy record')
    check('and $0.015$ percentage' in text,'reported largest-sample discrepancy')
    events=rows(ROOT/'results/clean_closure/event_table.csv')
    for method in ('R2_global','OR_kNN','mutual_kNN'):
        for event in ('first_connected','pre_completion','completion'):
            record=next(r for r in events if r['method']==method and r['event']==event)
            check(int(record['attempted'])==int(record['connected'])==60 and
                  int(record['empty'])==(0 if event=='completion' else 60) and
                  int(record['nonempty'])==(60 if event=='completion' else 0),
                  f'clean full-diagram check {method}/{event}')
    audit=sweep.audit_local(29100000,24)
    check(audit['samples']==24 and audit['full_connectivity_checks']==48,'local/full rank replay')

    clean_protocol = section(text, 'For the uniform samples we also counted', r'\FloatBarrier')
    consecutive = config['clean_circle']['n'] - 1
    nonconsecutive = []
    for method in ('OR_kNN', 'mutual_kNN', 'R2_global'):
        budgets = [int(r['M']) for r in event_records if r['method'] == method
                   and r['profile'] == 'uniform' and r['event'] == 'first_connected']
        check(len(budgets) == 20, f'Appendix F uniform prefix counts: {method}')
        share = median((budget - consecutive) / budget for budget in budgets)
        nonconsecutive.append(f'{100 * share:.1f}')
    check(re.findall(r'\$([0-9.]+)\\%\$', clean_protocol) == nonconsecutive,
          'Appendix F median non-consecutive edge shares')

    noise_protocol = config['radial_noise']
    noise_text = section(text, r'\label{sec:radial-noise}', r'\section{')
    noise_bounds = re.search(r'independent uniform draws on\s+\$\[(-?\d+),\s*(-?\d+)\]\$', noise_text)
    check(noise_bounds and noise_protocol['distribution'] ==
          {'family': 'uniform', 'interval': list(map(int, noise_bounds.groups()))},
          'radial-noise distribution in the manuscript and configuration')
    check(noise_protocol['n'] == noise_study.N
          and noise_protocol['amplitudes'] == list(noise_study.AMPLITUDES)
          and noise_protocol['seeds_inclusive'] == [noise_study.SEEDS.start, noise_study.SEEDS.stop - 1]
          and noise_protocol['uniform_draw_seed_offset'] == noise_study.graph.NOISE_SEED_OFFSET
          and noise_protocol['cutoffs_of_noise_free_reference_lifetime'] == [.25, .5],
          'radial-noise reproduction protocol')
    median_limit = re.search(r'median\s+longest interval is less than \$([0-9.]+)\\%\$', noise_text)
    check(median_limit, 'reported radial-noise median lifetime bound')
    noise_seeds = set(noise_study.SEEDS)
    noise=rows(ROOT/'results/radial_noise_knn.csv')
    table=section(text,r'\label{tab:radial-noise}',r'\end{table}')
    observed={}
    for line in table.splitlines():
        cells=[v.strip().strip('$') for v in line.rstrip('\\ ').split('&')]
        if len(cells)==7 and all(v.isdecimal() for v in cells[1:]):
            observed[float(cells[0])]=list(map(int,cells[1:]))
    check(len(noise)==3200 and len(observed)==4,'radial-noise inventory')
    for record in noise:
        ratio=float(record['longest_ratio'])
        check(int(record['empty'])==int(ratio==0) and
              int(record['at_least_quarter'])==int(ratio>=.25) and
              int(record['at_least_half'])==int(ratio>=.5),'radial-noise lifetime flags')
        check(record['at_least_quarter'] == record['at_least_half'],
              'Table 2 quarter- and half-lifetime cutoffs agree')
        if float(record['delta'])==0:
            check(int(record['empty'])==1-int(record['all_consecutive']),'clean radial-noise baseline')
    winding=rows(ROOT/'results/radial_noise_winding.csv')
    expected_winding={(r['seed'],r['method']):r for r in noise
                      if float(r['delta'])==.01 and int(r['at_least_half'])==1}
    check(len(winding)==len(expected_winding)==258 and
          len({(r['seed'],r['method']) for r in winding})==258,'noise winding inventory')
    for r in winding:
        original=expected_winding[r['seed'],r['method']]
        check(r['k']==original['k'] and original['all_consecutive']=='0' and
              int(r['winding_one_fundamental_cycles'])>0 and
              3<=int(r['shortest_winding_one_cycle_vertices'])<512,'noise winding witness')
    baseline={(r['seed'],r['method']):r for r in noise if float(r['delta'])==0}
    check(set(observed) == set(noise_study.AMPLITUDES), 'Table 2 noise amplitudes')
    for delta in noise_study.AMPLITUDES:
        actual=[]
        for method in sweep.METHODS:
            group=[r for r in noise if float(r['delta'])==delta and r['method']==method]
            check(len(group) == len(noise_seeds) and {int(r['seed']) for r in group} == noise_seeds,
                  'radial pairing')
            actual += [sum(int(r[k]) for r in group) for k in ('empty','at_least_half','all_consecutive')]
            if delta == .002:
                ratio = median(float(r['longest_ratio']) for r in group)
                check(100 * ratio < float(median_limit[1]),
                      f'radial-noise median lifetime ratio: {method}')
            if delta==.01:
                changed=sum(int(r['at_least_half']) and not int(baseline[r['seed'],method]['at_least_half']) for r in group)
                check(changed==({'OR_kNN':80,'mutual_kNN':86}[method]),'paired noise transitions')
        check(observed[delta]==actual,f'Table 2 at {delta}')

    peanut=rows(ROOT/'results/first_nonempty/trials.csv')
    closure=rows(ROOT/'results/peanut_neck/trials.csv')
    for a in (0.4,0.6,0.8,0.9):
        for method in ('R2_global','OR_kNN','mutual_kNN'):
            first=[r for r in peanut if float(r['alpha'])==a and r['method']==method]
            near=[r for r in closure if float(r['alpha'])==a and r['method']==method]
            check(len(first)==len(near)==20,'peanut sample counts')
            expected=20 if a<.9 else {'R2_global':20,'OR_kNN':13,'mutual_kNN':18}[method]
            check(sum(r['at_closure']=='True' for r in first)==expected,'peanut first detection')
            check(all(r['short_hops_first']=='False' for r in first if r['early']=='True'),
                  'peanut early stops')
            supported=20 if a<.9 else {'R2_global':0,'OR_kNN':0,'mutual_kNN':2}[method]
            check(sum(r['closure_supported']=='True' for r in near)==supported,'peanut support')
            if a==.9:
                count={'R2_global':20,'OR_kNN':16,'mutual_kNN':20}[method]
                check(sum(int(r['closure_bars'])==1 for r in near)==count,'peanut closure intervals')
    # Cross-check every saved first and closure diagram and its event record.
    closure_by_key={(float(r['alpha']),r['seed'],r['method']):r for r in closure}
    check(len(peanut)==len(closure_by_key)==300,'peanut record inventory')
    for r in peanut:
        key=(float(r['alpha']),r['seed'],r['method'])
        near=closure_by_key[key]
        conn,close,first=map(int,(r['M_conn'],r['M_close'],r['M_first']))
        check((conn,close)==(int(near['M_conn']),int(near['M_close'])) and conn<=first<=close,
              'peanut event budgets')
        check(int(r['prefix_queries'])==first-conn+1 and
              (r['at_closure']=='True')==(first==close) and (r['early']=='True')==(first<close) and
              r['short_hops_first']==near['short_hops_first'],'peanut event flags')
        stem=f"a{float(r['alpha']):g}_seed{r['seed']}_{r['method']}.npz"
        with np.load(ROOT/'results/first_nonempty/checkpoints'/stem,allow_pickle=False) as z:
            first_bars=z['first_bars'].copy()
            check(len(first_bars)==int(r['first_bars']) and len(first_bars)>0,'saved first diagram count')
            check(np.all(first_bars[:,1]>first_bars[:,0]) and
                  abs(float(np.max(first_bars[:,1]-first_bars[:,0]))-float(r['first_longest_lifetime']))<1e-14,
                  'saved first diagram lifetimes')
            if r['method']=='R2_global':
                check(len(z['double_bars'])==int(r['double_bars']),'saved doubled-radius diagram count')
        if first==close:
            stem=f"a{float(r['alpha']):g}_seed{r['seed']}.npz"
            with np.load(ROOT/'results/peanut_neck/checkpoints'/stem,allow_pickle=False) as z:
                close_bars=z[f"{r['method']}_closure_bars"]
                order=lambda a: a[np.lexsort((a[:,1],a[:,0]))]
                check(first_bars.shape==close_bars.shape and
                      np.allclose(order(first_bars),order(close_bars),rtol=0,atol=5e-8),
                      'first and closure diagrams agree at the same event')
    for r in closure:
        stem=f"a{float(r['alpha']):g}_seed{r['seed']}.npz"
        with np.load(ROOT/'results/peanut_neck/checkpoints'/stem,allow_pickle=False) as z:
            for event in ('connection','closure'):
                bars=z[f"{r['method']}_{event}_bars"]
                last_death=float(np.max(bars[:,1])) if len(bars) else 0.
                check(len(bars)==int(r[event+'_bars']) and np.all(bars[:,1]>bars[:,0]) and
                      abs(last_death-float(r[event+'_last_death']))<1e-14,
                      f'saved peanut {event} diagram')
            birth=float(np.min(bars[:,0])) if len(bars) else 0.
            check(abs(birth-float(r['closure_birth']))<1e-14 and
                  (int(r['closure_unsupported_edges'])==0)==(r['closure_supported']=='True'),
                  'saved peanut closure details')
    radius=[r for r in peanut if r['method']=='R2_global']
    check(len(radius)==100 and all(r['double_ring_present']=='True' and int(r['double_bars'])>0 for r in radius),
          'peanut and circle radius doubling')

    circle={(int(r['n']),round(1+float(r['epsilon']),2)):r for r in rows(ROOT/'results/radius_enlargement/summary.csv')}
    extra={(int(r['n']),float(r['alpha'])):r for r in rows(ROOT/'results/radius_enlargement_extra/summary.csv')}
    curves={(int(r['n']),float(r['alpha']),r['curve']):r for r in rows(ROOT/'results/smooth_curve_enlargement/summary.csv')}
    radius_protocol = config['radius_enlargement']
    check(set(extra) == {(n, 2.) for n in radius_protocol['large_circle_alpha_2_sizes']}
          and all(int(r['samples']) == radius_protocol['large_circle_alpha_2_samples']
                  for r in extra.values()), 'Table 3 larger circle cohorts')
    samples_per_curve_and_size = {(n, curve, int(r['samples'])) for (n, _, curve), r in curves.items()}
    check(samples_per_curve_and_size == {
        (n, curve, radius_protocol['ellipse_and_peanut_samples_per_size'])
        for n in radius_protocol['sizes'] for curve in ('ellipse', 'peanut')
    }, 'Table 3 smooth-curve cohorts')
    total_samples = sum(samples for _, _, samples in samples_per_curve_and_size)
    smooth_protocol = json.loads((ROOT / 'results/smooth_curve_enlargement/config.json').read_text())
    check(total_samples == radius_protocol['ellipse_and_peanut_total_samples']
          == smooth_protocol['new_clouds']
          and 2 * smooth_protocol['new_clouds_per_curve'] == total_samples,
          'smooth-curve cohort totals')
    appendix_e = section(text, r'\label{app:radius-enlargement}', r'\section{')
    check(f'${total_samples:,}$ point sets in all'.replace(',', '{,}') in appendix_e,
          'Appendix E smooth-curve cohort total')
    table=section(text,r'\label{tab:smooth-enlargement}',r'\end{table}').replace('{,}', ',')
    observed={}
    for line in table.splitlines():
        cells=[v.strip().strip('$').replace(',','') for v in line.rstrip('\\ ').split('&')]
        if len(cells)==7:
            try: observed[int(cells[1]),float(cells[0])]=[float(x) for x in cells[2:]]
            except ValueError: pass
    check(len(observed)==12 and set(observed)==set(circle),'Table 3 inventory')
    for (n,a), vals in observed.items():
        c=circle[n,a]; e=curves[n,a,'ellipse']; p=curves[n,a,'peanut']
        truth=[100*float(extra.get((n,a),c)['miss_rate']),100*float(e['miss_rate']),
               100*float(p['miss_rate']),100*float(c['spacing_probability']),
               100*float(c['asymptotic_probability'])]
        check(all(abs(x-y)<.000051 for x,y in zip(vals,truth)), f'Table 3 at {n}/{a}')

    local=rows(ROOT/'results/noncircular_locality/locality_check_results.csv')
    local_summary=json.loads((ROOT/'results/noncircular_locality/locality_check_summary.json').read_text())
    check(len(local)==60,'locality sample inventory')
    for curve in ('ellipse','peanut'):
        group=[r for r in local if r['curve']==curve]
        check(len(group)==20,'locality curve samples')
        for method in ('radius','OR','mutual'):
            check(all(r[f'{method}_loc_ok']=='True' and r[f'{method}_shh_ok']=='True' and
                      r[f'{method}_closed_sup']=='True' and int(r[f'{method}_bars_close'])==1 and
                      int(r[f'{method}_bars_end'])==1 for r in group),'locality checks')
    for curve,lo,hi,tol in (('peanut',.1661840819,.1698344462,5.1e-11),
                            ('ellipse',.1942495,.2064328,1e-7)):
        vals=[float(r[f'{m}_fail_chord']) for r in local if r['curve']==curve for m in ('radius','OR','mutual')]
        check(len(vals)==60 and abs(min(vals)-lo)<tol and abs(max(vals)-hi)<tol,'locality failure range')
    check(abs(local_summary['r_loc']['ellipse']['r_loc']-.19188129)<1e-7 and
          abs(local_summary['r_loc']['peanut']['r_loc']-.16618368)<1e-7,'locality grid estimates')
    circle_control=[r for r in local if r['curve']=='circle']
    check(len(circle_control)==20 and all(int(r[f'{m}_M_end'])==512*511//2
          for r in circle_control for m in ('radius','OR','mutual')),'circle support control')

    spots=rows(ROOT/'results/smooth_curve_enlargement/spot_checks.csv')
    check(len(spots)==13050 and all(r['agrees']=='True' and r['connected']=='True'
          and r['r_conn_is_second_chord'] == 'True' for r in spots),
          'smooth-curve diagram spot checks')
    exception=rows(ROOT/'results/smooth_curve_checks/exceptional_diagrams.csv')
    check(len(exception)==123 and all(r['diagram_miss']=='True' and r['table_proxy_agrees']=='True' for r in exception),
          'exceptional smooth-curve diagrams')
    spot_keys = {(r['curve'], int(r['n']), int(r['cloud'])) for r in spots}
    exceptional_keys = {(r['curve'], int(r['n']), int(r['cloud'])) for r in exception}
    check(spot_keys.isdisjoint(exceptional_keys) and '$13{,}050$ full-diagram spot checks' in appendix_e,
          'Appendix E separate full-diagram spot checks')
    for curve,expected in (('ellipse',2),('peanut',39)):
        group=[r for r in exception if r['curve']==curve]
        check(len({r['cloud'] for r in group})==expected and len(group)==3*expected,'exception sample counts')
        for cloud in {r['cloud'] for r in group}:
            with np.load(ROOT/f'results/smooth_curve_checks/exceptional_diagrams/{curve}_{cloud}.npz',allow_pickle=False) as z:
                check(all(len(z[f'bars_{a}'])==0 for a in ('1.25','1.5','2')),'exception saved bars')
    interp=json.loads((ROOT/'results/smooth_curve_checks/interpolation.json').read_text())['results']
    check(len(interp)==2 and all(r['max_coordinate_error']<5.4e-12 and r['max_chord_error']<4.8e-12 and
          r['changed_enlargement_classifications']==0 and r['comparisons']==3000 for r in interp),
          'interpolation checks')
    for result,script in (('exceptional_summary.json','check_exceptional.py'),
                          ('interpolation.json','check_interpolation.py')):
        metadata=json.loads((ROOT/'results/smooth_curve_checks'/result).read_text())
        expected=hashlib.sha256((ROOT/'code'/script).read_bytes()).hexdigest()
        check(metadata['verification_script_sha256']==expected and
              re.fullmatch(r'[0-9a-f]{64}',metadata['generation_script_sha256']) is not None,
              f'{result} script provenance')
    for ratio, lower, upper in [(4/9, .68503, .68514), (.5, .72134, .72136)]:
        period = math.log(1/ratio)
        frequency = 2*math.pi/period
        mean = (1-ratio)/period
        first = math.sqrt(math.pi*frequency/math.sinh(math.pi*frequency))
        radius = 2*mean*first/(1-math.sqrt(2)*math.exp(-math.pi*frequency/2))
        check(lower < mean-radius and mean+radius < upper, 'miss probability envelope')
    examples = json.loads((ROOT / 'results/small_examples_verification.json').read_text())
    checker = ROOT / 'code/verify_small_examples.py'
    check(examples['checker_sha256'] == hashlib.sha256(checker.read_bytes()).hexdigest(),
          'small-example checker provenance')
    shortcut = examples['dense_shortcut']
    check(shortcut['gap'] == '282/1000' and shortcut['betti_one'] == 2 and
          shortcut['tetrahedra'] == [[1, 2, 7, 8], [5, 6, 7, 8], [3, 4, 5, 6]] and
          shortcut['extra_edges'] == [[1, 4], [2, 3]], 'eight-point shortcut certificate')
    pdf=paper/'combined.pdf'
    manuscript_record=json.loads((ROOT/'manuscript_record.json').read_text())
    check(hashlib.sha256((paper/'combined.tex').read_bytes()).hexdigest()==manuscript_record['source_sha256'],
          'manuscript source identity')
    check(hashlib.sha256(pdf.read_bytes()).hexdigest()==manuscript_record['pdf_sha256'],
          'manuscript PDF identity')
    check(hashlib.sha256((paper/'figures/arc_support_fan.png').read_bytes()).hexdigest()==manuscript_record['figure_sha256'],
          'manuscript figure identity')
    pdfinfo=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
    pdftext=subprocess.check_output(['pdftotext',str(pdf),'-'],text=True)
    check(int(re.search(r'^Pages:\s+(\d+)',pdfinfo,re.M)[1]) == manuscript_record['pages'],'PDF page count')
    check('Loop detection in graph-geodesic persistent homology' in pdftext and
          '0.68503' in pdftext and '0.72136' in pdftext,'PDF identity')
    package=ROOT.parent
    overleaf=package/'focused-overleaf'
    check((overleaf/'main.tex').read_bytes()==(paper/'combined.tex').read_bytes(),'Overleaf source synchronization')
    check((overleaf/'figures/arc_support_fan.png').read_bytes()==(paper/'figures/arc_support_fan.png').read_bytes(),
          'Overleaf figure synchronization')
    with zipfile.ZipFile(package/'focused-overleaf.zip') as z:
        check(len(z.namelist())==3 and
              set(z.namelist())=={'main.tex','README.md','figures/arc_support_fan.png'},'Overleaf zip contents')
        for name in z.namelist():
            check(z.read(name)==(overleaf/name).read_bytes(),'Overleaf zip synchronization')
    manifest=json.loads((ROOT/'manifest.json').read_text())
    actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p!=ROOT/'manifest.json'
            and '__pycache__' not in p.parts and p.suffix!='.pyc'}
    check(set(manifest['sha256'])==actual,'release manifest inventory')
    for name,want in manifest['sha256'].items():
        check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want,f'release checksum: {name}')
    print('PASS: references, figure and appendix inventory, 8 Table 1 cells, 12 seed replays,')
    print('      24 local/full-rank rechecks, clean-circle and 100,000-sample checks,')
    print('      radial-noise, locality, peanut and radius checks, PDF identity,')
    print('      Overleaf synchronization and the complete release manifest.')
    print('This checks the reported evidence and document consistency, not the proofs.')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--manuscript-dir',type=Path,default=ROOT.parent/'main-paper')
    main(parser.parse_args().manuscript_dir.resolve())
