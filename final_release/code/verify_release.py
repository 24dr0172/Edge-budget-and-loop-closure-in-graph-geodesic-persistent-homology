"""Verify the release against its saved data and, optionally, the manuscript.

Modes:
  default        all saved-data checks, plus the manuscript if it is found
  --data-only    all saved-data checks, no manuscript
  --paper-only   only the checks behind the paper's claims, plus every file hash
Recorded seeds and rank-only results are replayed. No sample cohort is
generated and the full persistence studies are not rerun. Small independent
replays cover the peanut/trough failure cases, a tiny positive-edge case and
the seven-point example (by F2 reduction).
"""
from pathlib import Path
import argparse,hashlib,json,sys
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import graph_methods as m
from finite_diagram import reduction_h1,matching_distance
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from external_manuscript import locate_manuscript, manuscript_text

# Figure assets kept for the additional studies (see results/additional_figures/README.md).
ADDITIONAL_FIGURES=('clean_budget_curves','noise_budget_curves','extra_neighbor_distribution',
                    'radial_first_connectivity','noise_amplitude_errors','unsupported_shortcut')


def verify_manuscript_identity(root, require, paper):
    """Check the recorded hashes against the supplied manuscript source, PDF and inputs."""
    from pypdf import PdfReader
    data = root / 'results'
    record_file = json.loads((data / 'manuscript_record.json').read_text())
    require('manuscript_portable_path_base', record_file.get('path_base') == 'manuscript_directory')
    for kind, filename in [('source','combined.tex'),('pdf','combined.pdf')]:
        path=paper/filename
        require('manuscript_'+kind+'_present',path.is_file())
        record=record_file[kind]
        require('manuscript_'+kind+'_record_path',record['file']==filename and (paper/record['file']).resolve()==path.resolve())
        require('manuscript_'+kind+'_sha256',record['sha256']==hashlib.sha256(path.read_bytes()).hexdigest())
    for kind in ('inputs','graphics'):
        for record in record_file.get(kind,[]):
            path=(paper/record['file']).resolve()
            require('manuscript_'+kind+'_sha256_'+path.name,path.is_file() and paper.resolve() in path.parents and
                    record['sha256']==hashlib.sha256(path.read_bytes()).hexdigest())
    pages=len(PdfReader(paper/'combined.pdf').pages)
    require('manuscript_identity',record_file['pdf']['pages']==pages,{'pdf_pages':pages})


def verify_noise_coordinates(root, config, require):
    """Replay the seeded noise inputs independently of the exporter."""
    n = config['n']
    path = root / 'results/matched_noise/coordinates.npz'
    expected = {f'{seed}_{field}' for seed in config['seeds']
                for field in ['latent', 'clean', 'normal', 'signed_draw']}
    expected |= {f'{seed}_delta{delta:g}_observed' for seed in config['seeds']
                 for delta in config['delta']}
    max_difference = 0.
    with np.load(path, allow_pickle=False) as coords:
        require('noise_coordinate_inventory', set(coords.files) == expected)
        for seed in config['seeds']:
            s, clean = m.sample_cloud(seed, 'uniform', n=n)
            normal = clean * (2 * np.pi)
            draw = np.random.default_rng(seed + m.NOISE_SEED_OFFSET).uniform(-1., 1., n)
            for field, want in [('latent', s), ('clean', clean),
                                ('normal', normal), ('signed_draw', draw)]:
                got = coords[f'{seed}_{field}']
                require(f'noise_input_replay_{seed}_{field}',
                        got.shape == want.shape and np.array_equal(got, want))
            for delta in config['delta']:
                got = coords[f'{seed}_delta{delta:g}_observed']
                want = clean + delta * draw[:, None] * normal
                require(f'noise_observed_replay_{seed}_{delta:g}',
                        got.shape == (n, 2) and np.array_equal(got, want))
                # Also compare with the sampler's direct normal displacement.
                _, sampled = m.sample_cloud(seed, 'uniform', n=n, sigma=delta)
                max_difference = max(max_difference, float(np.max(np.abs(got - sampled))))
    require('noise_coordinates', max_difference < 1e-15,
            {'seeds': len(config['seeds']), 'points_per_sample': n,
             'amplitudes': len(config['delta']), 'arrays': len(expected),
             'observed_clouds': len(config['seeds']) * len(config['delta']),
             'maximum_direct_sampler_difference': max_difference})


def verify_curve_order_records(root, require):
    """Re-read the curve-order check (Example 1 remark) and replay its key graphs."""
    from verify_curve_order import COLUMNS, TOL, analyse
    data = root / 'results'
    summary = json.loads((data / 'curve_order/summary.json').read_text())
    rows = pd.read_csv(data / 'curve_order/curve_order_checks.csv', dtype={'seed': str},
                       float_precision='round_trip')
    script = root / 'code/verify_curve_order.py'
    require('curve_order_script_hash', summary['script'] == 'code/verify_curve_order.py' and
            summary['script_sha256'] == hashlib.sha256(script.read_bytes()).hexdigest())
    require('curve_order_tolerance', summary['tolerance'] == TOL == 1e-10)
    require('curve_order_columns', list(rows.columns) == COLUMNS)
    # Scope: every saved post-closure circle prefix and every ellipse cell.
    prefixes = pd.read_csv(data / 'useful_diagram/saved_prefixes.csv', dtype={'seed': str})
    cells = pd.read_csv(data / 'ellipse/cells.csv', dtype={'seed': str})
    key = ['seed', 'profile', 'method', 'M']
    circle, ellipse = rows[rows.curve == 'circle'], rows[rows.curve == 'ellipse']
    joined = circle.merge(prefixes, on=key, how='left', suffixes=('', '_saved'), validate='many_to_one')
    require('curve_order_circle_scope', len(circle) == 3439 and (circle.kind == 'closure').sum() == 180 and
            (circle.kind == 'grid').sum() == int(prefixes.is_saved_budget.sum()) and
            joined.death_saved.notna().all() and (joined.M_close == joined.M_close_saved).all() and
            np.allclose(joined.birth, joined.birth_saved, rtol=0, atol=1e-12) and
            np.allclose(joined.death, joined.death_saved, rtol=0, atol=1e-12))
    require('curve_order_ellipse_scope', len(ellipse) == len(cells) == 3516 and
            sorted(map(tuple, ellipse[key].values.tolist())) == sorted(map(tuple, cells[key].values.tolist())))
    require('curve_order_n', (rows.n == 512).all() and (rows.pairs == 512 * 511 // 2).all())
    # Internal consistency of the per-graph records.
    cert, before, alive = rows.certified_failure, rows.certified_before_death, rows.certified_while_alive
    require('curve_order_record_consistency',
            (cert == (rows.certified_pairs > 0)).all() and (cert == rows.first_certified_threshold.notna()).all() and
            (before <= cert).all() and (alive <= before).all() and
            (rows.first_certified_before_death[before] < rows.death[before]).all() and
            (rows.first_certified_while_alive[alive] >= rows.birth[alive]).all() and
            (rows.last_certified_while_alive[alive] <= rows.death[alive]).all() and
            (rows.certified_alive_measure[alive] > 0).all() and
            (rows.certified_alive_measure[alive] <= (rows.death - rows.birth)[alive]).all() and
            (rows.certified_alive_windows[~alive] == 0).all() and
            (rows.shorter_arc_failure == (rows.shorter_arc_failing_pairs > 0)).all())
    def span(values):
        values = values.dropna()
        return [float(values.min()), float(values.max())] if len(values) else None
    for curve, g in [('circle', circle), ('ellipse', ellipse)]:
        want = summary[curve]
        got = dict(graphs=len(g), certified_failure_graphs=int(g.certified_failure.sum()),
                   certified_before_death_graphs=int(g.certified_before_death.sum()),
                   certified_while_alive_graphs=int(g.certified_while_alive.sum()),
                   shorter_arc_failure_graphs=int(g.shorter_arc_failure.sum()),
                   shorter_arc_before_death_graphs=int(g.shorter_arc_before_death.sum()))
        ranges = dict(first_certified_threshold_range=span(g.first_certified_threshold),
                      first_certified_before_death_range=span(g.first_certified_before_death),
                      death_range=span(g.death), shorter_arc_first_threshold_range=span(g.shorter_arc_first_threshold))
        require(f'curve_order_summary_{curve}',
                all(want[k] == v for k, v in got.items()) and
                all((want[k] is None and v is None) or np.allclose(want[k], v, rtol=0, atol=1e-15)
                    for k, v in ranges.items()), {**got, **ranges})
    listed = summary['ellipse_shorter_arc_before_death']
    flagged = ellipse[ellipse.shorter_arc_before_death]
    require('curve_order_listed_graphs', sorted((r['seed'], r['profile'], r['method'], r['M']) for r in listed) ==
            sorted(map(tuple, flagged[key].values.tolist())))
    # Replay the arc-supported single-bar ellipse graphs from coordinates.
    replay = flagged[(flagged.arc_supported == True) & (flagged.bars == 1)]
    with np.load(data / 'ellipse_checks/coordinates.npz', allow_pickle=False) as coords:
        for r in replay.itertuples():
            s, xy = coords[f'{r.seed}_{r.profile}_latent'], coords[f'{r.seed}_{r.profile}_xy']
            _, iu, ju, length, orders, _ = m.build_orders(xy)
            D = m.graph_metric(orders[r.method], iu, ju, length, int(r.M), len(xy))
            out = analyse(D, s, r.birth, r.death)
            require(f'curve_order_replay_{r.seed}_{r.profile}_{r.method}_M{r.M}',
                    out['certified_pairs'] == r.certified_pairs and
                    abs(out['first_certified_threshold'] - r.first_certified_threshold) < 1e-15 and
                    out['certified_while_alive'] == r.certified_while_alive and
                    out['shorter_arc_failing_pairs'] == r.shorter_arc_failing_pairs)
    if 'two_sat' in summary:
        exact = summary['two_sat']['graphs']
        require('curve_order_two_sat', len(exact) == len(replay) and
                all(not sat for g in exact for sat, local in zip(g['satisfiable'], g['locally_certified']) if local),
                {'graphs': len(exact), 'thresholds': sum(len(g['thresholds']) for g in exact),
                 'unsatisfiable': sum(not s for g in exact for s in g['satisfiable'])})
    require('curve_order_check', summary['status'] == 'passed',
            {k: {f: summary[k][f] for f in ('graphs', 'certified_failure_graphs', 'certified_before_death_graphs',
                                            'certified_while_alive_graphs', 'shorter_arc_failure_graphs')}
             for k in ('circle', 'ellipse')} | {'replayed_ellipse_graphs': len(replay)})


def verify_payload(root,require,required_schemas=None):
    data=root/'results'
    # Keep mutable execution logs and this report outside the scientific manifest.
    manifest=data/'file_manifest.json'
    require('manifest_present',manifest.is_file())
    if manifest.exists():
        record=json.loads(manifest.read_text());entries=record['sha256']
        require('manifest_exclusion_policy',record['excluded']==[
            'results/file_manifest.json','results/verification_report.json','results/run_logs/'])
        actual={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file()
                and f not in {manifest,data/'verification_report.json'}
                and data/'run_logs' not in f.parents and '__pycache__' not in f.parts and f.suffix!='.pyc'}
        require('manifest_file_inventory',set(entries)==actual)
        for name,want in entries.items():
            require('checksum_'+name,hashlib.sha256((root/name).read_bytes()).hexdigest()==want)
    import csv
    schemas=json.loads((root/'code/result_schemas.json').read_text())
    selected = schemas if required_schemas is None else {name:schemas[name] for name in required_schemas}
    for name,columns in selected.items():
        path=data/name
        require('schema_file_'+name,path.is_file())
        with path.open() as stream:require('schema_columns_'+name,next(csv.reader(stream))==columns)
    for path in list((root/'code').glob('*.py'))+list((root/'tools').glob('*.py')):
        compile(path.read_text(),str(path),'exec')
    require('code_syntax',True)


def verify_documentation(root,require):
    """Every relative Markdown link in the guides resolves inside this release."""
    import re
    guides=sorted(root.rglob('*.md'))
    require('documentation_guides_present',{str(g.relative_to(root)) for g in guides}>={'README.md','code/README.md','results/README.md','tools/README.md','manuscript_figures/README.md','results/comparisons/README.md'})
    links=0
    for guide in guides:
        text=guide.read_text()
        for target in re.findall(r'\]\((?:<([^>]+)>|([^)\s]+))\)',text):
            target=(target[0] or target[1]).split('#',1)[0]
            if not target or re.match(r'[a-z]+:',target):continue
            path=(guide.parent/target).resolve()
            require('documentation_link_'+str(guide.relative_to(root))+'_'+target,path.exists() and (path==root or root in path.parents))
            links+=1
    require('documentation_links_checked',links>0,{'guides':len(guides),'relative_links':links})


def verify(root,manuscript_dir=None,data_only=False,paper_only=False):
    checks={}
    def require(name,condition,detail=None):
        if not bool(condition):raise AssertionError(name)
        checks[name]={'passed':True,**({'detail':detail} if detail is not None else {})}
    data=root/'results';native=set(m.METHODS)
    require('release_structure',{p.name for p in root.iterdir()}=={'code','tools','results','manuscript_figures','README.md','LICENSE','CITATION.cff','Dockerfile'})
    verify_documentation(root,require)
    require('release_no_latex_documents',not list(root.rglob('*.tex')) and not list(root.rglob('*.bib')))
    require('release_no_paper_or_source_bundle',not any(p.name in ('combined.pdf','main.pdf','combined.tex','main.tex','overleaf.zip','manuscript-png.zip') for p in root.rglob('*')))
    from reproduce import stages
    commands=stages(root,1)
    stage_paths=[root/name if name.startswith('tools/') else root/'code'/name for name,args in commands]
    require('generation_stage_files_present',all(path.is_file() for path in stage_paths))
    final_args=dict(commands)
    require('generation_data_only_defaults',final_args['tools/record_current.py']==['--data-only'] and
            final_args['verify_release.py']==['--data-only'])
    paper=None if data_only else locate_manuscript(root,manuscript_dir)
    if paper is not None:
        verify_manuscript_identity(root, require, paper)
        layout=json.loads((data/'manuscript_structure.json').read_text())
        manuscript=manuscript_text(paper)
        require('manuscript_structure_record',layout['status']=='passed',layout)
        environments = ('definition','theorem','proposition','corollary','structurallemma',
                        'appendixlemma','lemma','remark','example','proof')
        for name in environments:
            require('manuscript_environment_'+name,
                    layout['environment_counts'][name]==manuscript.count(r'\begin{'+name+'}'))
        from record_current import manuscript_records
        _, _, replayed_layout = manuscript_records(root,paper)
        require('manuscript_structure_replay',layout==replayed_layout)
        if not paper_only:
            comparison_record=json.loads((data/'comparisons/completion.json').read_text())
            require('comparison_manuscript_source_identity',
                    comparison_record['manuscript_source_sha256']==layout['source_sha256'])
            require('comparison_manuscript_pdf_identity',
                    comparison_record['manuscript_pdf_sha256']==layout['pdf_sha256'])
        require('manuscript_structure_source_identity',hashlib.sha256((paper/'combined.tex').read_bytes()).hexdigest()==layout['source_sha256'])
        require('manuscript_structure_pdf_identity',hashlib.sha256((paper/'combined.pdf').read_bytes()).hexdigest()==layout['pdf_sha256'])
        result_labels=['thm:loop-transition','thm:arc-supported','thm:circle-parameters']
        positions=[manuscript.index('\\label{'+label+'}') for label in result_labels]
        require('manuscript_result_order',positions==sorted(positions) and manuscript.count(r'\begin{theorem}')==3 and
                all(layout['result_numbers']['theorem'][label]==i for i,label in enumerate(result_labels,1)))
        from verify_manuscript_tables import verify_manuscript_tables
        verify_manuscript_tables(root,require,paper)

    from verify_paper_evidence import verify_paper_evidence
    paper_evidence=verify_paper_evidence(root,require,include_additional=not paper_only)
    require('paper_evidence',paper_evidence['status']=='passed',paper_evidence)
    if paper_only:
        verify_payload(root,require,['rank_table.csv','noncircular_locality/locality_check_results.csv',
                                    'retention_extra_layers/retention.csv'])
        return dict(status='passed',requested_scope='paper_saved_data',
            manuscript_comparison=dict(status='passed' if paper is not None else 'not_run',
                                       location='external' if paper is not None else None),
            performed_checks=['Table 3 raw rank counts and export','Small examples',
                              'Fixed-k, bridge-free, heuristic-model and diagram-cost claims',
                              'Extra-layer counts, radius beta integral and independent finite certificates',
                              'Every distributed payload hash and complete inventory','Code syntax']+
                             (['Manuscript identity, structure and tables'] if paper is not None else []),
            additional_scientific_checks='not_run',measurement_generation='not_run',
            checks_passed=len(checks),checks=checks)

    require('additional_figure_assets',all((root/'results/additional_figures'/(name+'.'+ext)).exists() for name in ADDITIONAL_FIGURES for ext in ['pdf','png']))
    cells=pd.read_csv(data/'clean_closure/cells.csv');clouds=pd.read_csv(data/'clean_closure/clouds.csv')
    require('clean_native_panel',len(cells)==4740 and set(cells.method)==native and set(cells.noise)=={'clean'})
    require('saved_sample_scope',len(clouds)==60 and (clouds.n==512).all() and clouds.seed.nunique()==20)
    table=pd.read_csv(data/'clean_closure/event_table.csv')
    require('eleven_event_rows',len(table)==11 and (table.attempted==60).all() and (table.connected==60).all())
    expected=[60,60,0,60,38,60,0,60,40,60,0]
    require('event_empty_counts',table['empty'].to_list()==expected,expected)
    knn_events=pd.read_csv(data/'clean_closure/knn_event_checks.csv')
    for method in ['OR_kNN','mutual_kNN']:
        for event,landmark in [('pre_completion','pre_coverage_last'),('completion','own_coverage')]:
            g=knn_events[(knn_events.method==method)&(knn_events.landmark==landmark)]
            row=table[(table.method==method)&(table.event==event)].iloc[0]
            require(f'knn_closure_table_{method}_{event}',len(g)==60 and g.connected.all() and row['empty']==int((g.est_h1_bars==0).sum()) and row.nonempty==int((g.est_h1_bars>0).sum()))
    control=pd.read_csv(data/'order_checks/shuffled_control.csv')
    require('permuted_scale_control',len(control)==60 and not control.short_hops_first_passes.any() and (control.M_first_connected<control.M_close).all() and (control.first_connected_H1_bars>0).sum()==6 and control.first_nonempty_matches_first_closed_grid.sum()==16,{'samples':60,'condition_passes':0,'loops_before_closure_at_connectivity':6,'grid_nonempty_matches_closure':16})
    require('control_grid_matches',(control.first_nonempty_matches_first_closed_grid==(control.first_nonempty_grid_M==control.first_closed_grid_M)).all())
    for row in control.itertuples():
        with np.load(data/f'order_checks/checkpoints/{row.profile}_{row.seed}.npz',allow_pickle=False) as z:
            require(f'shuffled_saved_connectivity_{row.profile}_{row.seed}',
                    len(z['R3_permuted_first_connected'])==row.first_connected_H1_bars)
            candidates=[(int(key.rsplit('_M',1)[1]),len(z[key])) for key in z.files if key.startswith('R3_permuted_M')]
            positive=[count for count,bars in candidates if bars>0]
            require(f'shuffled_saved_grid_{row.profile}_{row.seed}',min(positive)==row.first_nonempty_grid_M)
            for key in z.files:
                if key.endswith('_first_connected') or '_M' in key:
                    bars=z[key]
                    require(f'shuffled_finite_bars_{row.profile}_{row.seed}_{key}',
                            bars.ndim==2 and bars.shape[1]==2 and np.isfinite(bars).all() and (bars[:,1]>bars[:,0]).all())
    selected=pd.read_csv(data/'clean_closure/closure_selection.csv')
    require('180_first_closed_series',len(selected)==180 and selected.certificate_holds.all() and selected.nonempty_matches_covered.all())
    ft=pd.read_csv(data/'clean_closure/first_closed_errors.csv')
    require('closure_table_reports_all_errors',len(ft)==9 and (ft.samples==20).all() and ft[['median_e_inf','median_e_2','median_e_PH']].notna().all().all())
    for _,row in ft.iterrows():
        g=selected[(selected.method==row.method)&(selected.profile==row.profile)]
        for out,source in [('median_e_inf','first_covered_e_inf'),('median_e_2','first_covered_e_2'),('median_e_PH','first_covered_p1')]:
            require(f'closure_median_{row.profile}_{row.method}_{out}',abs(row[out]-g[source].median())<1e-14)
    mono=pd.read_csv(data/'clean_closure/covered_monotonicity.csv')
    require('180_saved_monotone_series',len(mono)==180 and (mono[['e_inf_decreases','e_2_decreases','p1_decreases']]==0).all().all(),{'series':180,'stored_loss_tolerance':5e-8})
    own=pd.read_csv(data/'clean_closure/own_closure_cells.csv')
    require('own_closure_not_global_split',len(own)==4740 and (own.covered==(own.M>=own.M_cov)).all())
    joined=cells.merge(own,on=['seed','profile','method','budget_label','M'],suffixes=('_a','_b'),validate='one_to_one')
    require('own_closure_rows_align',len(joined)==4740 and (joined.connected_a==joined.connected_b).all())
    grid=cells[cells.budget_label.str.startswith('a')];budget=pd.read_csv(data/'clean_closure/budget_summary.csv')
    require('grid_scope',len(grid)==4320 and len(budget)==216)
    for _,row in budget.iterrows():
        g=grid[(grid.profile==row.profile)&(grid.method==row.method)&(grid.M==row.M)]
        require(f'all_attempts_{row.profile}_{row.method}_{row.M}',len(g)==20 and row.connected==int(g.connected.sum()))
        for out,source in [('median_e_inf','e_inf'),('median_e_2','e_2'),('median_e_PH','p1')]:
            want=float(np.median(np.where(g.connected,g[source],np.inf)));got=row[out]
            require(f'all_sample_median_{row.profile}_{row.method}_{row.M}_{out}',(np.isinf(want) and np.isinf(got)) or abs(want-got)<1e-13)
    ranks=pd.read_csv(data/'clean_closure/rank_landmarks.csv');gaps=pd.read_csv(data/'extra_k/samples.csv')
    require('120_paired_extra_k_rows',len(gaps)==120 and (gaps.extra_k==gaps.k_close-gaps.k_conn).all() and (gaps.extra_k>=0).all())
    require('miss_means_unclosed_integer_layer',(gaps.miss_at_connectivity==(gaps.extra_k>0)).all())
    require('miss_counts_38_and_40',gaps.groupby('method').miss_at_connectivity.sum().to_dict()=={'OR_kNN':38,'mutual_kNN':40})
    es=pd.read_csv(data/'extra_k/summary.csv');expected_gaps={
        ('uniform','OR_kNN'):(13,1.,3),('strong_cosine','OR_kNN'):(15,1.,3),('localized_trough','OR_kNN'):(10,.5,3),
        ('uniform','mutual_kNN'):(10,.5,4),('strong_cosine','mutual_kNN'):(14,1.,3),('localized_trough','mutual_kNN'):(16,7.,26)}
    require('six_extra_k_groups',len(es)==6)
    for _,row in es.iterrows():
        g=gaps[(gaps.profile==row.profile)&(gaps.method==row.method)];want=expected_gaps[(row.profile,row.method)]
        require(f'extra_k_summary_{row.profile}_{row.method}',(row.misses,row.median_extra_k,row.max_extra_k)==want and row.samples==len(g)==20)
    # Replay the saved coordinates from their seeds; this is not a new sample-size experiment.
    max_coordinate_difference=0.;rank_checks=0
    with np.load(data/'clean_closure/latent_coordinates.npz') as latents:
        for profile in m.PROFILES:
            for seed in range(20310101,20310121):
                with np.load(data/f'useful_diagram/checkpoints/{profile}_{seed}.npz') as z:
                    s=z['latent'];xy=z['observed'];reference=z['reference']
                reconstructed,clean=m.sample_cloud(seed,profile)
                diff=max(float(np.max(np.abs(reconstructed-s))),float(np.max(np.abs(clean-xy))))
                max_coordinate_difference=max(max_coordinate_difference,diff)
                require(f'saved_coordinate_replay_{profile}_{seed}',diff<1e-13 and np.max(np.abs(s-latents[f'{seed}_{profile}']))<1e-13)
                require(f'positive_reference_{profile}_{seed}',len(reference)==1 and reference[0,1]>reference[0,0])
                chord,iu,ju,length,orders,keys=m.build_orders(xy)
                cyc=(ju-iu==1)|((iu==0)&(ju==len(s)-1))
                for method in ['OR_kNN','mutual_kNN']:
                    row=ranks[(ranks.seed==seed)&(ranks.profile==profile)&(ranks.method==method)].iloc[0]
                    conn=m.first_connected(orders[method],iu,ju,len(s));closed=m.closure(orders[method],iu,ju,len(s))
                    sorted_keys=np.sort(keys[method][cyc]);k_conn=int(keys[method][orders[method][conn-1]])
                    require(f'rank_landmark_replay_{profile}_{seed}_{method}',conn==row.M_conn==row.M_second_to_last_consecutive and closed==row.M_own_cover and k_conn==row.k_conn==sorted_keys[-2] and sorted_keys[-1]==row.k_largest_consecutive_key)
                    rank_checks+=1
    checks['coordinate_and_rank_replay']={'passed':True,'detail':{'saved_clouds':60,'rank_cases':rank_checks,'max_coordinate_difference':max_coordinate_difference}}
    cert=pd.read_csv(data/'order_checks/certificates.csv')
    require('certificate_counts',cert.groupby('method').holds.sum().to_dict()=={'R1_local':59,'R2_global':60,'OR_kNN':60,'mutual_kNN':60})
    require('CkNN_certificate_only',len(cert)==240 and 'R1_local' not in set(cells.method))
    # Match stored diagrams and their scores, without assuming single-bar reduction.
    prefixes=pd.read_csv(data/'useful_diagram/saved_prefixes.csv');trans=pd.read_csv(data/'useful_diagram/saved_transitions.csv')
    require('useful_diagram_scope',len(prefixes)==3556 and len(trans)==3376 and set(prefixes.method)==native and (prefixes.M>=prefixes.M_close).all())
    require('single_bar_fixed_birth', (prefixes.bars==1).all() and prefixes.birth_deviation.max()<1e-7, float(prefixes.birth_deviation.max()))
    require('small_gap_sufficient_condition',prefixes.g.max()<2/np.pi-.5, float(prefixes.g.max()))
    require('saved_diagram_transitions',not trans[['raw_drop','clipped_drop','crosses_below_baseline']].any().any() and (trans.delta>=-1e-7).all() and (trans.death_delta<=1e-7).all())
    max_score_difference=0.;stored_diagrams=0
    for (profile,seed),g in prefixes.groupby(['profile','seed']):
        with np.load(data/f'useful_diagram/checkpoints/{profile}_{seed}.npz') as z:
            ref=z['reference']
            for _,row in g.iterrows():
                bar=z[f'{row.method}_M{row.M}']
                require(f'stored_bar_{profile}_{seed}_{row.method}_{row.M}',len(bar)==1 and np.max(np.abs(bar[0]-[row.birth,row.death]))<1e-13)
                score=matching_distance(bar,ref)/row.diameter
                max_score_difference=max(max_score_difference,abs(score-row.e_PH));stored_diagrams+=1
    require('stored_diagram_matching',max_score_difference<1e-12,{'diagrams':stored_diagrams,'max_error_difference':max_score_difference})
    witness=json.loads((data/'useful_diagram/rational_counterexample.json').read_text());s=np.array(witness['coordinates'],float)
    true=m.intrinsic(s);ref,tiny= reduction_h1(true)
    require('seven_point_reference_no_filtered_positive_bars',tiny==0.)
    xy=np.c_[np.cos(2*np.pi*s),np.sin(2*np.pi*s)]/(2*np.pi)
    chord,iu,ju,length,orders,keys=m.build_orders(xy)
    for method,record in witness['methods'].items():
        scores=[]
        for stage in ['before','after']:
            row=record[stage];d=m.graph_metric(orders[method],iu,ju,length,row['M'],len(s));bar,tiny=reduction_h1(d)
            require(f'seven_point_{method}_{stage}_no_filtered_positive_bars',tiny==0.)
            score=matching_distance(bar,ref)/true.max();scores.append(score)
            require(f'seven_point_{method}_{stage}',np.max(np.abs(bar-np.array(record['diagrams'][0 if stage=='before' else 1])))<1e-12 and abs(score-row['e_PH'])<1e-12)
        require(f'capped_witness_{method}',scores[1]<scores[0] and abs(min(scores[0],record['before']['P_star'])-min(scores[1],record['before']['P_star']))<1e-12)
    noise=pd.read_csv(data/'matched_noise/cells.csv');ns=pd.read_csv(data/'matched_noise/summary.csv');config=json.loads((data/'matched_noise/config.json').read_text())
    verify_noise_coordinates(root, config, require)
    require('matched_noise_scope',len(noise)==240 and set(noise.event)=={'first_connected'} and (noise.n==512).all() and len(ns)==12 and noise.seed.nunique()==20)
    require('no_long_first_connected_bars', (ns.long_at_half==0).all())
    for _,row in ns.iterrows():
        g=noise[(noise.method==row.method)&(noise.delta==row.delta)]
        want_empty=20 if row.delta==0 else 9 if row.delta==.0005 else 0
        require(f'noise_summary_{row.method}_{row.delta}',len(g)==20 and g.empty_H1.sum()==row.empty_H1==want_empty and abs(g.longest_over_reference.median()-row.median_longest_over_reference)<1e-13)
    for seed in config['seeds']:
        g=noise[noise.seed==seed]
        with np.load(data/f'matched_noise/diagrams/diagrams_{seed}.npz') as z:
            require(f'noise_diagram_scope_{seed}',len(z.files)==13)
            ref=z[f'oracle_{seed}'];life=float(np.max(ref[:,1]-ref[:,0]))
            for _,row in g.iterrows():
                bar=z[f'{seed}_delta{row.delta:g}_{row.method}_first_connected'];l=float(np.max(bar[:,1]-bar[:,0])) if len(bar) else 0.
                require(f'noise_diagram_{seed}_{row.method}_{row.delta}',len(bar)==row.H1_bars and abs(l/life-row.longest_over_reference)<1e-12)
    closure_noise=pd.read_csv(data/'matched_noise/closure_comparison.csv')
    require('noise_closure_comparison',len(closure_noise)==80 and (closure_noise.delta==.01).all() and closure_noise.connected.all() and not closure_noise.empty_H1.any() and closure_noise.stability_ok.all() and (closure_noise.longest_over_reference>=.5).all(),{'observations':80,'before_closure_long_per_method':20,'closure_nonempty_per_method':20})
    require('noise_closure_event_positions',((closure_noise.M==closure_noise.M_cov)&(closure_noise.event=='coverage') | (closure_noise.M==closure_noise.M_cov-1)&(closure_noise.event=='pre_coverage')).all())
    with np.load(data/'matched_noise/closure_diagrams.npz') as z:
        require('noise_closure_diagram_inventory',len(z.files)==100)
        for row in closure_noise.itertuples():
            reference=z[f'oracle_{row.seed}'];bars=z[f'{row.seed}_delta0.01_{row.method}_{row.event}']
            ratio=float(np.max(bars[:,1]-bars[:,0])/np.max(reference[:,1]-reference[:,0]))
            require(f'noise_closure_diagram_{row.seed}_{row.method}_{row.event}',len(bars)==row.H1_bars and abs(ratio-row.longest_over_reference)<1e-12)
    noise_summary=pd.read_csv(data/'matched_noise/closure_summary.csv')
    require('noise_closure_summary',len(noise_summary)==4 and (noise_summary.attempted==20).all() and (noise_summary.empty_H1==0).all() and (noise_summary.long_at_half==20).all())
    ellipse=pd.read_csv(data/'ellipse_checks/checks.csv');constants=json.loads((data/'ellipse_checks/analytic_constants.json').read_text())
    from experiment_core import ellipse_points
    with np.load(data/'ellipse_checks/coordinates.npz',allow_pickle=False) as coordinates:
        for seed,profile in ellipse[['seed','profile']].drop_duplicates().itertuples(index=False,name=None):
            latent,_=m.sample_cloud(seed,profile)
            expected=ellipse_points(latent)
            require(f'ellipse_seed_coordinates_{profile}_{seed}',
                    np.array_equal(coordinates[f'{seed}_{profile}_latent'],latent) and
                    np.max(np.abs(coordinates[f'{seed}_{profile}_xy']-expected))<1e-12)
    require('ellipse_scope',len(ellipse)==180 and set(ellipse.method)==native and ellipse.certificate_passes.all() and ellipse.saved_bar_spans_bound.all())
    require('ellipse_constants',abs(constants['global_chord_arc_lower_factor']-1/np.pi)<1e-15 and constants['normal_chord_lower_after_normalization']>.18 and (ellipse.locality_radius_bound<.18).all() and (ellipse.B<1/(3*np.pi)).all())
    require('ellipse_loop_interval',(ellipse.analytic_lifetime_lower>0).all() and np.allclose(ellipse.guaranteed_death_lower,1/(3*np.pi)))
    with np.load(data/'ellipse_checks/diagrams.npz') as z:
        expected_diagrams={f'{r.seed}_{r.profile}_{r.method}_M{M}' for r in ellipse.itertuples() for M in [r.M_conn,r.M_close-1,r.M_close]}
        require('ellipse_saved_diagram_inventory',set(z.files)==expected_diagrams,{'distinct_diagrams':len(z.files),'event_checks':540,'coincident_events':1})
        for _,row in ellipse.iterrows():
            prefix=f'{row.seed}_{row.profile}_{row.method}_M'
            require(f'ellipse_preclosure_{row.seed}_{row.profile}_{row.method}',row.M_conn==row.M_second_to_last_consecutive and row.M_conn<row.M_close and len(z[prefix+str(row.M_conn)])==0 and len(z[prefix+str(row.M_close-1)])==0)
            bars=z[prefix+str(row.M_close)];tol=row.endpoint_tolerance
            require(f'ellipse_closure_{row.seed}_{row.profile}_{row.method}',np.any((bars[:,0]<=row.guaranteed_birth_upper+tol)&(bars[:,1]>=row.guaranteed_death_lower-tol)))
    fixed=pd.read_csv(data/'extra_k/rank_sweep/fixed_k_summary.csv')
    fixed_source=pd.read_csv(data/'extra_k/rank_sweep/samples.csv')
    expected_scope={('uniform',method,k,n) for method in ['OR_kNN','mutual_kNN']
                    for k in [10,15] for n in [128,512,4096,32768]}
    require('fixed_k_summary_scope',len(fixed)==16 and
            set(zip(fixed.profile,fixed.method,fixed.k,fixed.n))==expected_scope)
    for row in fixed.itertuples():
        group=fixed_source[(fixed_source.profile==row.profile)&(fixed_source.method==row.method)&(fixed_source.n==row.n)]
        count=sum(int(k)<=row.k for k in group.k_close)
        require(f'fixed_k_count_{row.method}_{row.k}_{row.n}',
                len(group)==row.samples and row.samples==(500 if row.n==32768 else 1500) and
                not group.seed.duplicated().any() and (group.largest_gap<1/3).all() and
                count==row.closed and abs(row.closure_percent-100*count/len(group))<1e-11)
    from rank_sweep import verify_rank_sweep
    candidate=verify_rank_sweep(root)
    require('rank_sweep_audit',candidate['status']=='passed',candidate)
    arc=json.loads((data/'arc_support_verification.json').read_text())
    require('arc_support_checker_hash',arc['checker_sha256']==hashlib.sha256((root/'code/verify_arc_support.py').read_bytes()).hexdigest())
    require('arc_support_small_graph_scope',arc['status']=='passed' and arc['small_graphs']['samples']==80 and arc['small_graphs']['supported_budgets']==3796 and arc['small_graphs']['empty_supported']>0 and arc['small_graphs']['small_gap_budgets']>0 and arc['small_graphs']['odd_field_checks']>0,arc['small_graphs'])
    require('arc_support_ellipse_closure',arc['ellipses']['coordinate_sets']==60 and arc['ellipses']['arc_supported']==180 and arc['ellipses']['saved_single_bar_diagrams']==180,arc['ellipses'])
    require('arc_support_six_point_witness',arc['witness']['beta1_at_point33_over_F2_F3_F5']==2 and len(arc['witness']['extra_edge_bars'])==2)
    from verify_studies import verify_studies
    studies=verify_studies(root)
    require('noise_and_ellipse_studies',studies['status']=='passed',studies)
    from audit_rank_scaling import audit_rank_scaling
    rank_audit=audit_rank_scaling(root)
    require('rank_scaling_raw_rows',rank_audit['status']=='passed',rank_audit)
    from verify_comparisons import verify_comparisons
    comparisons=verify_comparisons(root)
    require('study_comparison_evidence',comparisons['status']=='passed',comparisons)
    require('study_comparison_report_replay',json.loads((data/'comparisons/verification.json').read_text())==comparisons)
    verify_curve_order_records(root,require)
    verify_payload(root,require)
    summary_names={'release_structure','clean_native_panel','saved_sample_scope','eleven_event_rows','180_first_closed_series','closure_table_reports_all_errors','own_closure_not_global_split','grid_scope','120_paired_extra_k_rows','miss_means_unclosed_integer_layer','miss_counts_38_and_40','six_extra_k_groups','certificate_counts','CkNN_certificate_only','useful_diagram_scope','saved_diagram_transitions','matched_noise_scope','no_long_first_connected_bars','ellipse_scope','ellipse_constants','ellipse_loop_interval','ellipse_saved_diagram_inventory','manifest_file_inventory','manifest_exclusion_policy','code_syntax','permuted_scale_control','control_grid_matches','noise_closure_comparison','noise_closure_event_positions','noise_closure_diagram_inventory','noise_closure_summary'}
    summary={name:item for name,item in checks.items() if name in summary_names or 'detail' in item}
    return {'status':'passed','requested_scope':'complete_saved_data','measurement_generation':'not_run','additional_scientific_checks':'passed','manuscript_comparison':{'status':'passed' if paper is not None else 'not_run','location':'external' if paper is not None else None},'scope':'Table 3, Figure 3, the small examples and all additional studies; the varying-neck records and 240 one-edge prefix replays; the additional clean, noise and ellipse diagrams; seed and selector replay; finite witnesses; raw rank-sweep counts, Wilson intervals and study comparisons. No new samples or complete persistence-study regeneration.','checks_passed':len(checks)+studies['checks_passed']+rank_audit['checks_passed']+comparisons['checks_passed'],'checks':summary}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--report',type=Path)
    ap.add_argument('--manuscript-dir',type=Path,help='Directory containing the external combined.tex and combined.pdf; defaults to the parent if present.')
    ap.add_argument('--data-only',action='store_true',help='Check the numerical release without comparing a manuscript.')
    ap.add_argument('--paper-only',action='store_true',help='Only the checks behind the paper, plus every file hash; the additional studies are reported as not run.')
    args=ap.parse_args()
    if args.paper_only and args.report is None:
        ap.error('With --paper-only, pass --report with a path outside the release; the saved report covers the complete check.')
    if args.paper_only and args.report is not None and args.release.resolve() in args.report.resolve().parents:
        ap.error('Write --paper-only reports outside the release folder.')
    if args.data_only and args.manuscript_dir:ap.error('Choose --data-only or --manuscript-dir.')
    report=verify(args.release.resolve(),args.manuscript_dir,args.data_only,args.paper_only)
    (args.report or args.release/('results/paper_verification_report.json' if args.paper_only else 'results/verification_report.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(f"Release verification passed: {report['checks_passed']} checks.")
