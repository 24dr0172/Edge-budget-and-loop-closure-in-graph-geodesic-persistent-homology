"""Compare the manuscript with its saved data and figure assets.

The analytical checks verify source structure and displayed evidence, not the
mathematical proofs. The source and PDF identities are checked by verify_release.
"""
from pathlib import Path
import argparse
import csv
import json
import re
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from external_manuscript import locate_manuscript, manuscript_text, manuscript_inputs
from record_current import result_numbers

FIGURE_LABELS = ['fig:closure-schematic','fig:arc-support-fan','fig:rate','fig:peanut-comparison']
GRAPHICS = ['arc_support_fan.png','post_closure_rate.png','peanut_comparison.png']
TABLES = ['tab:rank-scaling','tab:radial-noise','tab:peanut-summary','tab:smooth-enlargement']
RESULTS = {
    'theorem': ['thm:loop-transition','thm:arc-supported','thm:circle-parameters','thm:ties','thm:smooth-radius'],
    'corollary': ['cor:geometric-survival','cor:first-nonempty','cor:closure-optimal','cor:rate','cor:certificate'],
    'proposition': ['prop:locality'],
    'lemma': ['lem:edge-path','lem:robinson','lem:circle-bar','lem:edgewise-death',
              'lem:rank-containment','lem:circle-orders','lem:local-rank-keys','lem:chord-arc','lem:cost-ceiling'],
    'remark': ['rem:polygon','rem:unsupported-shortcut'],
    'example': ['ex:natural-ellipse'],
}

def flat(text):
    return re.sub(r'\s+', ' ', text.replace('{,}', ','))

def section(text, start, end):
    i=text.index(start)
    return text[i:text.index(end, i)]

def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))

def verify_manuscript_tables(root, require, manuscript_dir=None):
    root=Path(root)
    paper=locate_manuscript(root,manuscript_dir)
    if paper is None:
        raise FileNotFoundError('Supply --manuscript-dir.')
    text=manuscript_text(paper).replace('{,}',',')
    numbers=result_numbers(text)
    require('focused_figure_labels', re.findall(r'\\label\{(fig:[^}]+)\}', text)==FIGURE_LABELS)
    require('focused_table_labels', re.findall(r'\\label\{(tab:[^}]+)\}', text)==TABLES)
    for kind, labels in RESULTS.items():
        require('focused_result_numbers_'+kind,
                text.count(r'\begin{'+kind+'}')==len(labels) and
                all(numbers[kind].get(label)==i for i,label in enumerate(labels,1)))
    headings=['Smooth-curve locality','Ties in edge orderings','Proof of the tied-rank theorem',
              'A seven-point example outside the gap condition','Comparisons with existing graph frameworks',
              'Finite-size check of radius enlargement','Protocols and numerical tolerances']
    appendix=text[text.index(r'\appendix'):]
    require('focused_appendices', re.findall(r'\\section\{([^}]+)\}',appendix)==headings)
    require('focused_intrinsic_proof_placement',
            appendix.index(r'\label{app:ties}')<
            appendix.index(r'\begin{proof}[Proof of Lemma~\ref{lem:circle-bar}]')<
            appendix.index(r'\label{app:tied-rank-proof}'))
    require('focused_no_external_tex_inputs',manuscript_inputs(paper)==[] and
            not re.search(r'\\(?:input|include|bibliography|addbibresource)\s*\{',text))
    graphics=re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',text)
    require('focused_graphics',graphics==['figures/'+x for x in GRAPHICS])
    require('focused_figure1_inline',r'\begin{picture}' in section(text,r'\begin{figure}',r'\end{figure}'))
    require('focused_figure_inventory',sorted(p.name for p in (root/'manuscript_figures').glob('*.png'))==sorted(GRAPHICS))
    for name in GRAPHICS:
        require('focused_asset_'+name,(paper/'figures'/name).read_bytes()==(root/'manuscript_figures'/name).read_bytes())
    require('focused_data_only_release',not list(root.rglob('*.tex')) and not list(root.rglob('*.bib')) and
            not any(p.name in ('combined.pdf','main.pdf','combined.tex','main.tex') for p in root.rglob('*')))
    require('focused_scope',r'\label{thm:general-first-nonempty}' not in text and
            r'\label{lem:winding-survival}' not in text and r'\label{fig:retention-extra}' not in text and
            r'\label{cor:first-nonempty}' in text and r'\label{thm:smooth-radius}' in text)
    abstract=flat(section(text,r'\begin{abstract}',r'\end{abstract}'))
    require('focused_abstract_scope',all(x in abstract for x in [
        'noise-free', 'consecutive edges', 'nearest-neighbor rules',
        'local counts', 'distances along the curve', 'stated sampling conditions',
        'On peanut curves, graph distances at cycle closure',
        r'can be connected and still have an empty shortest-path $H_1$ persistence diagram',
        'the sample locations and their cyclic order']))
    abstract_body=section(text,r'\begin{abstract}',r'\end{abstract}').split('}',1)[1]
    require('focused_abstract_nonvanishing', all(x in abstract for x in [
        'uniform circle samples', 'smallest connecting neighbor count',
        'is bounded away from zero']))
    require('focused_smooth_abstract_connected', all(x in abstract for x in [
        'sampled independently and uniformly in arc length', 'fixed smooth closed curve',
        'doubling the connecting radius leads to a miss probability asymptotic to two over the sample size']))
    abstract_words=len(re.sub(r'\$[^$]*\$', 'MATH', abstract_body).split())
    require('focused_abstract_plain_language',not re.search(
        r'short-hops-first|arc support|\bring\b|edge budget|winding triple',abstract_body,re.I))
    require('focused_abstract_no_citations',not re.search(r'\\cite\w*',abstract_body))
    require('focused_abstract_math_notation',
            all(symbol == 'H_1' for symbol in re.findall(r'\$([^$]*)\$',abstract_body)))
    require('focused_abstract_length',abstract_words<=180,{'words_with_each_math_expression_as_one':abstract_words})
    citations={key.strip() for group in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}',text) for key in group.split(',')}
    bibliography=set(re.findall(r'\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}',text))
    require('focused_bibliography_cited',citations==bibliography,{'cited':len(citations),'listed':len(bibliography)})
    introduction=flat(section(text,r'\section{Introduction}',r'\section{Related work}'))
    require('focused_intro_claims',all(x in introduction for x in [
        'connectivity precedes cycle closure', 'persistence in positive degrees is empty before closure',
        'probability is bounded away from zero as the sample size increases',
        r'$67.1\%$ of the uniform samples in Table',
        r'$\Gamma(\alpha+1)n^{1-\alpha}$', '$2/n$ for doubling',
        'edge to minimize diagram error']))
    require('focused_tie_proof_in_appendix',text.index(r'\appendix')<text.index(r'\label{app:tied-rank-proof}')<
            text.index(r'\begin{proof}[Proof of Theorem~\ref{thm:ties}]'))
    protocols=section(text,r'\label{app:protocols}',r'\begin{thebibliography}')
    require('focused_protocols_in_appendix',text.index(r'\appendix')<text.index(r'\label{app:protocols}') and
            all(r'\path{'+f+'}' in protocols for f in ['code/README.md','code/config.json']))
    protocol_guide=(root/'code/README.md').read_text()
    require('focused_tolerance_text',all(x in protocol_guide for x in ['1e-6','1e-10','5e-8','5e-12']))
    for name in re.findall(r'\\path\{([^}]+)\}',text):
        require('focused_release_path_'+name,(root/name).exists())

    # Every cell in the circle-retention table must agree with both its compact export and raw counts.
    from verify_retention_cost_witness import practical_counts, wilson
    counts=practical_counts(root,require)
    display={}
    profiles={'Uniform':'uniform','Cosine':'strong_cosine','Trough':'localized_trough'}
    methods={'Symmetric $k$NN':'OR_kNN','Mutual $k$NN':'mutual_kNN'}
    table=section(text,r'\label{tab:rank-scaling}',r'\end{table}')
    for line in table.splitlines():
        cell=[x.strip() for x in line.rstrip('\\ ').split('&')]
        if len(cell)==6 and cell[0] in profiles and cell[1] in methods:
            for n,value in zip((128,512,4096,32768),cell[2:]):
                display[profiles[cell[0]],methods[cell[1]],n]=value
    compact={}
    for row in rows(root/'results/retention_table/retention_compact.csv'):
        for n in (128,512,4096,32768):compact[row['profile'],row['method'],n]=row[str(n)]
    require('focused_retention_inventory',set(display)==set(compact)==set(counts) and len(display)==24)
    for key, value in counts.items():
        require('focused_retention_'+str(key),display[key]==compact[key]==f"{100*value['k_conn']/value['samples']:.1f}")
    for size, limit in ((1500,2.5),(500,4.1)):
        half=[50*(wilson(v['k_conn'],size)[1]-wilson(v['k_conn'],size)[0]) for v in counts.values() if v['samples']==size]
        require('focused_wilson_'+str(size),bool(half) and max(half)<limit)
    # Evaluate the approximation reported near Table 1.
    from check_tied_ranks import tie_approximations
    tie_statement=section(text,r'\label{thm:ties}',r'\end{theorem}')
    tie_comparison=flat(section(text,r'\label{thm:ties}','A separate clean-circle study'))
    require('focused_tie_asymptotic_remainder',r'e^{-n\pi_m}+o(1)' in tie_statement and
            'The $o(1)$ remainder is unquantified.' in tie_comparison)
    for rule, approximations in tie_approximations((128,512,4096,32768)).items():
        values=[f'{100*rate:.1f}'+r'\%' for rate in approximations]
        display_values=', '.join('$'+value+'$' for value in values[:-1])+', and $'+values[-1]+'$'
        require('focused_tie_table_approximation_'+rule,display_values in tie_comparison)
    from radial_noise_knn import AMPLITUDES, METHODS, N, SEEDS
    noise=rows(root/'results/radial_noise_knn.csv')
    expected_keys={(seed,delta,method) for seed in SEEDS for delta in AMPLITUDES for method in METHODS}
    actual_keys={(int(row['seed']),float(row['delta']),row['method']) for row in noise}
    require('focused_radial_pairing',len(noise)==len(expected_keys) and actual_keys==expected_keys)
    noise_table=section(text,r'\label{tab:radial-noise}',r'\end{table}')
    for row in noise:
        ratio=float(row['longest_ratio'])
        require('focused_radial_row_'+str((row['seed'],row['delta'],row['method'])),
                ratio>=0 and 1<=int(row['k'])<N and row['all_consecutive'] in ('0','1')
                and int(row['empty'])==(ratio==0)
                and int(row['at_least_quarter'])==(ratio>=.25)
                and int(row['at_least_half'])==(ratio>=.5)
                and (float(row['delta'])>0 or int(row['all_consecutive'])==1-int(row['empty'])))
    noise_display={}
    for line in noise_table.splitlines():
        cells=[cell.strip().strip('$') for cell in line.rstrip('\\ ').split('&')]
        if len(cells)==7 and all(cell.isdecimal() for cell in cells[1:]):
            noise_display[float(cells[0])]=list(map(int,cells[1:]))
    require('focused_radial_table_inventory',set(noise_display)==set(AMPLITUDES))
    require('focused_radial_lifetime_caption',
            'Using a quarter of the lifetime gives the same counts for both rules at every amplitude.'
            in flat(section(text,r'\caption{Counts at the smallest connecting',r'\label{tab:radial-noise}')))
    for delta in AMPLITUDES:
        cells=[]
        for method in METHODS:
            group=[r for r in noise if float(r['delta'])==delta and r['method']==method]
            require('focused_radial_equal_cutoffs_'+str((delta,method)),
                    sum(int(r['at_least_quarter']) for r in group)==
                    sum(int(r['at_least_half']) for r in group))
            cells.extend(sum(int(r[name]) for r in group)
                         for name in ('empty','at_least_half','all_consecutive'))
        require('focused_radial_table_'+str(delta),noise_display.get(delta)==cells)

    # Peanut closure table: four noncircular waists, preserving per-order counts and pooled ranges.
    from build_peanut_comparison import table_rows
    peanut=json.loads((root/'results/peanut_neck/closure_error_bound.json').read_text())
    table=section(text,r'\label{tab:peanut-summary}',r'\end{table}')
    require('focused_peanut_inventory',peanut['status']=='passed' and len(peanut['groups'])==15)
    display_rows=table_rows(peanut)
    require('focused_peanut_display_scope',len(display_rows)==4 and
            sum('&' in line for line in table.splitlines())==5)
    raw_closure=rows(root/'results/peanut_neck/trials.csv')
    for alpha,row in zip((.4,.6,.8,.9),display_rows):
        require('focused_peanut_'+str(alpha),row in table)
        raw=[r for r in raw_closure if float(r['alpha'])==alpha]
        expected_counts=[];expected_single=[]
        for method in ('R2_global','OR_kNN','mutual_kNN'):
            group=[r for r in raw if r['method']==method]
            require('focused_peanut_raw_scope_'+str((alpha,method)),len(group)==20)
            expected_counts.append(sum(r['closure_supported']=='True' for r in group))
            expected_single.append(sum(int(r['closure_bars'])==1 for r in group))
        cells=[cell.strip() for cell in row.rstrip('\\ ').split('&')]
        require('focused_peanut_raw_counts_'+str(alpha),
                cells[1]=='/'.join(map(str,expected_counts)) and
                cells[2]=='/'.join(map(str,expected_single)))
        lo=min(float(r['closure_error']) for r in raw)*1000
        hi=max(float(r['closure_error']) for r in raw)*1000
        require('focused_peanut_raw_range_'+str(alpha),cells[-1]==f'${lo:.3g}$--${hi:.3g}$')
    trial=rows(root/'results/first_nonempty/trials.csv')
    require('focused_peanut_trial_scope',len(trial)==300)
    for a in (0,.4,.6,.8):
        require('focused_peanut_first_'+str(a),all(r['at_closure']=='True' for r in trial if float(r['alpha'])==a))
    late=[r for r in trial if float(r['alpha'])==.9]
    expected={'R2_global':20,'OR_kNN':13,'mutual_kNN':18}
    for method,count in expected.items():
        group=[r for r in late if r['method']==method]
        require('focused_peanut_first_.9_'+method,len(group)==20 and sum(r['at_closure']=='True' for r in group)==count)
    early=[r for r in late if r['early']=='True']
    require('focused_peanut_early',len(early)==9 and all(r['short_hops_first']=='False' for r in early))
    radius=[r for r in trial if r['method']=='R2_global']
    closure_trials=rows(root/'results/peanut_neck/trials.csv')
    empty_matches=[r for r in closure_trials if float(r['alpha'])>=.6]
    require('focused_peanut_identical_euclidean_ranges',len(empty_matches)==180 and
            all(float(r['euclidean_error'])==float(r['P_star']) for r in empty_matches) and
            r'For each $a\ge0.6$, Euclidean error equals $P_\star$' in flat(text) and
            '$0.31361$ to $0.32476$' in flat(text))
    euclidean_by_cloud={}
    for r in closure_trials:
        key=(float(r['alpha']),r['seed'])
        count=int(r['euclidean_bars'])
        if key in euclidean_by_cloud:
            require('focused_euclidean_count_shared_'+str(key)+'_'+r['method'],euclidean_by_cloud[key]==count)
        euclidean_by_cloud[key]=count
    import numpy as np
    for (alpha,seed),count in euclidean_by_cloud.items():
        with np.load(root/'results/peanut_neck/checkpoints'/f'a{alpha:g}_seed{seed}.npz',allow_pickle=False) as cloud:
            bars=cloud['euclidean_bars']
            require('focused_euclidean_diagram_count_'+str((alpha,seed)),
                    len(bars)==count and np.all(bars[:,1]>bars[:,0]))
            if alpha==.4:
                reference=cloud['reference_bars'][0]
                lifetime=bars[:,1]-bars[:,0];i=int(np.argmax(lifetime))
                death_displacement=abs(float(bars[i,1])-float(reference[1]))
                latent=cloud['s'];difference=np.abs(latent[:,None]-latent[None,:])
                diameter=float(np.minimum(difference,1-difference).max())
                row=next(r for r in closure_trials if float(r['alpha'])==alpha and r['seed']==seed)
                require('focused_waist04_dominant_death_'+seed,
                        death_displacement>=abs(float(bars[i,0])-float(reference[0])) and
                        death_displacement>=float(np.max(np.delete(lifetime,i)/2,initial=0)) and
                        abs(death_displacement/diameter-float(row['euclidean_error']))<1e-10)
    require('focused_euclidean_discarded_reference',all(r['euclidean_reference_discarded']=='True'
            for r in closure_trials if float(r['alpha'])>=.6))
    expected_bars={0.:(1,1),.4:(2,3),.6:(2,2),.8:(2,2),.9:(2,2)}
    for alpha,expected_range in expected_bars.items():
        counts=[value for (a,seed),value in euclidean_by_cloud.items() if a==alpha]
        require('focused_euclidean_bar_range_'+str(alpha),len(counts)==20 and (min(counts),max(counts))==expected_range)
    require('focused_narrow_waist_support',all(r['closure_supported']=='False' for r in closure_trials
            if float(r['alpha'])==.9 and r['method'] in ('R2_global','OR_kNN')))
    require('focused_peanut_doubling',len(radius)==100 and all(r['double_ring_present']=='True' and int(r['double_bars'])>0 for r in radius))
    prose=section(text,r'\label{fig:peanut-comparison}',r'\section{Limitations and future work}')
    require('focused_peanut_reported_counts',all(v in flat(prose) for v in ['$20/20$', '$13/20$', '$18/20$', 'nine cases with a nonempty diagram before cycle closure', '$100$ point sets']))

    # Check the high-precision first-failure ranges and their reported margin.
    local_rows=rows(root/'results/noncircular_locality/locality_check_results.csv')
    local_summary=json.loads((root/'results/noncircular_locality/locality_check_summary.json').read_text())
    local_text=flat(section(text,r'\paragraph{A numerical check on two curves.}',r'\section{Ties in edge orderings}'))
    for curve,lo,hi,digits,tolerance in (('peanut',.1661840819,.1698344462,10,5.1e-11),
                                         ('ellipse',.1942495,.2064328,7,1e-7)):
        vals=[float(row[method+'_fail_chord']) for row in local_rows if row['curve']==curve
              for method in ('radius','OR','mutual')]
        require('focused_local_failure_range_'+curve,len(vals)==60 and abs(min(vals)-lo)<tolerance and
                abs(max(vals)-hi)<tolerance and f'${lo:.{digits}f}$--${hi:.{digits}f}$' in local_text)
    threshold=local_summary['r_loc']['peanut']['r_loc']
    margin=min(float(row[method+'_fail_chord']) for row in local_rows if row['curve']=='peanut'
               for method in ('radius','OR','mutual'))-threshold
    require('focused_local_peanut_margin',4.005e-7<margin<4.015e-7 and r'$4.01\times10^{-7}$' in local_text)

    # Radius enlargement table: circle observations, two curve proxies, spacing law and asymptotic rate.
    circle={(int(r['n']),round(1+float(r['epsilon']),2)):r for r in rows(root/'results/radius_enlargement/summary.csv')}
    from radius_enlargement_extra import verify as verify_extra_circle
    verify_extra_circle(root,require)
    extra_circle={(int(r['n']),float(r['alpha'])):r for r in rows(root/'results/radius_enlargement_extra/summary.csv')}
    require('focused_extra_circle_rows',set(extra_circle)=={(4096,2.),(32768,2.)})
    smooth={(int(r['n']),float(r['alpha']),r['curve']):r for r in rows(root/'results/smooth_curve_enlargement/summary.csv')}
    table=section(text,r'\label{tab:smooth-enlargement}',r'\end{table}')
    displayed={}
    for line in table.splitlines():
        cell=[x.strip().strip('$').replace(',','') for x in line.rstrip('\\ ').split('&')]
        if len(cell)==7:
            try:displayed[int(cell[1]),float(cell[0])]=[float(x) for x in cell[2:]]
            except ValueError:pass
    require('focused_radius_table_inventory',set(displayed)==set(circle) and len(displayed)==12)
    for (n,a), values in displayed.items():
        c=circle[n,a];e=smooth[n,a,'ellipse'];p=smooth[n,a,'peanut']
        circle_rate=extra_circle.get((n,a),c)['miss_rate']
        expected=[100*float(circle_rate),100*float(e['miss_rate']),100*float(p['miss_rate']),
                  100*float(c['spacing_probability']),100*float(c['asymptotic_probability'])]
        require('focused_radius_table_'+str((n,a)),all(abs(x-y)<.000051 for x,y in zip(values,expected)),(values,expected))
    require('focused_proxy_scope','chord-event proxy' in text and '$13,050$' in text)
    cert=json.loads((root/'results/rational_certificates/certificates.json').read_text())
    require('focused_exact_certificate',cert['six_point_example']['status']=='passed' and cert['threshold_root']['status']=='passed')
    natural=json.loads((root/'results/natural_instance_verification.json').read_text())
    require('focused_ellipse_certificate',natural['status']=='passed' and natural['margins']['smallest_structural_margin']>.03)
    seven=json.loads((root/'results/useful_diagram/rational_counterexample.json').read_text())['methods']['R2_global']
    example=section(text,r'\label{app:seven-point}',r'\section{Comparisons with existing graph frameworks}')
    before,after=seven['before'],seven['after']
    require('focused_seven_point',all(v in example for v in [f"{before['B']:.6f}",f"{before['death']:.6f}",
            f"{after['death']:.6f}",f"{before['e_PH']:.7f}",f"{after['e_PH']:.7f}"]))
    readme=(root/'README.md').read_text();citation=(root/'CITATION.cff').read_text()
    title=' '.join(re.search(r'\\title\{([\s\S]*?)\}\s*\\author',text)[1].replace('\\\\',' ').split())
    require('focused_title',readme.splitlines()[0]=='# '+title and citation.count(f'title: "{title}"')==2)
    disclosure=flat(section(text,r'\subsection*{Use of generative-AI tools}',r'\subsection*{Author contributions}')).split('}',1)[1].strip()
    disclosure=disclosure.replace(r'\LaTeX{}','LaTeX')
    require('focused_ai_binding',disclosure in flat(readme))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--manuscript-dir',type=Path)
    args=p.parse_args(); checks=[]
    def require(name,condition,detail=None):
        if not condition:raise AssertionError((name,detail))
        checks.append(name)
    verify_manuscript_tables(args.release.resolve(),require,args.manuscript_dir)
    print(f'Focused manuscript comparison passed: {len(checks)} checks.')
