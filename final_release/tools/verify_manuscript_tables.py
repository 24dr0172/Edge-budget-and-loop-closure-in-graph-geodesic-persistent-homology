"""Compare the focused manuscript with raw records, certificates and release files.

This expands all active inputs. Additional studies are checked separately by the
complete numerical verifier; they impose no extra manuscript sections.
"""
from pathlib import Path
import argparse,csv,json,re,sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from external_manuscript import locate_manuscript,manuscript_text,manuscript_inputs
from record_current import result_numbers

FIGURES={'fig:closure-schematic':'closure_schematic.png','fig:arc-support':'arc_support_fan.png','fig:retention-extra':'retention_extra_layers.png'}
TABLES={'tab:ingredients':None,'tab:notation':None,'tab:rank-scaling':'retention_compact.csv'}

def flat(t):return re.sub(r'\s+',' ',t.replace('{,}',','))
def between(t,a,b):
    i=t.index(a);return t[i:t.index(b,i)]
def rows(path):
    with path.open() as f:return list(csv.DictReader(f))

def verify_manuscript_tables(root,require,manuscript_dir=None):
    root=Path(root);paper=locate_manuscript(root,manuscript_dir)
    if paper is None:raise FileNotFoundError('Supply --manuscript-dir.')
    text=manuscript_text(paper).replace('{,}',',')
    require('focused_figures',re.findall(r'\\label\{(fig:[^}]+)\}',text)==list(FIGURES))
    require('focused_tables',re.findall(r'\\label\{(tab:[^}]+)\}',text)==list(TABLES))
    numbers=result_numbers(text)
    for kind,labels in {'theorem':['thm:loop-transition','thm:arc-supported','thm:circle-parameters'],
                        'corollary':['cor:cost-ceiling','cor:closure-optimal'],
                        'proposition':['prop:radius-enlargement','prop:locality']}.items():
        require('focused_result_numbers_'+kind,all(numbers[kind].get(v)==i for i,v in enumerate(labels,1)) and text.count(r'\begin{'+kind+'}')==len(labels))
    require('focused_path_lemma',numbers['lemma']['lem:edge-path']==1 and text.count(r'\begin{lemma}')==8)
    require('focused_notation_placement',text.index(r'\label{tab:notation}')<text.index(r'\section{Before ring closure}'))
    require('focused_scope',all(v not in text for v in [r'\label{sec:selection}',r'\label{sec:peanut-neck}',r'\label{app:noise-pilot}']) and
            re.search(r'supplement|online resource',text,re.I) is None)
    require('focused_appendices',re.findall(r'\\label\{(app:[^}]+)\}',text)==['app:curve-bounds','app:curve-numerics','app:ties','app:seven-point'])
    require('focused_terminology','ring completion' not in text.lower())
    require('focused_declarations',r'\subsection*{Funding} NA' in flat(text) and r'\subsection*{Competing interests}' in text and 'early drafts of some proofs' in text)
    abstract=flat(between(text,r'\begin{abstract}',r'\end{abstract}'))
    require('focused_abstract_scope','without noise' in abstract and 'uniform' in abstract and 'supported edges' in abstract)
    require('focused_checked_prior_work',all(v in text for v in ['metricgluings2020','boissonnatreach2019','federer1959','laurentseminaroti2017']))
    require('focused_inputs',manuscript_inputs(paper)==[])
    require('focused_no_external_text',re.search(r'\\(?:input|include|bibliography|addbibresource)\s*\{',text) is None)
    graphics=re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',text)
    require('focused_graphics',graphics==['figures/'+name for name in FIGURES.values()])
    require('focused_figure_directory',sorted(p.name for p in (root/'manuscript_figures').iterdir())==sorted(['README.md',*FIGURES.values()]))
    for name in FIGURES.values():
        require('focused_asset_'+name,(paper/'figures'/name).read_bytes()==(root/'manuscript_figures'/name).read_bytes())
    require('focused_release_no_latex',not list(root.rglob('*.tex')) and not list(root.rglob('*.bib')))
    require('focused_release_no_paper',not any(p.name in ('combined.pdf','main.pdf','combined.tex','main.tex','overleaf.zip','manuscript-png.zip') for p in root.rglob('*')))
    availability=between(text,r'\subsection*{Data and code}',r'\subsection*{Funding}')
    require('focused_path_placement',r'\path{' not in text.replace(availability,''))
    for value in re.findall(r'\\path\{([^}]+)\}',text):require('focused_release_path_'+value,(root/value).exists())

    from verify_retention_cost_witness import practical_counts,wilson
    counts=practical_counts(root,require)
    table=between(text,r'\label{tab:rank-scaling}',r'\end{table}')
    names={'Uniform':'uniform','Cosine':'strong_cosine','Trough':'localized_trough','OR-$k$NN':'OR_kNN','Mutual-$k$NN':'mutual_kNN'}
    displayed={}
    for line in table.splitlines():
        cell=[v.strip() for v in line.rstrip('\\ ').split('&')]
        if len(cell)==6 and cell[0] in names and cell[1] in names:
            for n,value in zip((128,512,4096,32768),cell[2:]):displayed[names[cell[0]],names[cell[1]],n]=value
    require('focused_table_inventory',set(displayed)==set(counts))
    exported={}
    for row in rows(root/'results/retention_table/retention_compact.csv'):
        for n in (128,512,4096,32768):exported[row['profile'],row['method'],n]=row[str(n)]
    require('focused_export_inventory',set(exported)==set(counts))
    for key,value in counts.items():
        require('focused_retention_'+str(key),displayed[key]==exported[key]==f"{100*value['k_conn']/value['samples']:.1f}")
    for size,limit in [(1500,2.5),(500,4.1)]:
        half=[50*(wilson(v['k_conn'],size)[1]-wilson(v['k_conn'],size)[0]) for v in counts.values() if v['samples']==size]
        require('focused_wilson_'+str(size),max(half)<limit and f'${limit:.1f}$' in text)
    practical=json.loads((root/'results/practical_claims_verification.json').read_text())
    body=flat(between(text,r'\label{sec:practical-checks}',r'\subsection{Diagram accuracy}'))
    for v in practical['fixed_k']:
        if v['method']=='mutual_kNN' and v['k']==10 and v['n'] in (512,4096):
            require('focused_fixed_k_'+str(v['n']),f"{v['closed']:,}/{v['samples']:,}" in body and f"{v['closure_percent']:.1f}"+r'\%' in body)
    for v in practical['bridge_free_n512']:
        if (v['profile'],v['method']) in [('uniform','OR_kNN'),('localized_trough','mutual_kNN')]:
            require('focused_bridge_'+v['profile'],f"{v['closed']:,}/{v['samples']:,}" in body and f"{v['closure_percent']:.1f}"+r'\%' in body)
    extra=rows(root/'results/retention_extra_layers/retention.csv')
    for profile,method,j in [('uniform','OR_kNN',3),('uniform','OR_kNN',5),('localized_trough','mutual_kNN',3),('localized_trough','mutual_kNN',5)]:
        v=next(v for v in extra if int(v['n'])==512 and v['profile']==profile and v['method']==method and int(v['j'])==j)
        require('focused_extra_layers_'+str((profile,method,j)),f"{100*float(v['rate']):.1f}"+r'\%' in body)
    model=between(text,r'\label{sec:heuristic}',r'\subsection{Required inputs')
    for v in practical['independent_key_model']:
        if v['n']==512:require('focused_model_'+v['method'],f"{v['probability']:.3f}" in model)
    require('focused_explicit_heuristic',r'\label{eq:poisson-retention}' in model and 'heuristic rather than a retention law' in flat(model))
    observed=between(text,r'\label{sec:diagram-experiment}',r'\section{Discussion}')
    for field in ['maximum','empty_cost_minimum','empty_cost_maximum']:
        require('focused_error_'+field,f"{practical['diagram_error'][field]:.3f}" in observed)
    require('focused_grid_scope','stored budget grid' in flat(observed))

    natural=json.loads((root/'results/natural_instance_verification.json').read_text())
    example=between(text,r'\label{ex:natural-ellipse}',r'\end{example}')
    require('focused_ellipse_margins',natural['status']=='passed' and natural['margins']['smallest_structural_margin']>0.03 and
            natural['margins']['smallest_enclosure_margin']>2e-5 and 'rational interval arithmetic' in flat(example))
    for key,lo,hi in [('04',1.7807,1.7808),('14',2.0018,2.0019),('12',1.9698,1.9699)]:
        v=natural['chords'][key]
        require('focused_ellipse_bound_'+key,lo<float(v['lower'])<=float(v['upper'])<hi and str(lo) in example and str(hi) in example)
    cert=json.loads((root/'results/small_examples_verification.json').read_text())['shortcut_decimal_certificate']
    require('focused_shortcut_certificate',cert['status']=='passed' and cert['lower_bound_margin']>1e-5 and 'Rational sine enclosures' in text)
    appendix=between(text,r'\label{app:seven-point}',r'\begingroup')
    witness=json.loads((root/'results/useful_diagram/rational_counterexample.json').read_text())['methods']['R2_global']
    before,after=witness['before'],witness['after'];diam=before['diameter']
    values=[f'{diam:.3f}',f"{before['B']:.6f}",f"{before['death']:.6f}",f"{after['death']:.6f}",f"{before['e_PH']:.7f}",f"{after['e_PH']:.7f}"]
    values += [f"{max(before['true_birth']-v['birth'],abs(v['true_death']-v['death']))/diam:.6f}" for v in (before,after)]
    require('focused_seven_point',all(v in appendix for v in values) and '$(0,0.08,0.15,0.161,0.29,0.60,0.66)$' in appendix and '[0.34,0.37)' in appendix and before['M_close']==13)
    require('focused_radius_constant',r'\Gamma(\varepsilon+2)n^{-\varepsilon}' in text and r'\label{eq:gap-ratio}' in text)
    require('focused_reach',r'\nu\ge2\tau>0' in text and r'c_{\rm ant}\ge2\tau' in text and '0.1032' in text and '0.1919' in text)
    protocol=flat(between(text,r'\label{sec:experimental-setup}',r'\subsection{Connectivity and ring closure}'))
    config=json.loads((root/'code/config.json').read_text())
    require('focused_tolerances',config['arithmetic']['clean_normalized_error_tolerance']==5e-8 and config['arithmetic']['support_tolerance']==1e-10 and
            r'$5\times10^{-8}$' in protocol and '$10^{-10}$' in protocol and 'do not impose a lifetime cutoff' in protocol)
    cutoff=json.loads((root/'results/cutoff_validation.json').read_text())
    require('focused_cutoff',cutoff['status']=='passed' and cutoff['random_positive_weight_graphs']==150 and cutoff['n512_saved_unrestricted_diagrams']==126 and
            '$150$' in protocol and '$126$' in protocol and 'adaptive persistence cutoff' in protocol)
    clouds=rows(root/'results/clean_closure/clouds.csv')
    require('focused_cloud_scope',len(clouds)==60 and all(int(v['n'])==512 for v in clouds) and f"{max(float(v['g_max']) for v in clouds):.3f}" in protocol)
    locality=json.loads((root/'results/noncircular_locality/locality_check_summary.json').read_text())
    local=flat(between(text,r'\label{app:curve-numerics}',r'\section{Ties and parameter layers}'))
    for kind in ('ellipse','peanut'):require('focused_locality_'+kind,f"{locality['r_loc'][kind]['r_loc']:.4f}" in local)
    require('focused_estimate_scope','grid evaluation' in local and 'grid estimates' in text and 'derivative zeros without a sign change' in flat(text))

    readme=(root/'README.md').read_text();citation=(root/'CITATION.cff').read_text()
    title=' '.join(re.search(r'\\title\{([\s\S]*?)\}\s*\\author',text)[1].replace('\\\\',' ').split())
    require('focused_release_title',readme.splitlines()[0]=='# '+title and citation.count(f'title: "{title}"')==2)
    for name in ['Shivam Singh','Annesha Sen','Tiwari']:require('focused_author_'+name,name in text and name in readme)
    for kind,assets in [('Table',TABLES),('Figure',FIGURES)]:
        mapped=[v for v in readme.splitlines() if re.match(r'\| '+kind+r' \d+ ',v)]
        require('focused_readme_'+kind,len(mapped)==len(assets))
        for i,(label,asset) in enumerate(assets.items(),1):
            line=next((v for v in mapped if v.startswith(f'| {kind} {i} ')),'')
            require('focused_map_'+label,bool(line) and (asset in line if asset else 'written directly in the manuscript' in line or 'drawn inline' in line))
    disclosure=flat(between(text,r'\subsection*{Use of generative-AI tools}',r'\subsection*{Author contributions}')).split('}',1)[1].strip()
    require('focused_ai_binding',disclosure in flat(readme))
    require('focused_no_supplement_files',not any(re.search('supplement',v.name,re.I) for v in root.rglob('*') if v.is_file()))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--manuscript-dir',type=Path)
    args=parser.parse_args();checks=[]
    def require(name,condition,detail=None):
        if not condition:raise AssertionError((name,detail))
        checks.append(name)
    verify_manuscript_tables(args.release.resolve(),require,args.manuscript_dir)
    print(f'Focused manuscript comparison passed: {len(checks)} checks.')
