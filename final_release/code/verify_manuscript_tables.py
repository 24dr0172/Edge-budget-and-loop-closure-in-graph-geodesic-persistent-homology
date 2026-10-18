"""Check the current manuscript tables and Table 3 export against the saved scientific summaries."""
from pathlib import Path
import argparse,csv,re
from external_manuscript import locate_manuscript

def verify_manuscript_tables(root,require,manuscript_dir=None):
    root=Path(root);paper=locate_manuscript(root,manuscript_dir)
    if paper is None:raise FileNotFoundError('Supply --manuscript-dir to compare the external paper.')
    text=(paper/'combined.tex').read_text()
    def rows(name):
        with (root/'results'/name).open() as f:return list(csv.DictReader(f))
    def table(label,starts):
        start=text.index('\\label{'+label+'}')
        ends=[text.find(tok,start) for tok in (r'\end{tabular}', r'\end{longtable}')]
        block=text[start:min(e for e in ends if e!=-1)]
        return [line for line in block.splitlines() if line.startswith(starts)]
    def numbers(line):
        return [float(x) for x in re.findall(r'\d+(?:\.\d+)?',' & '.join(line.split(' & ')[2:]))]
    def rounded(value,places):return float(format(float(value),'.'+str(places)+'f'))
    profiles=['uniform','strong_cosine','localized_trough']
    methods=['R2_global','OR_kNN','mutual_kNN']
    bridge=rows('phase1/bridgeless/sweep_summary.csv')
    lines=table('tab:extra-k-main',('Uniform &','Cosine &','Trough &'))
    require('manuscript_bridge_table_scope',len(lines)==6)
    for line,(profile,method) in zip(lines,[(p,m) for p in profiles for m in methods[1:]]):
        row=next(r for r in bridge if int(r['n'])==512 and r['profile']==profile and r['method']==method)
        expected=[rounded(100*float(row[event+'_'+field]),1) for event in ['connection','bridgeless','plus3'] for field in ['rate','wilson_lo','wilson_hi']]
        require('manuscript_bridge_table_'+profile+'_'+method,numbers(line)==expected)
    exported=rows('extra_k/candidate_audit/main_text_summary.csv')
    ranks=rows('extra_k/candidate_audit/summary.csv')
    require('table3_export_scope',len(exported)==6 and
            [(r['profile'],r['method']) for r in exported]==[(p,m) for p in profiles for m in methods[1:]])
    for row in exported:
        rank=next(r for r in ranks if r['profile']==row['profile'] and r['method']==row['method'] and int(r['n'])==512)
        b=next(r for r in bridge if r['profile']==row['profile'] and r['method']==row['method'] and int(r['n'])==512)
        require('table3_export_sample_count_'+row['profile']+'_'+row['method'],
                int(row['n'])==512 and int(row['samples'])==int(rank['samples'])==int(b['samples']))
        for event, j in [('connection',0),('plus3',3),('bridgeless',None)]:
            for field, suffix in [('rate',''),('wilson_lo','_lo'),('wilson_hi','_hi')]:
                source=float(b[event+'_'+field]) if j is None else float(rank[f'cdf_{j}{suffix}'])
                require('table3_export_'+row['profile']+'_'+row['method']+'_'+event+'_'+field,
                        abs(float(row[event+'_'+field+'_percent'])-100*source)<1e-11)
    baseline=rows('phase1/euclidean/summary.csv')
    oracle=rows('phase1/grid_oracle_summary.csv')
    pairs=rows('phase1/euclidean/paired_summary.csv')
    lines=table('tab:noise-comparators',('Radial &','Isotropic &'))
    require('manuscript_comparator_table_scope',len(lines)==6)
    for line,(model,method) in zip(lines,[(m,k) for m in ['radial','isotropic'] for k in methods]):
        b=next(r for r in baseline if r['model']==model and float(r['delta'])==.01)
        o=next(r for r in oracle if r['model']==model and float(r['delta'])==.01 and r['method']==method)
        q=next(r for r in pairs if r['model']==model and float(r['delta'])==.01 and r['method']==method and r['selector']=='longest_lifetime')
        expected=[rounded(r[field],4) for r in [b,o] for field in ['median_e_PH','q1_e_PH','q3_e_PH']]+[float(q['selected_better']),float(q['samples'])]
        require('manuscript_comparator_table_'+model+'_'+method,numbers(line)==expected)
    clean_baselines=[r for r in baseline if float(r['delta'])==0]
    expected_baseline=(r'\multicolumn{2}{l}{Euclidean, $\delta=0$} & \multicolumn{6}{c}{'+
                       f"{float(clean_baselines[0]['median_e_PH']):.4f}"+r'} \\')
    require('manuscript_amplitude_clean_euclidean_baseline',len(clean_baselines)==2 and
            abs(float(clean_baselines[0]['median_e_PH'])-float(clean_baselines[1]['median_e_PH']))<1e-12 and
            expected_baseline in text)
    amplitudes=rows('phase1/all_amplitudes.csv')
    lines=table('tab:noise-amplitudes',('Radial &','Isotropic &'))
    require('manuscript_amplitude_table_scope',len(lines)==12)
    for line,(model,delta) in zip(lines,[(m,d) for m in ['radial','isotropic'] for d in [0,.0005,.001,.002,.005,.01]]):
        require('manuscript_amplitude_label_'+model+'_'+str(delta),float(line.split(' & ')[1])==delta)
        expected=[]
        for method in methods:
            r=next(r for r in amplitudes if r['model']==model and float(r['delta'])==delta and r['method']==method)
            expected.extend(rounded(r[rule+'_median_e_PH'],4) for rule in ['first_nonempty','longest_lifetime'])
        require('manuscript_amplitude_table_'+model+'_'+str(delta),numbers(line)==expected)
    amp=text.index('\\label{tab:noise-amplitudes}')
    require('manuscript_euclidean_zero_row','0.1154' in text[amp:text.index(r'\end{tabular}',amp)])
    # Compare the remaining numerical displays directly with measured rows.
    clean=rows('clean_closure/cells.csv');events=rows('clean_closure/knn_event_checks.csv')
    landmark_lines=table('tab:landmarks',('Radius &','OR-$k$NN &','Mutual-$k$NN &'))
    labels={'Radius':'R2_global','OR-$k$NN':'OR_kNN','Mutual-$k$NN':'mutual_kNN'}
    for line in landmark_lines:
        parts=line.split(' & ');method=labels[parts[0]];event=parts[1]
        for profile,cell in zip(profiles,parts[2:]):
            expected=list(map(int,re.findall(r'\d+',cell)))
            if event=='First connected / before closure':
                sets=[ [r for r in clean if r['profile']==profile and r['method']==method and r['budget_label']=='first_connected'] ]
                sets.append([r for r in (clean if method=='R2_global' else events) if r['profile']==profile and r['method']==method and r.get('budget_label',r.get('landmark'))==('pre_completion' if method=='R2_global' else 'pre_coverage_last')])
            elif event=='First connected integer $k$':
                sets=[[r for r in clean if r['profile']==profile and r['method']==method and r['budget_label']=='smallest_connected_k']]
            else:
                sets=[[r for r in (clean if method=='R2_global' else events) if r['profile']==profile and r['method']==method and r.get('budget_label',r.get('landmark'))==('completion' if method=='R2_global' else 'own_coverage')]]
            for index,group in enumerate(sets):
                got=[sum(float(r['est_h1_bars'])==0 for r in group),sum(float(r['est_h1_bars'])>0 for r in group)]
                require(f'manuscript_event_{method}_{profile}_{event}_{index}',len(group)==20 and got==expected)
    predicted=rows('extra_k/candidate_audit/iid_predictions.csv')
    for line in table('tab:race-model',('128 &','512 &','4,096 &','32,768 &')):
        n=int(line.split(' & ')[0].replace(',',''));expected=[]
        for method in methods[1:]:
            p=next(r for r in predicted if int(r['n'])==n and r['method']==method)
            q=next(r for r in ranks if int(r['n'])==n and r['method']==method and r['profile']=='uniform')
            expected.extend([rounded(100*float(p['cdf_0']),1),rounded(100*float(q['cdf_0']),1)])
        got=[float(x) for x in re.findall(r'\d+(?:\.\d+)?',' & '.join(line.split(' & ')[1:]))]
        require(f'manuscript_race_{n}',got==expected)
    selector_rows=rows('noise/selector_summary.csv')
    for line,(model,method) in zip(table('tab:noise-selectors',('Radial &','Isotropic &')),[(m,k) for m in ['radial','isotropic'] for k in methods]):
        expected=[]
        for rule in ['first_nonempty','longest_lifetime']:
            q=next(r for r in selector_rows if r['model']==model and r['method']==method and r['selector']==rule and float(r['delta'])==.01)
            from math import log
            expected.extend([float(q['before_closure']),float(q['samples']),rounded(float(q['median_M'])/(512*log(512)),2)])
            expected.extend(rounded(q[field],4) for field in ['median_e_PH','q1_e_PH','q3_e_PH'])
        require(f'manuscript_selectors_{model}_{method}',numbers(line)==expected)
    ellipses=rows('ellipse/summary.csv')
    for line,(profile,method) in zip(table('tab:ellipse-sweep',('Uniform &','Cosine &','Trough &')),[(p,m) for p in profiles for m in methods]):
        q=next(r for r in ellipses if r['profile']==profile and r['method']==method)
        expected=[float(q[k]) for k in ['median_first_nonlocal','cells','support_failing_cells','multiple_bar_cells']]
        require(f'manuscript_ellipse_{profile}_{method}',numbers(line.replace(',',''))==expected)
    screens=rows('noise/screen_summary.csv')
    for line,(model,method) in zip(table('tab:noise',('Radial &','Isotropic &')),[(m,k) for m in ['radial','isotropic'] for k in methods]):
        expected=[]
        for event in ['first_connected','before_closure']:
            q=next(r for r in screens if r['model']==model and r['method']==method and r['event']==event and float(r['delta'])==.01)
            for field in ['quarter','half','three_quarters']:expected.extend([float(q[field]),float(q['samples'])])
        require(f'manuscript_lifetime_screen_{model}_{method}',numbers(line)==expected)
    exact=rows('exact_closure_errors.csv');grid_closure=rows('clean_closure/first_closed_errors.csv')
    for line,(profile,method) in zip(table('tab:closed-errors',('Uniform &','Cosine &','Trough &')),[(p,m) for p in profiles for m in methods]):
        e=next(r for r in exact if r['profile']==profile and r['method']==method)
        g=next(r for r in grid_closure if r['profile']==profile and r['method']==method)
        expected=[rounded(e[k],5) for k in ['e_inf','e_2','p1']]+[rounded(g[k],5) for k in ['median_e_inf','median_e_2','median_e_PH']]
        require(f'manuscript_closed_errors_{profile}_{method}',numbers(line)==expected)
    for line in table('tab:rank-scaling',('Uniform &','Cosine &','Trough &')):
        parts=line.split(' & ');profile={'Uniform':'uniform','Cosine':'strong_cosine','Trough':'localized_trough'}[parts[0]];n=int(parts[1].replace(',',''));expected=[]
        for method in methods[1:]:
            q=next(r for r in ranks if r['profile']==profile and r['method']==method and int(r['n'])==n)
            expected.extend(rounded(100*float(q[k]),1) for k in ['cdf_0','cdf_3'])
        require(f'manuscript_rank_display_{profile}_{n}',numbers(line)==expected)
    maximum=max(float(r['birth_deviation']) for r in rows('useful_diagram/saved_prefixes.csv'))
    require('manuscript_birth_deviation_rounding',f'{maximum/1e-9:.2f}'+r'\times10^{-9}' in text)

if __name__=='__main__':
    count=0
    def require(name,condition):
        global count
        assert condition,name
        count+=1
    parser=argparse.ArgumentParser();parser.add_argument('--manuscript-dir',type=Path)
    verify_manuscript_tables(Path(__file__).resolve().parents[1],require,parser.parse_args().manuscript_dir)
    print(f'{count} manuscript table checks passed.')
