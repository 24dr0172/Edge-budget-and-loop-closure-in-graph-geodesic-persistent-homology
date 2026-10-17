"""Check the displayed Phase 3 comparisons against the saved scientific summaries."""
from pathlib import Path
import csv,re

def verify_manuscript_tables(root,require):
    root=Path(root);text=(root.parent/'combined.tex').read_text()
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

if __name__=='__main__':
    count=0
    def require(name,condition):
        global count
        assert condition,name
        count+=1
    verify_manuscript_tables(Path(__file__).resolve().parents[1],require)
    print(f'{count} manuscript table checks passed.')
