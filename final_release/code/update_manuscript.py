"""Refresh numerical table bodies in combined.tex from released summaries."""
from pathlib import Path
import csv,math
F=Path(__file__).resolve().parents[1];P=F.parent;R=F/'results'
METHODS={'R2_global':'Radius','OR_kNN':'OR-$k$NN','mutual_kNN':'Mutual-$k$NN'}
PROFILES={'uniform':'Uniform','strong_cosine':'Cosine','localized_trough':'Trough'}
def rows(file):return list(csv.DictReader(file.open()))
def main():
    path=P/'combined.tex';text=path.read_text()
    def body(label,lines):
        nonlocal text
        anchor=text.index('\\label{'+label+'}')
        long_end=text.find(r'\end{longtable}',anchor)
        tab_end=text.find(r'\end{tabular}',anchor)
        if long_end!=-1 and (tab_end==-1 or long_end<tab_end):
            start=text.index(r'\endlastfoot',anchor)+len(r'\endlastfoot')
            end=long_end
        else:
            start=text.index(r'\midrule',anchor)+len(r'\midrule')
            end=text.index(r'\bottomrule',start)
        text=text[:start]+'\n'+'\n'.join(lines)+'\n'+text[end:]
    noise=rows(R/'noise/selector_summary.csv');lines=[]
    for model in ['radial','isotropic']:
        for method,name in METHODS.items():
            values=[]
            for selector in ['first_nonempty','longest_lifetime']:
                r=next(r for r in noise if r['model']==model and r['method']==method and float(r['delta'])==.01 and r['selector']==selector)
                error=(r'\shortstack[r]{$'+f"{float(r['median_e_PH']):.4f}"+
                       r'$\\$['+f"{float(r['q1_e_PH']):.4f}"+r',\,'+
                       f"{float(r['q3_e_PH']):.4f}"+r']$}')
                values += [str(int(r['before_closure']))+'/20',f"{float(r['median_M'])/(512*math.log(512)):.2f}",error]
            lines.append(model.capitalize()+' & '+name+' & '+' & '.join(values)+r' \\')
    body('tab:noise-selectors',lines)
    bridge=rows(R/'phase1/bridgeless/sweep_summary.csv');lines=[]
    for profile,pname in PROFILES.items():
        for method in ['OR_kNN','mutual_kNN']:
            r=next(r for r in bridge if r['profile']==profile and r['method']==method and int(r['n'])==512)
            values=[]
            for event in ['connection','bridgeless','plus3']:
                values.append('$'+f"{100*float(r[event+'_rate']):.1f}"+
                              r'\ ['+f"{100*float(r[event+'_wilson_lo']):.1f}, {100*float(r[event+'_wilson_hi']):.1f}"+']$')
            lines.append(pname+' & '+METHODS[method]+' & '+' & '.join(values)+r' \\')
    body('tab:extra-k-main',lines)
    def error(r):
        return (r'\shortstack[r]{$'+f"{float(r['median_e_PH']):.4f}"+
                r'$\\$['+f"{float(r['q1_e_PH']):.4f}"+r',\,'+
                f"{float(r['q3_e_PH']):.4f}"+r']$}')
    euclidean=rows(R/'phase1/euclidean/summary.csv')
    paired=rows(R/'phase1/euclidean/paired_summary.csv')
    oracle=rows(R/'phase1/grid_oracle_summary.csv');lines=[]
    for model in ['radial','isotropic']:
        baseline=next(r for r in euclidean if r['model']==model and float(r['delta'])==.01)
        for method,name in METHODS.items():
            r=next(r for r in oracle if r['model']==model and r['method']==method and float(r['delta'])==.01)
            paired_row=next(r for r in paired if r['model']==model and r['method']==method and float(r['delta'])==.01 and r['selector']=='longest_lifetime')
            values=[error(baseline),error(r),str(int(paired_row['selected_better']))+'/20']
            lines.append(model.capitalize()+' & '+name+' & '+' & '.join(values)+r' \\')
    body('tab:noise-comparators',lines)
    amplitude=rows(R/'phase1/all_amplitudes.csv');lines=[]
    for model in ['radial','isotropic']:
        if lines:lines.append(r'\addlinespace')
        for delta in [0,.0005,.001,.002,.005,.01]:
            values=[]
            for method in METHODS:
                r=next(r for r in amplitude if r['model']==model and r['method']==method and float(r['delta'])==delta)
                values.extend(f"{float(r[selector+'_median_e_PH']):.4f}" for selector in ['first_nonempty','longest_lifetime'])
            lines.append(model.capitalize()+' & '+format(delta,'.4f').rstrip('0').rstrip('.')+' & '+' & '.join(values)+r' \\')
    clean_baselines=[r for r in euclidean if float(r['delta'])==0]
    assert len(clean_baselines)==2 and abs(float(clean_baselines[0]['median_e_PH'])-float(clean_baselines[1]['median_e_PH']))<1e-12
    lines.extend([r'\addlinespace',r'\multicolumn{2}{l}{Euclidean, $\delta=0$} & \multicolumn{6}{c}{'+f"{float(clean_baselines[0]['median_e_PH']):.4f}"+r'} \\'])
    body('tab:noise-amplitudes',lines)
    screen=rows(R/'noise/screen_summary.csv');lines=[]
    for model in ['radial','isotropic']:
        for method,name in METHODS.items():
            values=[]
            for event in ['first_connected','before_closure']:
                r=next(r for r in screen if r['model']==model and r['method']==method and float(r['delta'])==.01 and r['event']==event)
                values += [str(int(r[k]))+'/20' for k in ['quarter','half','three_quarters']]
            lines.append(model.capitalize()+' & '+name+' & '+' & '.join(values)+r' \\')
    body('tab:noise',lines)
    ellipse=rows(R/'ellipse/summary.csv');lines=[]
    for profile,pname in PROFILES.items():
        for method,mname in METHODS.items():
            r=next(r for r in ellipse if r['profile']==profile and r['method']==method)
            median=f"{float(r['median_first_nonlocal']):,.1f}".removesuffix('.0')
            values=[median]+[str(int(r[k])) for k in ['cells','support_failing_cells','multiple_bar_cells']]
            lines.append(pname+' & '+mname+' & '+' & '.join(values)+r' \\')
    body('tab:ellipse-sweep',lines)
    scaling=rows(R/'extra_k/candidate_audit/summary.csv');lines=[]
    for profile,pname in PROFILES.items():
        if lines:lines.append(r'\addlinespace')
        for n in [128,512,4096,32768]:
            ntex=f'{n:,}' if n!=1000000 else '$10^6$'
            values=[]
            for method in ['OR_kNN','mutual_kNN']:
                r=next(r for r in scaling if r['profile']==profile and r['method']==method and int(r['n'])==n)
                for j in [0,3]:
                    value=f"{100*float(r[f'cdf_{j}']):.1f}"
                    values.append(value)
            lines.append(pname+' & '+ntex+' & '+' & '.join(values)+r' \\')
    body('tab:rank-scaling',lines)
    if text!=path.read_text():path.write_text(text)
    print('Checked/refreshed the seven numerical table bodies; prose and proofs are untouched.')
if __name__=='__main__':main()
