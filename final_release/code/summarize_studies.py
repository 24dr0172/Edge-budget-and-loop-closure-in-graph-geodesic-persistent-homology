"""Build study tables and figures from saved per-sample measurements."""
from pathlib import Path
import json,sys
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
PAPER=Path(__file__).resolve().parents[1];SUPPORT=PAPER
RESULTS=SUPPORT/'results';FIGURES=SUPPORT/'manuscript_figures'
sys.path.insert(0,str(PAPER/'code'))
import paper1_methods as m
METHODS=list(m.METHODS)
NAMES=dict(zip(METHODS,['Radius','OR-kNN','Mutual-kNN']))
COLORS=dict(zip(METHODS,['#2166ac','#d6604d','#4d9221']))
PROFILES={'uniform':'Uniform','strong_cosine':'Cosine','localized_trough':'Trough'}

def ellipse_summary():
    rows=pd.concat([pd.read_csv(p) for p in sorted((RESULTS/'ellipse/checkpoints').glob('*.csv'))],ignore_index=True)
    sample=pd.DataFrame([json.loads(p.read_text()) for p in sorted((RESULTS/'ellipse/checkpoints').glob('*.json'))])
    assert len(sample)==60
    rows=rows.merge(sample[['seed','profile','g','b_star']],on=['seed','profile'],validate='many_to_one')
    # Exact same-sample intrinsic diameter, not a continuum substitute.
    with np.load(PAPER/'results/ellipse_checks/coordinates.npz') as z:
        diam={(int(r.seed),r.profile):float(m.intrinsic(z[f'{r.seed}_{r.profile}_latent']).max()) for r in sample.itertuples()}
    rows['P_star']=[(r.b_star-r.g)/(2*diam[(int(r.seed),r.profile)]) for r in rows.itertuples()]
    series=[]
    for (seed,profile,method),g in rows.groupby(['seed','profile','method']):
        g=g.sort_values('M')
        event=g[g.first_nonlocal].iloc[0]
        series.append(dict(seed=seed,profile=profile,method=method,cells=len(g),
                           first_nonlocal=int(g.M_first_nonlocal.iloc[0]),
                           first_nonlocal_supported=bool(event.arc_supported),first_nonlocal_bars=int(event.bars),
                           support_failing_cells=int((~g.arc_supported).sum()),
                           multiple_bar_cells=int((g.bars>1).sum()),maximum_bars=int(g.bars.max()),
                           ordinary_decreases=int((g.e_PH.diff() < -1e-6).sum()),
                           capped_decreases=int((np.minimum(g.e_PH,g.P_star).diff() < -1e-6).sum())))
    series=pd.DataFrame(series)
    summary=series.groupby(['profile','method']).agg(samples=('seed','count'),cells=('cells','sum'),
        median_first_nonlocal=('first_nonlocal','median'),minimum_first_nonlocal=('first_nonlocal','min'),
        maximum_first_nonlocal=('first_nonlocal','max'),support_failing_cells=('support_failing_cells','sum'),
        first_nonlocal_failing_series=('first_nonlocal_supported',lambda v:int((~v).sum())),
        first_nonlocal_multiple_bar_series=('first_nonlocal_bars',lambda v:int((v>1).sum())),
        support_failing_series=('support_failing_cells',lambda v:int((v>0).sum())),
        multiple_bar_cells=('multiple_bar_cells','sum'),maximum_bars=('maximum_bars','max'),
        ordinary_decreasing_series=('ordinary_decreases',lambda v:int((v>0).sum())),
        capped_decreasing_series=('capped_decreases',lambda v:int((v>0).sum()))).reset_index()
    assert rows[rows.maximum_chord<.18].arc_supported.all()
    assert (rows[rows.maximum_chord<.18].bars==1).all()
    rows.to_csv(RESULTS/'ellipse/cells.csv',index=False);sample.to_csv(RESULTS/'ellipse/samples.csv',index=False)
    series.to_csv(RESULTS/'ellipse/series.csv',index=False);summary.to_csv(RESULTS/'ellipse/summary.csv',index=False)
    report=dict(samples=len(sample),series=len(series),cells=len(rows),condition_5c=int(sample.condition_5c.sum()),
        nn_crust_exact=int(sample.nn_crust_exact_ring.sum()),tau_E_min=float(sample.tau_E.min()),tau_E_max=float(sample.tau_E.max()),
        unsupported_cells=int((~rows.arc_supported).sum()),multiple_bar_cells=int((rows.bars>1).sum()),
        maximum_bars=int(rows.bars.max()),ordinary_decreases=int(series.ordinary_decreases.sum()),capped_decreases=int(series.capped_decreases.sum()),
        first_nonlocal_support_failures=int((~series.first_nonlocal_supported).sum()),
        local_cells=int((rows.maximum_chord<.18).sum()),all_local_supported=True,
        failed_5c_samples=sample[~sample.condition_5c].to_dict('records'))
    (RESULTS/'ellipse/report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

def noise_summary():
    rows=pd.read_csv(RESULTS/'noise/cells.csv');selectors=pd.read_csv(RESULTS/'noise/selectors.csv')
    assert len(selectors)==20*2*6*3*2
    summary=selectors.groupby(['model','delta','method','selector']).agg(samples=('seed','count'),selected=('selected','sum'),
        before_closure=('before_closure','sum'),median_M=('M','median'),median_e_PH=('e_PH','median'),
        q1_e_PH=('e_PH',lambda v:v.quantile(.25,interpolation='linear')),
        q3_e_PH=('e_PH',lambda v:v.quantile(.75,interpolation='linear')),
        median_longest_over_reference=('longest_over_reference','median')).reset_index()
    summary.to_csv(RESULTS/'noise/selector_summary.csv',index=False)
    screens=[]
    for landmark in ['first_connected','before_closure','closure']:
        for key,g in rows[rows[landmark]].groupby(['model','delta','method']):
            assert len(g)==20
            screens.append(dict(zip(['model','delta','method'],key))|dict(event=landmark,samples=len(g),empty=int((g.bars==0).sum()),
                quarter=int((g.longest_over_reference>=.25).sum()),half=int((g.longest_over_reference>=.5).sum()),
                three_quarters=int((g.longest_over_reference>=.75).sum()),median_ratio=float(g.longest_over_reference.median()),
                median_e_PH=float(g.e_PH.median())))
    pd.DataFrame(screens).to_csv(RESULTS/'noise/screen_summary.csv',index=False)
    curves=rows[rows.grid].groupby(['model','delta','method','M']).agg(samples=('seed','count'),connected=('connected','sum'),
        median_e_PH=('e_PH','median'),median_e_inf=('e_inf','median'),median_e_2=('e_2','median')).reset_index()
    assert (curves.samples==20).all()
    curves.to_csv(RESULTS/'noise/curves.csv',index=False)
    clean=selectors[selectors.delta==0]
    assert clean.selected.all() and not clean.before_closure.any()
    for key,g in rows[(rows.delta==0)&rows.grid].groupby(['seed','model','method']):
        first_closed=int(g[g.M>=g.M_close].M.min())
        selected=clean[(clean.seed==key[0])&(clean.model==key[1])&(clean.method==key[2])]
        assert len(selected)==2 and (selected.M==first_closed).all()
    a=rows[(rows.delta==.01)&(rows.model=='radial')&rows.first_connected]
    return dict(cells=len(rows),grid_cells=int(rows.grid.sum()),selector_records=len(selectors),clean_selector_records=len(clean),distinct_clean_cases=len(clean.drop_duplicates(['seed','method'])),distinct_clean_selections=len(clean.drop_duplicates(['seed','method','selector'])),
                clean_premature=0,largest_delta=.01)

def curves_plot():
    d=pd.read_csv(RESULTS/'noise/curves.csv');norm=512*np.log(512)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(7.0,2.65),sharey=True)
    for ax,model in zip(axes,['radial','isotropic']):
        for method in METHODS:
            g=d[(d.model==model)&(d.delta==.01)&(d.method==method)].sort_values('M')
            y=g.median_e_PH.replace([np.inf,-np.inf],np.nan)
            ax.plot(g.M/norm,y,label=NAMES[method],color=COLORS[method],lw=1.6)
        ax.set(xscale='log',yscale='log',xlabel=r'$M/(n\log n)$',title=model.capitalize()+r', $\delta=0.01$')
        ax.grid(alpha=.2)
    axes[0].set_ylabel(r'Median $e_{\mathrm{PH}}$');axes[1].legend(frameon=False,fontsize=8)
    fig.tight_layout();fig.savefig(FIGURES/'noise-eph.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(FIGURES/'noise-eph.png',dpi=300);plt.close(fig)
    fig,axes=plt.subplots(3,2,figsize=(7.0,5.0),sharex=True)
    for col,model in enumerate(['radial','isotropic']):
        for row,(metric,label) in enumerate([('median_e_PH',r'$e_{\mathrm{PH}}$'),('median_e_inf',r'$e_\infty$'),('median_e_2',r'$e_2$')]):
            ax=axes[row,col]
            for delta,style in [(0.,':'),(.001,'--'),(.01,'-')]:
                for method in METHODS:
                    g=d[(d.model==model)&(d.delta==delta)&(d.method==method)].sort_values('M')
                    y=g[metric].replace([np.inf,-np.inf],np.nan)
                    ax.plot(g.M/norm,y,color=COLORS[method],ls=style,lw=1.2)
            ax.set(xscale='log',yscale='log');ax.grid(alpha=.2)
            if col==0:ax.set_ylabel('Median '+label)
            if row==0:ax.set_title(model.capitalize())
            if row==2:ax.set_xlabel(r'$M/(n\log n)$')
    from matplotlib.lines import Line2D
    handles=[Line2D([0],[0],color=COLORS[k],label=NAMES[k]) for k in METHODS]
    handles += [Line2D([0],[0],color='gray',ls=ls,label=r'$\delta='+str(v)+'$') for v,ls in [(0,':'),(.001,'--'),(.01,'-')]]
    fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.5,-.015),fontsize=8)
    fig.tight_layout(rect=(0,.055,1,1));fig.savefig(FIGURES/'noise-all-errors.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(FIGURES/'noise-all-errors.png',dpi=300);plt.close(fig)

def unsupported_panel():
    plt.rcParams.update({'font.family':'serif','font.serif':['Latin Modern Roman','DejaVu Serif'],'mathtext.fontset':'cm','pdf.fonttype':42})
    values=np.array([0,.051,.442,.491,.543,.953]);xy=np.c_[np.cos(2*np.pi*values),np.sin(2*np.pi*values)]
    fig,ax=plt.subplots(figsize=(2.05,1.90))
    angle=np.linspace(0,2*np.pi,200);ax.plot(np.cos(angle),np.sin(angle),color='.78',lw=.6)
    for i in range(6):
        j=(i+1)%6;ax.plot(xy[[i,j],0],xy[[i,j],1],color='.2',lw=1)
    ax.plot(xy[[0,3],0],xy[[0,3],1],color='#d6604d',lw=1.8)
    ax.scatter(xy[:,0],xy[:,1],s=12,color='black',zorder=3)
    ax.set(xlim=(-1.2,1.2),ylim=(-1.2,1.2));ax.set_aspect('equal');ax.axis('off');ax.set_title('Unsupported shortcut',fontsize=9,fontweight='bold',pad=6)
    ax.text(.5,-.12,'Two loop bars',transform=ax.transAxes,ha='center',fontsize=9)
    fig.subplots_adjust(top=.79,bottom=.20);fig.savefig(FIGURES/'unsupported-shortcut.pdf',bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});fig.savefig(FIGURES/'unsupported-shortcut.png',dpi=300,bbox_inches='tight');plt.close(fig)

def main(make_figures=True):
    FIGURES.mkdir(exist_ok=True)
    report={'ellipse':ellipse_summary()}
    if (RESULTS/'noise/selectors.csv').exists():
        report['noise']=noise_summary()
        if make_figures:curves_plot()
    if make_figures:unsupported_panel()
    minimum=float(report['ellipse']['tau_E_min'])
    minima=[float(json.loads((RESULTS/'number_audit.json').read_text())['shortest_original_stored_positive_lifetime'])]
    for kind in ['ellipse','noise']:
        for p in (RESULTS/kind/'checkpoints').glob('*.csv'):
            d=pd.read_csv(p);v=d.minimum_positive_lifetime.dropna()
            if len(v):minima.append(float(v.min()))
    report['shortest_stored_positive_lifetime']=min(minima)
    (RESULTS/'study_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
