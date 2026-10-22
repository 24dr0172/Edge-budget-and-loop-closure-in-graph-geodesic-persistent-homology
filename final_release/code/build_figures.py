"""Export additional-study figures from saved rows."""
from pathlib import Path
import argparse,os
import numpy as np
import pandas as pd

def build(root):
    os.environ.setdefault('MPLCONFIGDIR',str(Path(__import__('tempfile').gettempdir())/'edge-budget-matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    out=root/'results/additional_figures';out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':9.2,'axes.labelsize':9.2,'axes.titlesize':10,'xtick.labelsize':8.5,
                         'ytick.labelsize':8.5,'legend.fontsize':9,'pdf.fonttype':42,'ps.fonttype':42,'lines.linewidth':1.3})
    styles=[('R2_global','Radius','#0072B2','o','-'),('OR_kNN','OR kNN','#009E73','s','--'),
            ('mutual_kNN','Mutual kNN','#CC79A7','^','-.')]
    profiles=['uniform','strong_cosine','localized_trough'];titles=['Uniform','Cosine','Sparse trough']
    summary=pd.read_csv(root/'results/clean_closure/budget_summary.csv')
    fig,axes=plt.subplots(3,3,figsize=(6.6,6.3),sharex=True)
    for c,profile in enumerate(profiles):
        for method,label,color,marker,ls in styles:
            g=summary[(summary.profile==profile)&(summary.method==method)].sort_values('M')
            for row,field in [(0,'connected'),(1,'median_e_inf'),(2,'median_e_PH')]:
                y=(g.connected/g.attempted).to_numpy() if row==0 else g[field].to_numpy().copy()
                if row:y[~np.isfinite(y)|(y<=0)]=np.nan
                axes[row,c].plot(g.a_M,y,label=label,color=color,marker=marker,ls=ls,ms=3,markevery=3)
        axes[0,c].set_title(titles[c],pad=12)
        for row in range(3):
            ax=axes[row,c];ax.set_xscale('log');ax.set_xticks([.5,1,3,12],labels=['0.5','1','3','12'])
            ax.grid(alpha=.2,linewidth=.4);ax.spines[['top','right']].set_visible(False)
            ax.text(.02,1.015,chr(97+3*row+c),transform=ax.transAxes,fontweight='bold')
            if row:ax.set_yscale('log')
            else:ax.set_ylim(-.025,1.025);ax.set_yticks([0,.5,1])
        axes[2,c].set_xlabel(r'$a_M=M/(n\log n)$')
    axes[0,0].set_ylabel('Connected fraction');axes[1,0].set_ylabel(r'Median $e_\infty$');axes[2,0].set_ylabel(r'Median $e_{\mathrm{PH}}$')
    h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,ncol=3,loc='upper center',bbox_to_anchor=(.53,.995),frameon=False)
    fig.subplots_adjust(left=.12,right=.995,bottom=.085,top=.87,hspace=.33,wspace=.5)
    fig.savefig(out/'clean_budget_curves.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(out/'clean_budget_curves.png',dpi=300);plt.close(fig)
    noise=pd.read_csv(root/'results/matched_noise/summary.csv');fig,axes=plt.subplots(1,2,figsize=(6.6,2.7))
    for method,label,color,marker,ls in styles[:2]:
        g=noise[noise.method==method].sort_values('delta')
        axes[0].plot(g.delta_n,g.empty_H1/g.attempted,label=label,color=color,marker=marker,ls=ls,ms=4)
        axes[1].plot(g.delta_n,g.median_longest_over_reference,label=label,color=color,marker=marker,ls=ls,ms=4)
    axes[0].set_ylabel('Empty-diagram fraction');axes[0].set_ylim(-.03,1.03)
    axes[1].set_ylabel('Median longest lifetime /\nreference lifetime')
    for i,ax in enumerate(axes):
        ax.set_xlabel(r'Noise / mean spacing: $\delta n$');ax.grid(alpha=.2,linewidth=.4)
        ax.spines[['top','right']].set_visible(False);ax.text(.02,1.03,chr(97+i),transform=ax.transAxes,fontweight='bold')
    axes[0].legend(frameon=False,loc='upper right');fig.subplots_adjust(left=.095,right=.985,bottom=.23,top=.94,wspace=.5)
    fig.savefig(out/'radial_first_connectivity.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(out/'radial_first_connectivity.png',dpi=300);plt.close(fig)
    gaps=pd.read_csv(root/'results/extra_k/samples.csv');fig,axes=plt.subplots(1,3,figsize=(6.6,2.7))
    for c,profile in enumerate(profiles):
        for method,label,color,marker,ls in styles[1:]:
            values=gaps[(gaps.profile==profile)&(gaps.method==method)].extra_k.to_numpy()
            unique,counts=np.unique(values,return_counts=True)
            x=np.r_[-.25,unique,unique[-1]+.5];y=np.r_[0,np.cumsum(counts)/len(values),1]
            axes[c].step(x,y,where='post',label=label,color=color,ls=ls)
            axes[c].plot(unique,np.cumsum(counts)/len(values),color=color,marker=marker,ls='none',ms=3)
        maximum=int(gaps[gaps.profile==profile].extra_k.max());axes[c].set_xlim(-.25,maximum+.5)
        axes[c].set_ylim(-.025,1.025);axes[c].set_yticks([0,.5,1]);axes[c].set_title(titles[c])
        axes[c].xaxis.set_major_locator(MaxNLocator(nbins=5,integer=True));axes[c].set_xlabel(r'$k_{\mathrm{close}}-k_{\mathrm{conn}}$')
        axes[c].grid(alpha=.2,linewidth=.4);axes[c].spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Fraction at or below');h,l=axes[0].get_legend_handles_labels()
    fig.legend(h,l,ncol=2,loc='upper center',bbox_to_anchor=(.53,.995),frameon=False)
    fig.subplots_adjust(left=.09,right=.99,bottom=.23,top=.78,wspace=.30)
    fig.savefig(out/'extra_neighbor_distribution.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(out/'extra_neighbor_distribution.png',dpi=300);plt.close(fig)
    from summarize_studies import curves_plot,unsupported_panel
    curves_plot(root);unsupported_panel(root)
    print('Exported additional-study figure assets from saved observations.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1]);args=ap.parse_args();build(args.release.resolve())
