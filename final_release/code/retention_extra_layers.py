"""Export retention at k_conn+j, j=0..5, from the saved rank records.

No new sample or persistence diagram is generated. The main figure uses
n=512; the CSV contains all four primary sample sizes and pointwise Wilson
intervals. Methods and profiles share seeds, as in Table 3.
"""
from pathlib import Path
import csv,hashlib,json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
SIZES=(128,512,4096,32768)
PROFILES=('uniform','strong_cosine','localized_trough')
METHODS=('OR_kNN','mutual_kNN')
def wilson(k,n):
    p=k/n;z=1.959963984540054;d=1+z*z/n;c=(p+z*z/(2*n))/d
    h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d;return c-h,c+h
def main():
    source=ROOT/'results/extra_k/rank_sweep/samples.csv'
    records=list(csv.DictReader(source.open()));rows=[]
    for n in SIZES:
        for profile in PROFILES:
            for method in METHODS:
                group=[v for v in records if int(v['n'])==n and v['profile']==profile and v['method']==method]
                assert len(group)==(500 if n==32768 else 1500)
                assert len({v['seed'] for v in group})==len(group)
                for j in range(6):
                    k=sum(int(v['k_close'])<=int(v['k_conn'])+j for v in group);lo,hi=wilson(k,len(group))
                    rows.append(dict(n=n,profile=profile,method=method,j=j,closed=k,samples=len(group),rate=k/len(group),wilson_lo=lo,wilson_hi=hi))
    out=ROOT/'results/retention_extra_layers';out.mkdir(exist_ok=True)
    with (out/'retention.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':300})
    fig,axes=plt.subplots(1,2,figsize=(6.8,2.65),sharey=True,layout='constrained')
    colors=('#0072B2','#D55E00','#009E73');names=('Uniform','Cosine','Trough')
    for ax,method,title in zip(axes,METHODS,('OR kNN','Mutual kNN')):
        for profile,color,label in zip(PROFILES,colors,names):
            group=[v for v in rows if v['n']==512 and v['profile']==profile and v['method']==method]
            ax.plot([v['j'] for v in group],[100*v['rate'] for v in group],'-o',color=color,label=label,ms=3.5,lw=1.3)
            ax.fill_between([v['j'] for v in group],[100*v['wilson_lo'] for v in group],[100*v['wilson_hi'] for v in group],color=color,alpha=.12,lw=0)
        ax.set_title(title);ax.set_xlabel(r'Extra layers $j$ beyond $k_{\rm conn}$');ax.set_xticks(range(6));ax.set_xlim(-.1,5.1);ax.set_ylim(0,102);ax.grid(axis='y',alpha=.22)
    axes[0].set_ylabel('Ring retention (%)');axes[1].legend(loc='lower right',frameon=False,fontsize=9)
    (ROOT/'results/figure_exports').mkdir(parents=True,exist_ok=True)
    (ROOT/'manuscript_figures').mkdir(exist_ok=True)
    fig.savefig(ROOT/'results/figure_exports/retention_extra_layers.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(ROOT/'manuscript_figures/retention_extra_layers.png');plt.close(fig)
    config=dict(event='k_close <= k_conn+j',j=list(range(6)),sizes=list(SIZES),figure_n=512,
        source='results/extra_k/rank_sweep/samples.csv',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        interval='Pointwise two-sided 95% Wilson; profiles and methods paired within seed.',rows=len(rows),new_samples=0,new_diagrams=0)
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    print(json.dumps(dict(rows=len(rows),figure_n=512,j0_to_j5=[v for v in rows if v['n']==512 and v['j'] in [0,3,5]]),indent=2))
if __name__=='__main__':main()
