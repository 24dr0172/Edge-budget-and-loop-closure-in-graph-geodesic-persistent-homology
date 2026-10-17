"""Rebuild study comparisons from saved per-sample measurements."""
from pathlib import Path
import hashlib,json,sys
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
F=Path(__file__).resolve().parents[1]
OUT=F/'results/phase1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def oracle():
    cells=pd.read_csv(F/'results/noise/cells.csv')
    cfg=json.loads((F/'config.json').read_text())
    grid=cells[cells.grid]
    assert len(grid)==17280 and set(grid.M)==set(cfg['noise']['grid_M'])
    keys=['seed','model','delta','method']
    valid=grid[grid.connected&np.isfinite(grid.e_PH)]
    best=valid.sort_values(['e_PH','M'],kind='stable').groupby(keys,sort=True).head(1)
    assert len(best)==720
    columns=keys+['M','e_PH','bars','longest','M_close']
    best=best[columns].rename(columns={'M':'oracle_M','e_PH':'oracle_e_PH'})
    best['oracle_before_latent_closure']=best.oracle_M<best.M_close
    best.sort_values(keys).to_csv(OUT/'grid_oracle_samples.csv',index=False)
    summaries=[]
    for key,g in best.groupby(['model','delta','method']):
        summaries.append(dict(zip(['model','delta','method'],key),samples=len(g),
                              median_e_PH=g.oracle_e_PH.median(),
                              q1_e_PH=g.oracle_e_PH.quantile(.25,interpolation='linear'),
                              q3_e_PH=g.oracle_e_PH.quantile(.75,interpolation='linear'),
                              median_M=g.oracle_M.median(),
                              median_a=float(np.median(g.oracle_M/(512*np.log(512)))),
                              before_latent_closure=int(g.oracle_before_latent_closure.sum())))
    pd.DataFrame(summaries).to_csv(OUT/'grid_oracle_summary.csv',index=False)
    selected=pd.read_csv(F/'results/noise/selectors.csv')
    paired=selected.merge(best,on=keys,validate='many_to_one',suffixes=('','_oracle'))
    paired['selected_minus_oracle']=paired.e_PH-paired.oracle_e_PH
    assert (paired.selected_minus_oracle>=-1e-12).all()
    paired.to_csv(OUT/'selected_vs_oracle.csv',index=False)
    raw=pd.read_csv(F/'results/noise/selector_summary.csv')
    combined=raw.merge(pd.DataFrame(summaries),on=['model','delta','method'],
                       suffixes=('_selected','_oracle'),validate='many_to_one')
    combined.to_csv(OUT/'selector_oracle_summary.csv',index=False)


def amplitudes():
    selected=pd.read_csv(F/'results/noise/selectors.csv')
    rows=[]
    for key,g in selected.groupby(['model','delta','method']):
        row=dict(zip(['model','delta','method'],key));row['samples']=g.seed.nunique()
        for rule,h in g.groupby('selector'):
            row.update({rule+'_median_e_PH':h.e_PH.median(),
                        rule+'_q1_e_PH':h.e_PH.quantile(.25,interpolation='linear'),
                        rule+'_q3_e_PH':h.e_PH.quantile(.75,interpolation='linear'),
                        rule+'_median_a':float(np.median(h.M/(512*np.log(512))))})
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT/'all_amplitudes.csv',index=False)


def rank_counts():
    source=F/'results/extra_k/candidate_audit/samples.csv'
    rows=pd.read_csv(source)
    subsets={'primary_display_through_32768':rows[rows.n<=32768],'full_release':rows}
    records=[]
    for label,g in subsets.items():
        records.append(dict(subset=label,sizes=sorted(map(int,g.n.unique())),
                            seed_size_combinations=len(g[['seed','n']].drop_duplicates()),
                            profile_seed_samples=len(g[['profile','seed','n']].drop_duplicates()),
                            method_records=len(g),
                            repeats_per_profile_by_size={str(int(n)):int(len(h)//6) for n,h in g.groupby('n')}))
    assert [r['method_records'] for r in records]==[30000,30780]
    (OUT/'rank_display_counts.json').write_text(json.dumps(dict(status='passed',source_sha256=sha(source),subsets=records),indent=2)+'\n')


def literature():
    # Source locations and decisions are explicit inputs, not automatically inferred.
    # Rerunning reproduces the record; empirical claims remain in the study files.
    record=dict(
        reviewed_on='2026-10-16',
        connectivity_stopping=dict(
            decision='Choosing the smallest connecting parameter is described as a simple heuristic.',
            stronger_practice_claim_supported=False,
            reason='The checked software documentation uses a fixed neighbor default; connectivity requirements alone do not establish choosing the smallest connecting parameter.',
            source='https://scikit-learn.org/stable/modules/generated/sklearn.manifold.Isomap.html',
            scope='This is a wording decision from a bounded literature check, not a claim that no such reference exists.'),
        spectral_graph=dict(
            citation_key='damrich2024',
            source='https://papers.nips.cc/paper_files/paper/2024/file/4a32a646254d2e37fc74a38d65796552-Paper-Conference.pdf',
            locations=['Section 5','Appendix I, effective resistance and diffusion distance'],
            graph='Symmetric OR kNN; the reported spectral comparisons use unweighted graphs without self-loops.',
            parameters=dict(k=[15,100],diffusion_time=[8,64]),
            scope='Reported robustness to shortcut edges is empirical; no immunity or transferred shortest-path theorem is asserted.'),
        cyclic_graphs=dict(
            citation_key='adams2019',
            source='https://arxiv.org/html/1812.03374v2',
            locations=['Section 3 definition','Theorem 1','Section 4.2'],
            comparison='An inclusion filtration in one cyclic order with consistent directed edges; the Example 1 certificate rules out even varying directions in a fixed order.'))
    (OUT/'literature_checks.json').write_text(json.dumps(record,indent=2)+'\n')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    oracle();amplitudes();rank_counts();literature()
    source_names=['config.json','results/noise/cells.csv','results/noise/selectors.csv',
                  'results/noise/selector_summary.csv','results/extra_k/candidate_audit/samples.csv']
    (OUT/'summary_inputs.json').write_text(json.dumps({name:sha(F/name) for name in source_names},indent=2)+'\n')
    print('Saved: 720 sample grid oracles, 36 oracle groups, 36 all-amplitude rows, paired rule gaps, sweep counts and literature decisions.')


if __name__=='__main__':main()
