"""Optional recomputation of retained trials; never used by saved-table extraction.
Run from the release with --study clean/noise/useful/orders/ellipse and --output PATH.
Output must be inside paper 1 and outside final_release.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import paper1_methods as m

def regenerate(root,output,study,seed):
    paper=root.parent
    first=20320101 if study=='noise' else 20310101
    if seed not in range(first,first+20):
        raise ValueError('Choose one of the retained 20 seed streams for this study')
    output=output.resolve()
    if not output.is_relative_to(paper) or output==root or output.is_relative_to(root):
        raise ValueError('Use scratch output inside paper 1, outside final_release')
    output.mkdir(parents=True,exist_ok=True);rows=[]
    if study in ['clean','orders','useful']:
        for profile in m.PROFILES:
            path=root/f'results/useful_diagram/checkpoints/{profile}_{seed}.npz'
            with np.load(path) as z:s=z['latent'];xy=z['observed'];ref=z['reference']
            true=m.intrinsic(s);chord,iu,ju,length,orders,keys=m.build_orders(xy,include_cknn=study=='orders')
            if study=='orders':
                for method,order in orders.items():
                    cover,holds,checked,witness=m.certificate(s,chord,order,iu,ju)
                    rows.append(dict(seed=seed,profile=profile,method=method,M_cov=cover,holds=holds,checked=checked,**witness))
                continue
            source=root/('results/clean_closure/cells.csv' if study=='clean' else 'results/useful_diagram/saved_prefixes.csv')
            cells=pd.read_csv(source);cells=cells[(cells.seed==seed)&(cells.profile==profile)]
            for method in m.METHODS:
                for M in sorted(set(cells[cells.method==method].M)):
                    d=m.graph_metric(orders[method],iu,ju,length,int(M),len(s))
                    rows.append(dict(seed=seed,profile=profile,method=method,M=int(M),**m.score(d,true,ref)))
    elif study=='ellipse':
        checks=pd.read_csv(root/'results/ellipse_checks/checks.csv')
        with np.load(root/'results/ellipse_checks/coordinates.npz') as z:
            for profile in m.PROFILES:
                s=z[f'{seed}_{profile}_latent'];xy=z[f'{seed}_{profile}_xy']
                true=m.intrinsic(s);ref=m.diagram(true)
                chord,iu,ju,length,orders,keys=m.build_orders(xy)
                for method in m.METHODS:
                    row=checks[(checks.seed==seed)&(checks.profile==profile)&(checks.method==method)].iloc[0]
                    for M in sorted(set([int(row.M_conn),int(row.M_close)-1,int(row.M_close)])):
                        d=m.graph_metric(orders[method],iu,ju,length,M,len(s))
                        rows.append(dict(seed=seed,profile=profile,method=method,M=M,**m.score(d,true,ref)))
    else:
        config=json.loads((root/'results/matched_noise/config.json').read_text());s,clean=m.sample_cloud(seed,'uniform')
        true=m.intrinsic(s);ref=m.diagram(true);normal=clean*2*np.pi
        vector=np.random.default_rng(seed+m.NOISE_SEED_OFFSET).uniform(-1,1,len(s))
        for delta in config['delta']:
            xy=clean+delta*vector[:,None]*normal;chord,iu,ju,length,orders,keys=m.build_orders(xy)
            for method in config['methods']:
                M=m.first_connected(orders[method],iu,ju,len(s));d=m.graph_metric(orders[method],iu,ju,length,M,len(s))
                rows.append(dict(seed=seed,delta=delta,method=method,M=M,**m.score(d,true,ref)))
    pd.DataFrame(rows).to_csv(output/f'{study}_{seed}.csv',index=False)
    print('Recomputed retained trials:',len(rows))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--release',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--study',choices=['clean','noise','useful','orders','ellipse'],required=True);ap.add_argument('--seed',type=int,required=True)
    ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();regenerate(a.release.resolve(),a.output,a.study,a.seed)
