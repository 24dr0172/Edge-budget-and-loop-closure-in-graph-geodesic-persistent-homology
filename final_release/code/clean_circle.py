"""Recompute the clean-circle full-diagram checks at the reported events.

By default replay one paired sample. Use --all for all 60 point sets; all
three orderings and both nearest-neighbor threshold rules are checked.
"""
from pathlib import Path
import argparse
import csv
import numpy as np
import graph_methods as graph

ROOT=Path(__file__).resolve().parents[1]
PROFILES=('uniform','strong_cosine','localized_trough')
METHODS=('R2_global','OR_kNN','mutual_kNN')


def evaluate(seed,profile,n=512):
    s,xy=graph.sample_cloud(seed,profile,n=n)
    _,ii,jj,lengths,orders,ranks=graph.build_orders(xy)
    results=[]
    for method in METHODS:
        order=orders[method]
        conn=graph.first_connected(order,ii,jj,n)
        close=graph.closure(order,ii,jj,n)
        budgets={'first_connected':conn,'pre_completion':close-1,'completion':close}
        if method in ranks:
            keys=ranks[method][order]; k=int(keys[conn-1])
            budgets['smallest_connected_k']=int(np.searchsorted(keys,k,side='right'))
        cache={}
        for event,M in budgets.items():
            if M not in cache:
                D=graph.graph_metric(order,ii,jj,lengths,M,n)
                if not np.isfinite(D).all():
                    raise AssertionError(f'{method} {event} graph is disconnected')
                bars=graph.finite_graph_diagram(D,float(lengths[order[:M]].max()))
                cache[M]=len(bars)
            results.append(dict(seed=seed,profile=profile,method=method,event=event,
                                M=M,connected=True,h1_bars=cache[M]))
    return results,float(np.max(np.diff(np.r_[s,s[0]+1])))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=int,default=20310101)
    p.add_argument('--profile',choices=PROFILES,default='uniform')
    p.add_argument('--all',action='store_true',help='Replay all 60 saved point sets')
    args=p.parse_args()
    path=ROOT/'results/clean_closure/events.csv'
    with path.open(newline='') as f: saved=list(csv.DictReader(f))
    jobs=[(seed,profile) for profile in PROFILES for seed in range(20310101,20310121)] if args.all else [(args.seed,args.profile)]
    checked=0
    for seed,profile in jobs:
        actual,gap=evaluate(seed,profile)
        expected={(r['method'],r['event']):r for r in saved if int(r['seed'])==seed and r['profile']==profile}
        assert len(expected)==len(actual)==11 and gap<1/3
        for row in actual:
            other=expected[row['method'],row['event']]
            assert row['M']==int(other['M']) and row['h1_bars']==int(float(other['h1_bars']))
            checked+=1
        print(f'{profile} seed {seed}: 11 event checks passed',flush=True)
    print(f'Clean-circle event replay passed: {checked} event checks.')

if __name__=='__main__':main()
