"""Export connectivity-based choices and fixed-k failure partitions from saved ranks."""
from pathlib import Path
import csv, json
ROOT = Path(__file__).resolve().parents[1]
PROFILES = ('uniform', 'strong_cosine', 'localized_trough')
METHODS = ('OR_kNN', 'mutual_kNN')

def records(root):
    root = Path(root)
    with (root/'results/extra_k/rank_sweep/samples.csv').open() as f:
        ranks = list(csv.DictReader(f))
    with (root/'results/comparisons/bridgeless/sweep_samples.csv').open() as f:
        bridge = {(v['profile'],v['method'],v['n'],v['seed']):v for v in csv.DictReader(f)}
    rows = []
    def count(group, choice):
        connected, closed = [], []
        for v in group:
            kc, kl = int(v['k_conn']), int(v['k_close'])
            b = bridge[v['profile'],v['method'],v['n'],v['seed']]
            assert int(b['k_conn']) == kc and int(b['k_close']) == kl
            k = {'connection':kc, 'bridge_free':int(b['k_br']),
                 'plus3':kc+3, 'plus5':kc+5, 'fixed10':10}[choice]
            connected.append(k >= kc); closed.append(k >= kl)
        assert all(not c or z for c,z in zip(closed,connected))
        return dict(profile=group[0]['profile'], n=int(group[0]['n']), method=group[0]['method'],
                    choice=choice, samples=len(group), disconnected=sum(not z for z in connected),
                    connected_unclosed=sum(z and not c for z,c in zip(connected,closed)),
                    closed=sum(closed), closed_percent=100*sum(closed)/len(group))
    for profile in PROFILES:
        for method in METHODS:
            group = [v for v in ranks if v['profile']==profile and v['method']==method and int(v['n'])==512]
            for choice in ('connection','bridge_free','plus3','plus5'):
                rows.append(count(group,choice))
    for n in (512,4096):
        for method in METHODS:
            group = [v for v in ranks if v['profile']=='uniform' and v['method']==method and int(v['n'])==n]
            rows.append(count(group,'fixed10'))
    assert len(rows)==28
    assert all(v['samples']==1500 and sum(v[k] for k in ('disconnected','connected_unclosed','closed'))==1500 for v in rows)
    return rows

def verify_practical_choices(root, require):
    root = Path(root)
    expected = records(root)
    with (root/'results/practical_choices/choices.csv').open() as f:
        actual = list(csv.DictReader(f))
    require('practical_choice_partitions', len(actual)==len(expected)==28)
    for i,(a,b) in enumerate(zip(actual,expected)):
        require(f'practical_choice_partition_{i}', all(a[k]==str(v) for k,v in b.items()), b)
    fixed = {}
    with (root/'results/retention_table/retention_fixed_k.csv').open() as f:
        for v in csv.DictReader(f):
            fixed[v['profile'],v['method'],int(v['n']),int(v['k'])]=v
    for v in expected:
        if v['choice']=='fixed10':
            other = fixed[v['profile'],v['method'],v['n'],10]
            require(f"practical_choice_fixed_crosscheck_{v['n']}_{v['method']}",
                    v['closed']==int(other['closed']) and
                    v['disconnected']==v['samples']-int(other['connected']))
    return dict(status='passed', rows=28, partitions='Disconnected, connected but unclosed, and closed; all attempted samples.')

def main():
    rows = records(ROOT)
    out = ROOT/'results/practical_choices';out.mkdir(exist_ok=True)
    with (out/'choices.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (out/'summary.json').write_text(json.dumps(dict(status='passed',rows=len(rows),measurement_generation=False),indent=2)+'\n')
    print('Practical choices: 24 connectivity-based cells and four fixed-k partitions.')

if __name__=='__main__':
    main()
