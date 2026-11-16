"""Reproduce the independent uniform-circle tied-rank comparison."""
from pathlib import Path
import concurrent.futures, csv, hashlib, json, subprocess, sys, time
import numpy as np

RELEASE = Path(__file__).resolve().parents[1]
HERE = RELEASE/'results/detection_rank_sweep/uniform_100k'
sys.path.insert(0, str(RELEASE/'code'))
from rank_sweep import neighboring_keys, tails, wilson

BIN = '/tmp/targeted-rank-sweep'
SIZES = (128, 512, 4096, 32768)
REPS = 100000

def main():
    started = time.monotonic()
    rng = np.random.default_rng(2026111501)
    cases = [rng.exponential(size=n) for n in (8, 32, 128, 512, 4096) for _ in range(20)]
    cases = [z for z in cases if z.max()/z.sum() < .5]
    inputs = ''.join(str(len(z))+' '+' '.join(format(x,'.17g') for x in z)+'\n' for z in cases)
    out = subprocess.run([BIN,'--keys'],input=inputs,text=True,capture_output=True,check=True).stdout.splitlines()
    for z, line in zip(cases,out):
        s = np.r_[0.,np.cumsum(z[:-1])]/z.sum()
        expected,_,_ = neighboring_keys(s)
        actual = np.fromstring(line,sep=' ',dtype=int)
        assert np.array_equal(actual[:len(z)],expected['OR_kNN'])
        assert np.array_equal(actual[len(z):],expected['mutual_kNN'])
    assert len(out)==len(cases)
    print('Independent C++/released rank comparisons:',len(cases),'point sets',flush=True)
    probs=tails(120)
    jobs=[(n,REPS,202611150000+n) for n in SIZES]
    def one(job):
        line=subprocess.run([BIN,*map(str,job)],capture_output=True,text=True,check=True).stdout.strip()
        print('Completed',line,flush=True)
        return line
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        raw=list(pool.map(one,jobs))
    rows=[]
    for line in raw:
        n,reps,seed,ts,tm,bad,g=line.split(',')
        n,reps,seed,ts,tm,bad=map(int,(n,reps,seed,ts,tm,bad))
        for method,count in [('OR_kNN',ts),('mutual_kNN',tm)]:
            pi=probs[method]
            approx=1-float(np.sum(n*(pi[:-1]-pi[1:])*np.exp(-n*pi[:-1])))
            lo,hi=wilson(count,reps)
            rows.append(dict(n=n,samples=reps,seed=seed,method=method,ties=count,
                tie_rate=count/reps,formula=approx,difference_percentage_points=100*(count/reps-approx),
                standard_error=np.sqrt((count/reps)*(1-count/reps)/reps),wilson_lo=lo,wilson_hi=hi,
                gap_at_least_third=bad,largest_observed_gap=float(g)))
    with (HERE/'rank_sweep.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (HERE/'rank_protocol.json').write_text(json.dumps(dict(samples_per_size=REPS,sizes=SIZES,
        independent_kernel_checks=len(cases),generator='std::mt19937_64; open-interval uniform to Exp(1); normalized circular spacings',
        seed_rule='202611150000+n; independent streams by size; paired rules within each point set',
        source_sha256=hashlib.sha256((RELEASE/'code/uniform_rank_sweep.cpp').read_bytes()).hexdigest(),
        runtime_seconds=time.monotonic()-started,results=rows),indent=2)+'\n')
    print(json.dumps(rows,indent=2),flush=True)

if __name__=='__main__':main()
