"""Generate clean, noise, ellipse and Euclidean measurements from seeds."""
import argparse,concurrent.futures,json,time
from experiment_core import CONFIG,ROOT
from generate_clean import clean_job,aggregate
from run_studies import noise_job,ellipse_job

def main(workers=8):
    jobs=[]
    for seed in CONFIG['clean']['seeds']:
        for profile in CONFIG['clean']['profiles']: jobs.append(('clean',(seed,profile)))
    for seed in CONFIG['ellipse']['seeds']:
        for profile in CONFIG['ellipse']['profiles']: jobs.append(('ellipse',(seed,profile)))
    for seed in CONFIG['noise']['seeds']:
        for model in CONFIG['noise']['models']:
            for delta in CONFIG['noise']['amplitudes']:
                if delta==0 and model=='isotropic': continue
                jobs.append(('noise',(seed,model,delta)))
    functions={'clean':clean_job,'ellipse':ellipse_job,'noise':noise_job}
    start=time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(functions[kind],args) if kind=='clean' else pool.submit(functions[kind],*args):(kind,args) for kind,args in jobs}
        for i,f in enumerate(concurrent.futures.as_completed(futures),1):
            kind,args=futures[f]; result=f.result()
            print(i,len(jobs),kind,args,result,flush=True)
    aggregate()
    report={'completed_jobs':len(jobs),'workers':workers,'elapsed_seconds':time.monotonic()-start,'started_with':'Empty results; only this run’s checkpoints are reused.'}
    (ROOT/'results/generation_timing.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=8)
    main(p.parse_args().workers)
