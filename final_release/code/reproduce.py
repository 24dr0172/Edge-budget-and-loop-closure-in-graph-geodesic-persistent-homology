"""Run every experiment in a new three-folder release, starting with no results."""
from pathlib import Path
import argparse,hashlib,json,os,platform,shutil,subprocess,sys,time
sys.dont_write_bytecode=True
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
ROOT=Path(__file__).resolve().parents[1]

def source_hashes(root):
    return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((root/'code').rglob('*')) if p.is_file() and '__pycache__' not in p.parts}

def pipeline(root,workers,manuscript_dir=None):
    commands=[
        ('generate_ellipse.py',['--workers',str(workers)]),
        ('check_knn_miss_candidate.py',[]),
        ('generate_measurements.py',['--workers',str(workers)]),
        ('run_studies.py',['--aggregate-only']),
        ('generate_matched_noise.py',[]),
        ('generate_witness.py',[]),
        ('verify_separation.py',[]),
        ('verify_arc_support.py',[]),
        ('run_studies.py',['--validate-only']),
        ('reproduce_results.py',['--base-only']),
        ('phase1_bridgeless.py',['--validate','--clean','--sweep']),
        ('phase1_noise.py',[]),
        ('check_cknn_precision.py',[]),
        ('check_rank_pairing.py',[]),
        ('validate_cutoff.py',[]),
        ('benchmark_diagrams.py',[]),
        ('reproduce_results.py',[]),
        ('build_figures.py',[]),
        ('write_guides.py',[]),
        ('record_current.py',[]),
        ('verify_release.py',[])]
    if manuscript_dir is not None:
        commands=[(name,args+['--manuscript-dir',str(manuscript_dir)] if name in ('record_current.py','verify_release.py') else args) for name,args in commands]
    logdir=root/'results/run_logs';logdir.mkdir(parents=True,exist_ok=True)
    record_path=root/'results/reproduction_run.json'
    record=json.loads(record_path.read_text()) if record_path.exists() else dict(
        started_with_empty_results=True,scope='All current-paper scientific experiments, summaries, figures and checks.',
        source_sha256=source_hashes(root),python=platform.python_version(),stages=[])
    if record['source_sha256']!=source_hashes(root):
        raise ValueError('Code or protocol changed since this run began. Start in a new output folder.')
    completed={r['script'] for r in record['stages']}
    # The two summary calls and both run_studies calls have distinct arguments.
    completed_keys={r['script']+' '+ ' '.join(r['arguments']) for r in record['stages']}
    start=time.monotonic()
    for index,(name,args) in enumerate(commands,1):
        key=name+' '+' '.join(args)
        if key in completed_keys:
            print(f'[{index}/{len(commands)}] Already complete: {key}',flush=True);continue
        print(f'[{index}/{len(commands)}] {key}',flush=True)
        stage_start=time.monotonic();log=logdir/f'{index:02d}_{Path(name).stem}.log'
        with log.open('w') as stream:
            command=[sys.executable,'-B','-u',str(root/'code'/name),*args]
            with subprocess.Popen(command,cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                  text=True,bufsize=1) as process:
                for line in process.stdout:
                    print(line,end='',flush=True);stream.write(line);stream.flush()
                if process.wait():raise subprocess.CalledProcessError(process.returncode,command)
        record['stages'].append(dict(script=name,arguments=args,elapsed_seconds=time.monotonic()-stage_start,
                                     log=str(log.relative_to(root))))
        # The manifest and verifier must see the final, stable run record.
        if name not in ('record_current.py','verify_release.py'):
            record_path.write_text(json.dumps(record,indent=2)+'\n')
    # Do not rewrite evidence after manifest creation.
    print(f'Complete. Fresh results and six PDF/PNG figure pairs are in {root}',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True,help='New output folder; it will contain code, results and manuscript_figures.')
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--manuscript-dir',type=Path,help='Optional external manuscript directory for the final paper comparison.')
    p.add_argument('--resume',action='store_true',help='Continue the same run in this output folder.')
    args=p.parse_args();target=args.output.expanduser().resolve()
    if args.workers<1: p.error('--workers must be positive')
    if target==ROOT or ROOT in target.parents:p.error('Choose an output folder outside this release.')
    if args.resume:
        if not (target/'results/reproduction_run.json').is_file():p.error('No reproduction run exists at this output path.')
        if source_hashes(ROOT)!=source_hashes(target):p.error('The output code differs from the supplied release.')
    else:
        if target.exists() and any(target.iterdir()):p.error('Output is not empty; choose a new folder or use --resume.')
        target.mkdir(parents=True,exist_ok=True)
        shutil.copytree(ROOT/'code',target/'code',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        (target/'results').mkdir();(target/'manuscript_figures').mkdir()
        (target/'results/reproduction_run.json').write_text(json.dumps(dict(started_with_empty_results=True,
            scope='All current-paper scientific experiments, summaries, figures and checks.',
            source_sha256=source_hashes(target),python=platform.python_version(),stages=[]),indent=2)+'\n')
    # Execute from the copied code so all paths resolve inside the new release.
    script='from reproduce import pipeline; from pathlib import Path; pipeline(Path(__import__("sys").argv[1]),int(__import__("sys").argv[2]),__import__("sys").argv[3] or None)'
    subprocess.run([sys.executable,'-B','-u','-c',script,str(target),str(args.workers),str(args.manuscript_dir.resolve()) if args.manuscript_dir else ''],cwd=target/'code',check=True)
if __name__=='__main__':main()
