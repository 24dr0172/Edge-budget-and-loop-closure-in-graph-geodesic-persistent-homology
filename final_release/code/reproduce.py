"""Regenerate the whole release from seeds in an empty folder, or check that this is possible (--check)."""
from pathlib import Path
import argparse,hashlib,json,os,platform,shutil,subprocess,sys,time
from importlib import metadata
sys.dont_write_bytecode=True
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
ROOT=Path(__file__).resolve().parents[1]

def source_hashes(root):
    return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted([p for folder in ['code','tools'] for p in (root/folder).rglob('*')]) if p.is_file() and '__pycache__' not in p.parts}

def stages(root,workers,manuscript_dir=None):
    commands=[
        ('generate_ellipse.py',['--workers',str(workers)]),
        ('rank_sweep.py',[]),
        ('generate_measurements.py',['--workers',str(workers)]),
        ('run_studies.py',['--aggregate-only']),
        ('generate_matched_noise.py',[]),
        ('generate_witness.py',[]),
        ('verify_separation.py',[]),
        ('verify_natural_instance.py',[]),
        ('verify_small_examples.py',[]),
        ('verify_arc_support.py',['--report',str(root/'results/arc_support_verification.json')]),
        ('verify_curve_order.py',['--workers',str(workers),'--two-sat']),
        ('run_studies.py',['--validate-only']),
        ('reproduce_results.py',['--base-only']),
        ('bridge_free.py',['--validate','--clean','--sweep']),
        ('euclidean_baseline.py',[]),
        ('check_cknn_precision.py',[]),
        ('check_rank_pairing.py',[]),
        ('validate_cutoff.py',[]),
        ('benchmark_diagrams.py',[]),
        ('reproduce_results.py',[]),
        ('build_figures.py',[]),
        ('check_practical_claims.py',[]),
        ('build_retention_table.py',[]),
        ('plot_post_closure_error.py',[]),
        ('check_seven_point_example.py',[]),
        ('verify_retention_cost_witness.py',['--write']),
        ('check_noncircular_locality.py',['--workers',str(workers)]),
        ('retention_extra_layers.py',[]),
        ('check_peanut_trough.py',['--workers',str(workers)]),
        ('verify_extra_layers_and_peanut.py',[]),
        ('build_practical_choices.py',[]),
        ('peanut_neck_study.py',['--workers',str(workers)]),
        ('tools/write_guides.py',[]),
        ('verify_comparisons.py',[]),
        ('tools/record_current.py',[]),
        ('verify_release.py',[])]
    comparison_args = (['--manuscript-dir', str(manuscript_dir)]
                       if manuscript_dir is not None else ['--data-only'])
    return [(name, args+comparison_args if name in ('tools/record_current.py','verify_release.py') else args)
            for name,args in commands]


def preflight(root,workers,manuscript_dir=None):
    """Check archive files and public dependencies without generating data."""
    root=Path(root)
    if sys.version_info[:2] != (3,12):
        raise RuntimeError('Use Python 3.12 with code/requirements.txt.')
    commands=stages(root,workers,manuscript_dir)
    required=[root/name if name.startswith('tools/') else root/'code'/name for name,_ in commands]
    required += [root/name for name in ('code/config.json','code/requirements.txt',
                  'README.md','LICENSE','CITATION.cff','Dockerfile',
                  'code/bridgeless_circle.cpp','code/result_schemas.json',
                  'code/graph_methods.py','code/experiment_core.py','code/finite_diagram.py',
                  'code/rank_table.py','code/verify_paper_evidence.py','code/verify_noncircular_locality.py',
                  'code/verify_reported_claims.py','code/comparison_report.py','code/comparison_summaries.py',
                  'code/verify_peanut_neck.py',
                  'code/generate_clean.py','code/summarize_studies.py','code/verify_studies.py',
                  'code/audit_rank_scaling.py','code/audit_numbers.py','code/export_noise_coordinates.py',
                  'tools/external_manuscript.py','tools/verify_manuscript_tables.py')]
    required += [root/'manuscript_figures'/f'{name}.png' for name in ('closure_schematic','arc_support_fan')]
    required += [root/'results/figure_exports'/f'{name}.pdf' for name in ('closure_schematic','arc_support_fan')]
    required += [root/'results/additional_figures'/f'{name}.{ext}' for name in ('arc_support_figure','ellipse_figure') for ext in ('pdf','png')]
    missing=sorted({str(path.relative_to(root)) for path in required if not path.is_file()})
    if missing:raise FileNotFoundError('Release is incomplete: '+', '.join(missing))
    versions={}
    for line in (root/'code/requirements.txt').read_text().splitlines():
        line=line.split('#',1)[0].strip()
        if not line:continue
        package,expected=line.split('==',1)
        try:actual=metadata.version(package)
        except metadata.PackageNotFoundError as exc:
            raise RuntimeError(f'Missing {package}; install code/requirements.txt.') from exc
        if actual != expected:
            raise RuntimeError(f'{package}=={expected} is required, found {actual}; use a clean virtual environment.')
        versions[package]=actual
    compiler=shutil.which('g++')
    if compiler is None:raise RuntimeError('Complete generation requires g++ with C++17 support (or use the Dockerfile).')
    if manuscript_dir is not None:
        paper=Path(manuscript_dir).expanduser().resolve()
        if paper==root.resolve() or root.resolve() in paper.parents:
            raise ValueError('Keep the optional manuscript outside the release.')
        for name in ('combined.tex','combined.pdf'):
            if not (paper/name).is_file():raise FileNotFoundError(f'Optional manuscript is missing {name}.')
    return dict(status='passed',stages=len(commands),packages=versions,compiler=compiler,
                manuscript_comparison='requested' if manuscript_dir is not None else 'not_requested')


def pipeline(root,workers,manuscript_dir=None):
    preflight(root,workers,manuscript_dir)
    commands=stages(root,workers,manuscript_dir)
    logdir=root/'results/run_logs';logdir.mkdir(parents=True,exist_ok=True)
    record_path=root/'results/reproduction_run.json'
    record=json.loads(record_path.read_text()) if record_path.exists() else dict(
        started_with_empty_results=True,scope='Complete release: evidence for the paper and the additional studies, generated from seeds.',
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
            command=[sys.executable,'-B','-u',str(root/name) if name.startswith('tools/') else str(root/'code'/name),*args]
            with subprocess.Popen(command,cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                  text=True,bufsize=1) as process:
                for line in process.stdout:
                    print(line,end='',flush=True);stream.write(line);stream.flush()
                if process.wait():raise subprocess.CalledProcessError(process.returncode,command)
        record['stages'].append(dict(script=name,arguments=args,elapsed_seconds=time.monotonic()-stage_start,
                                     log=str(log.relative_to(root))))
        # The run record is sealed before manifest creation. Mutable run_logs
        # are explicitly excluded from that manifest by recorder and verifier.
        if name not in ('tools/record_current.py','verify_release.py'):
            record_path.write_text(json.dumps(record,indent=2)+'\n')
    # Do not rewrite evidence after manifest creation.
    print(f'Complete. Regenerated data and figure assets are in {root}',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,help='New output folder for a complete regenerated release.')
    p.add_argument('--check',action='store_true',help='Validate all pipeline files, pinned dependencies and compiler without generating measurements.')
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--manuscript-dir',type=Path,help='Optional external manuscript directory for the manuscript comparison.')
    p.add_argument('--resume',action='store_true',help='Continue the same run in this output folder.')
    args=p.parse_args()
    if args.workers<1: p.error('--workers must be positive')
    try: readiness=preflight(ROOT,args.workers,args.manuscript_dir)
    except (RuntimeError,FileNotFoundError,ValueError) as exc:p.error(str(exc))
    if args.check:
        print(json.dumps(readiness,indent=2));return
    if args.output is None:p.error('--output is required unless --check is used')
    target=args.output.expanduser().resolve()
    if target==ROOT or ROOT in target.parents:p.error('Choose an output folder outside this release.')
    if args.resume:
        if not (target/'results/reproduction_run.json').is_file():p.error('No reproduction run exists at this output path.')
        if source_hashes(ROOT)!=source_hashes(target):p.error('The output code differs from the supplied release.')
    else:
        if target.exists() and any(target.iterdir()):p.error('Output is not empty; choose a new folder or use --resume.')
        target.mkdir(parents=True,exist_ok=True)
        shutil.copytree(ROOT/'code',target/'code',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        shutil.copytree(ROOT/'tools',target/'tools',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        for name in ['README.md','LICENSE','CITATION.cff','Dockerfile']:shutil.copy2(ROOT/name,target/name)
        (target/'results').mkdir();(target/'manuscript_figures').mkdir()
        # Native drawings are fixed exposition, not generated measurements.
        for name in ('closure_schematic','arc_support_fan'):
            shutil.copy2(ROOT/'manuscript_figures'/f'{name}.png',target/'manuscript_figures'/f'{name}.png')
        for folder,names in [('figure_exports',('closure_schematic','arc_support_fan')),('additional_figures',('arc_support_figure','ellipse_figure'))]:
            (target/'results'/folder).mkdir()
            for name in names:
                for path in (ROOT/'results'/folder).glob(name+'.*'):shutil.copy2(path,target/'results'/folder/path.name)
        (target/'results/reproduction_run.json').write_text(json.dumps(dict(started_with_empty_results=True,
            scope='Complete release: evidence for the paper and the additional studies, generated from seeds.',
            source_sha256=source_hashes(target),python=platform.python_version(),stages=[]),indent=2)+'\n')
    # Execute from the copied code so all paths resolve inside the new release.
    script='from reproduce import pipeline; from pathlib import Path; pipeline(Path(__import__("sys").argv[1]),int(__import__("sys").argv[2]),__import__("sys").argv[3] or None)'
    subprocess.run([sys.executable,'-B','-u','-c',script,str(target),str(args.workers),str(args.manuscript_dir.resolve()) if args.manuscript_dir else ''],cwd=target/'code',check=True)
if __name__=='__main__':main()
