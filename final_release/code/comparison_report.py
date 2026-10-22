"""Write results/comparisons/README.md and results/comparisons/completion.json from saved summaries."""
from pathlib import Path
import argparse, hashlib, json, platform, sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from external_manuscript import locate_manuscript
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

F = Path(__file__).resolve().parents[1]
OUT = F/'results/comparisons'

def main(manuscript_dir=None,data_only=False):
    separation = json.loads((F/'results/separation_verification.json').read_text())['filtration']
    sweep = pd.read_csv(OUT/'bridgeless/sweep_summary.csv')
    timings = pd.read_csv(OUT/'bridgeless/sweep_timings.csv')
    euclidean = pd.read_csv(OUT/'euclidean/summary.csv')
    paired = pd.read_csv(OUT/'euclidean/paired_summary.csv')
    oracle = pd.read_csv(OUT/'grid_oracle_summary.csv')
    raw_oracle = pd.read_csv(OUT/'grid_oracle_samples.csv')
    lines = [
        '# Study comparisons',
        '',
        'This folder holds the bridge-free study, the Euclidean baseline and the grid oracle. The paper cites two bridge-free proportions at n=512 (Section 7.2). The other comparisons are additional studies that the paper does not rely on. The exact certificate for Example 1 is linked below for convenience. Links are relative to this folder.',
        '',
        '## Evidence',
        '',
        '- [Exact Example 1 certificate](../separation_verification.json): all 120 cyclic orders at 14 distinct threshold graphs. Each graph admits an order separately; no order works throughout the filtration. The certificate records distances and admissible orders, including the thresholds 3/2 and 19/10 inside the positive interval. Table 3 is exported separately to [rank_table.csv](../rank_table.csv); it is not the bridge-free summary.',
        '- [Clean bridge events](bridgeless/clean_samples.csv): 60 samples and two kNN rules, checked against full Cartesian ranks and an independent bridge search. The [summary](bridgeless/clean_summary.csv) gives rates and Wilson intervals.',
        '- [Bridge sweep](bridgeless/sweep_samples.csv): 15,390 profile/seed/size samples and 30,780 method records through n=1,000,000. See the [summary](bridgeless/sweep_summary.csv), [validation](bridgeless/validation.json) and [timings](bridgeless/sweep_timings.csv).',
        '- [Euclidean baseline](euclidean/samples.csv): 220 distinct observed clouds. Reusing clean clouds across models gives 240 model/amplitude/seed rows. The baseline is also shared across graph methods. See [error quartiles](euclidean/summary.csv) and [paired comparisons](euclidean/paired_summary.csv).',
        '- [Grid oracle](grid_oracle_samples.csv): the lowest error among connected points on the prescribed 24-budget grid, chosen separately for each of 720 cases. Extra diagnostic budgets are excluded. See [group summaries](grid_oracle_summary.csv) and [selector gaps](selected_vs_oracle.csv). This benchmark uses the clean reference.',
        '- [All amplitudes](all_amplitudes.csv): both rules, all six amplitudes, error medians/quartiles and median normalized budgets. Samples and displacements remain paired.',
        '- [Sample counts](rank_display_counts.json): the primary display has 5,000 shared seed/size combinations, 15,000 profile/seed samples and 30,000 method records. The full release has 30,780 method records.',
        '- [Literature sources](literature_checks.json): source details for the connectivity, spectral-distance and cyclic-graph comparisons in the paper.',
        '',
        '## Connected bridge-free graphs',
        '',
        'The test inserts a whole integer-k layer, then checks that the graph is connected and has no bridge. A bridge is an edge whose removal disconnects the graph. Every record satisfies k_conn ≤ k_br ≤ k_close. The bridge-free test can still stop before closure.',
        '',
        'At n=512, each profile has 1,500 paired seeds. Closure rates at connection and at the first connected bridge-free layer are:',
        '']
    display = {'uniform':'Uniform', 'strong_cosine':'Cosine', 'localized_trough':'Trough',
               'OR_kNN':'OR kNN', 'mutual_kNN':'mutual kNN', 'radial':'Radial',
               'isotropic':'Isotropic', 'R2_global':'radius'}
    for profile in ('uniform', 'strong_cosine', 'localized_trough'):
        for method in ('OR_kNN', 'mutual_kNN'):
            row = sweep[(sweep.n==512)&(sweep.profile==profile)&(sweep.method==method)].iloc[0]
            lines.append(f'- {display[profile]}, {display[method]}: {100*row.connection_rate:.1f}% at connection; {100*row.bridgeless_rate:.1f}% at the bridge-free layer.')
    lines += [
        '',
        'The CSV summaries give pointwise 95% Wilson intervals. Profiles and rules reuse seeds, so pooled rows are not independent repetitions.',
        '',
        '## Noise comparison at delta=0.01',
        '',
        'The clean Euclidean median error is about 0.1154 against the intrinsic reference. To compare noise sensitivity, subtract each pipeline’s own clean error from its noisy error for each seed before computing medians and paired counts. The longest-lifetime rule selects its budget separately at each amplitude. Euclidean increments are shared across graph methods.',
        '']
    for row in oracle[oracle.delta==.01].itertuples():
        baseline = euclidean[(euclidean.model==row.model)&(euclidean.delta==row.delta)].iloc[0]
        choice = paired[(paired.model==row.model)&(paired.delta==row.delta)&
                        (paired.method==row.method)&(paired.selector=='longest_lifetime')].iloc[0]
        lines.append(f'- {display[row.model]}, {display[row.method]}: median noise increment {choice.median_selected_noise_increment:.5f} for longest-lifetime selection and {choice.median_euclidean_noise_increment:.5f} for Euclidean; selection has the smaller increment in {int(choice.increment_selected_better)}/{int(choice.samples)} paired samples. Grid-oracle absolute error: {row.median_e_PH:.5f}.')
    lines += [
        '',
        'These results describe the tested models and grid. The selection rules have no general guarantee under noise. The grid oracle uses the clean reference and serves as a benchmark.',
        '',
        '## Reproduction',
        '',
        'Run commands from the release folder with the environment and complete sequence in the [code guide](../../code/README.md). To rebuild these summaries and this guide from saved measurements, run:',
        '',
        '```sh',
        'python3 code/reproduce_results.py',
        'python3 code/comparison_report.py --data-only',
        'python3 code/verify_comparisons.py',
        'python3 tools/record_current.py --data-only',
        'python3 code/verify_release.py --data-only',
        '```',
        '',
        'To rerun the bridge comparison and Euclidean baseline, use:',
        '',
        '```sh',
        'python3 code/bridge_free.py --validate --clean --sweep',
        'python3 code/euclidean_baseline.py',
        'python3 code/reproduce_results.py',
        'python3 code/comparison_report.py --data-only',
        'python3 code/verify_comparisons.py',
        'python3 tools/record_current.py --data-only',
        'python3 code/verify_release.py --data-only',
        '```',
        '',
        'The bridge calculation requires g++ with C++17 support. The helper is compiled in a temporary folder. Verification uses saved evidence and needs no C++ compiler. The [code guide](../../code/README.md) explains checkpoint reuse and fresh runs.',
        '',
        '## Timing definitions',
        '',
        'Individual Euclidean times include distances and persistence, excluding checkpoint file I/O. The per-checkpoint times cover all 220 diagrams.',
        '',
        'Bridge times include sampling, helper startup, graph reconstruction and event calculations, excluding saved-row writing. Child RSS is the cumulative maximum across helper processes. The pilot also includes compiler memory when compilation was needed.',
        '',
        'Timing records describe the runs that produced the saved measurements and depend on the computer and its load. Rebuilding this guide does not measure time again.',
        '']
    (OUT/'README.md').write_text('\n'.join(lines).replace('`', chr(96)))
    record = dict(
        status='complete', scope='Bridge-free study (two n=512 proportions are cited in the paper), Example 1 certificate and the additional comparisons; this summary step reads saved measurements and does not regenerate them.',
        largest_bridgeless_n=int(timings.n.max()),
        bridgeless_profile_seed_samples=len(timings), bridgeless_method_records=int(sweep.samples.sum()),
        clean_samples=60, distinct_euclidean_diagrams=220,
        euclidean_model_amplitude_seed_rows=240, grid_oracle_samples=len(raw_oracle),
        cyclic_orders=separation['cyclic_orders_checked'],
        distinct_threshold_graphs=separation['threshold_graphs'],
        summary_environment=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__))
    paper=None if data_only else locate_manuscript(F,manuscript_dir)
    previous_path=OUT/'completion.json'
    if paper is None and previous_path.is_file():
        previous=json.loads(previous_path.read_text())
        for key in ('manuscript_source_sha256','manuscript_pdf_sha256'):
            if key in previous:record[key]=previous[key]
    if paper is not None:
        record['manuscript_source_sha256']=hashlib.sha256((paper/'combined.tex').read_bytes()).hexdigest()
        record['manuscript_pdf_sha256']=hashlib.sha256((paper/'combined.pdf').read_bytes()).hexdigest()
    (OUT/'completion.json').write_text(json.dumps(record, indent=2)+'\n')
    print('Study-comparison guide and evidence record written.')

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--manuscript-dir',type=Path)
    group.add_argument('--data-only',action='store_true',help='Rebuild the comparison guide from release records only.')
    args=parser.parse_args()
    main(args.manuscript_dir,args.data_only)
