"""Write the current Phase 1 completion record from its saved outputs."""
from pathlib import Path
import hashlib,json,platform,sys
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
F=Path(__file__).resolve().parents[1];OUT=F/'results/phase1'


def main():
    separation=json.loads((F/'results/separation_verification.json').read_text())['filtration']
    sweep=pd.read_csv(OUT/'bridgeless/sweep_summary.csv')
    timings=pd.read_csv(OUT/'bridgeless/sweep_timings.csv')
    euclidean=pd.read_csv(OUT/'euclidean/summary.csv')
    paired=pd.read_csv(OUT/'euclidean/paired_summary.csv')
    oracle=pd.read_csv(OUT/'grid_oracle_summary.csv')
    raw_oracle=pd.read_csv(OUT/'grid_oracle_samples.csv')
    lines=[
        '# Phase 1: completed computations and checks',
        '',
        'All nine Phase 1 tasks are complete. The records below supply the computational evidence for the revision plan. Manuscript restructuring and the integration of these results are tracked separately in that plan.',
        '',
        '- **1.1, Example 1:** all 120 cyclic orders were checked at all 14 distinct threshold graphs, including the empty graph. Exact radical arithmetic and outward integer bounds identify ties and certify distance comparisons. Each graph has an admissible order separately; the intersection over thresholds is empty. Direct enumeration independently checks all orientations for all 120 orders at each of the two proof graphs. See [the certificate](../separation_verification.json) and [the exact argument](../arc_support_example/exact_calculation.md).',
        '- **1.2, original clean study:** all 60 samples, both graph rules, were replayed. Sparse event results agree with full Cartesian ranks and independent depth-first bridge traversal. See [raw events](bridgeless/clean_samples.csv) and [summary and Wilson intervals](bridgeless/clean_summary.csv).',
        '- **1.3, full planned sweep:** 15,000 profile/seed/size samples and 30,000 method records through n = 32,768 are complete, using every corresponding released seed. See [raw events](bridgeless/sweep_samples.csv), [summary and intervals](bridgeless/sweep_summary.csv), [validation](bridgeless/validation.json), and [timings](bridgeless/sweep_timings.csv).',
        '- **1.4, Euclidean baseline:** full unrestricted finite H1 diagrams were computed for all 220 distinct observed clouds. Reusing the shared clean cloud gives 240 model/amplitude/seed rows, with each baseline shared across graph methods. See [baseline errors](euclidean/samples.csv), [quantiles](euclidean/summary.csv), and [paired rule comparisons](euclidean/paired_summary.csv).',
        '- **1.5, grid oracle:** 720 sample-specific minima are computed over connected points of the prescribed 24-budget grid; extra diagnostic budgets are excluded. Minimize within each sample before taking group medians or quartiles. See [sample choices](grid_oracle_samples.csv), [36 group summaries](grid_oracle_summary.csv), and [paired selector gaps](selected_vs_oracle.csv). This comparator uses the clean reference.',
        '- **1.6, amplitudes:** [all_amplitudes.csv](all_amplitudes.csv) contains both rules at all six amplitudes, separately by noise model and method, including medians, IQRs and median normalized budgets. Shared clean samples and displacement vectors remain paired.',
        '- **1.7, displayed sweep counts:** [rank_display_counts.json](rank_display_counts.json) confirms 5,000 seed/size combinations, 15,000 profile/seed samples and 30,000 method records in the primary display. The full released rank study retains all 30,780 method records, including the two larger sizes.',
        '- **1.8, connectivity wording:** retain “a simple heuristic”. The bounded primary-source check does not establish that choosing the smallest connecting parameter is standard practice. The checked [Isomap documentation](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.Isomap.html) uses a fixed default; that does not support the stronger claim.',
        '- **1.9, spectral graph:** the [Damrich et al. primary paper](https://papers.nips.cc/paper_files/paper/2024/file/4a32a646254d2e37fc74a38d65796552-Paper-Conference.pdf) specifies unweighted symmetric (OR) kNN graphs in the reported spectral comparisons. Keep the qualified robustness comparison and its parameter dependence. Source locations and the wording decision are in [literature_checks.json](literature_checks.json).',
        '',
        '## Connected bridgeless events',
        '',
        'Every complete integer-k layer is inserted before testing. The graph must be connected and have no bridges; minimum degree two is not substituted. For every record, k_conn ≤ k_br ≤ k_close. Closure at k_br is a necessary-condition comparison, not a loop guarantee.',
        '',
        'At n = 512 in the 1,500-seed-per-profile sweep, the closure percentages at connection and at the first connected bridgeless layer are:',
        ''
    ]
    display={'uniform':'Uniform','strong_cosine':'Cosine','localized_trough':'Trough',
             'OR_kNN':'OR','mutual_kNN':'mutual'}
    for row in sweep[sweep.n==512].itertuples():
        lines.append(f'- {display[row.profile]}, {display[row.method]}: {100*row.connection_rate:.1f}% at connection; {100*row.bridgeless_rate:.1f}% at the connected bridgeless layer.')
    lines+=['','All rates and pointwise 95% Wilson intervals are in the CSV summaries. Profiles and rules reuse seeds; pooled rows are not independent repetitions. The bridgeless rule leaves a substantial fraction of samples unclosed.',
            '',
            '## Noise comparison at the largest amplitude',
            '',
            'The following are medians across the same 20 seeds within each group. Euclidean errors are shared across graph methods. A difference of medians is not a median paired improvement; the latter and per-sample win counts are stored separately.',
            '']
    for row in oracle[oracle.delta==.01].itertuples():
        e=euclidean[(euclidean.model==row.model)&(euclidean.delta==row.delta)].iloc[0]
        q=paired[(paired.model==row.model)&(paired.delta==row.delta)&(paired.method==row.method)&(paired.selector=='longest_lifetime')].iloc[0]
        lines.append(f'- {row.model}, {row.method}: Euclidean {e.median_e_PH:.5f}; longest-lifetime selection {q.median_selected_e_PH:.5f}; grid oracle {row.median_e_PH:.5f}.')
    lines+=['',
            'These comparisons describe the sampled noise models and grid. Neither a long bar nor the connected bridgeless check is a noisy-data theorem.',
            '',
            '## Reproduction',
            '',
            'Run from final_release with the documented Python environment. Recomputing the connected bridgeless study additionally requires a C++17 compiler (g++); the helper is compiled into the system temporary directory, and no binary is distributed.',
            '',
            chr(96)*3+'sh',
            'python3 code/verify_separation.py',
            'python3 code/phase1_bridgeless.py --validate --clean --pilot',
            'python3 code/phase1_bridgeless.py --sweep',
            'python3 code/phase1_noise.py',
            'python3 code/phase1_summaries.py',
            'python3 code/phase1_report.py',
            'python3 code/verify_phase1.py',
            'python3 code/record_current.py',
            'python3 code/verify_release.py',
            chr(96)*3,
            '',
            'The bridge and Euclidean scripts resume verified current checkpoints. To recompute an existing Euclidean diagram, remove its output checkpoint first; inputs remain in results/noise/checkpoints. Do not delete or regenerate the original study records. The regular verifier needs no C++ compiler and does not rerun the large persistence study.',
            '',
            'Individual Euclidean timings include distances and persistence, excluding checkpoint I/O. Its pilot and full timing records describe separate invocations; per-checkpoint elapsed times cover all 220 diagrams. Bridge timings include sampling, helper startup, graph reconstruction and the event calculation; saved-row writing is excluded. The child RSS field is a cumulative maximum across helper processes, not a per-case allocation; the pilot invocation also includes compiler memory if compilation was needed. Bridge timings were measured while the independent serial Euclidean job was running, so they are feasibility measurements rather than isolated benchmark claims.',
            '',
            'The existing Figure PDF/PNG pairs, manuscript and prior experimental data are preserved. No deposition, registration or manuscript release URL is part of this phase.',
            '']
    (OUT/'README.md').write_text('\n'.join(lines))
    status=dict(status='complete',tasks=[f'1.{i}' for i in range(1,10)],
                manuscript_changed_in_phase1=False,largest_bridgeless_n=int(timings.n.max()),
                bridgeless_profile_seed_samples=len(timings),bridgeless_method_records=30000,
                original_clean_samples=60,distinct_euclidean_diagrams=220,
                euclidean_model_amplitude_seed_rows=240,grid_oracle_samples=len(raw_oracle),
                cyclic_orders=separation['cyclic_orders_checked'],
                distinct_threshold_graphs=separation['threshold_graphs'],
                python=platform.python_version(),numpy=np.__version__,
                manuscript_source_sha256=hashlib.sha256((F.parent/'combined.tex').read_bytes()).hexdigest(),
                manuscript_pdf_sha256=hashlib.sha256((F.parent/'combined.pdf').read_bytes()).hexdigest(),
                original_study_inputs_preserved=True)
    (OUT/'completion.json').write_text(json.dumps(status,indent=2)+'\n')
    print('Phase 1 completion record and current results guide written.')


if __name__=='__main__':main()
