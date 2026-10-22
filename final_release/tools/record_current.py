"""Record the manuscript identity and structure, output counts and the hash of every release file."""
from pathlib import Path
import argparse, hashlib, json, os, re, sys
from external_manuscript import locate_manuscript, manuscript_text, manuscript_inputs
sys.dont_write_bytecode = True
from pypdf import PdfReader

F = Path(__file__).resolve().parents[1]
ENVIRONMENTS = ('definition', 'theorem', 'proposition', 'corollary',
                'structurallemma', 'appendixlemma', 'lemma', 'remark', 'example', 'proof')
DISPLAY_ENVIRONMENTS = ('equation', 'equation*', 'align', 'align*',
                        'gather', 'gather*', 'multline', 'multline*', 'displaymath')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def result_numbers(text):
    numbers = {}
    for name in ('theorem', 'proposition', 'corollary', 'lemma', 'example', 'remark'):
        blocks = re.findall(r'\\begin\{'+name+r'\}([\s\S]*?)\\end\{'+name+r'\}', text)
        numbers[name] = {
            label:i for i, block in enumerate(blocks, 1)
            for label in re.findall(r'\\label\{([^}]+)\}', block)}
    return numbers

def manuscript_records(root,manuscript_dir=None):
    paper = locate_manuscript(root,manuscript_dir)
    if paper is None:
        raise FileNotFoundError('Supply --manuscript-dir to compare an external manuscript.')
    for name in ('combined.tex', 'combined.pdf'):
        if not (paper/name).is_file():
            raise FileNotFoundError(f'Place {name} in the external manuscript directory before recording it.')
    # Expanded source (active \input files inlined, comments removed): the compiled text.
    text = manuscript_text(paper)
    labels = re.findall(r'\\label\{([^}]+)\}', text)
    references = re.findall(r'\\(?:ref|eqref)\{([^}]+)\}', text)
    citations = {k for group in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}', text) for k in group.split(',')}
    bibliography = set(re.findall(r'\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}', text))
    if len(labels) != len(set(labels)):
        raise ValueError('The manuscript has duplicate labels.')
    if set(references) - set(labels):
        raise ValueError('The manuscript has unresolved source references.')
    if citations - bibliography:
        raise ValueError('The manuscript has unresolved source citations.')
    source = dict(file='combined.tex', sha256=sha(paper/'combined.tex'))
    pdf = dict(file='combined.pdf', sha256=sha(paper/'combined.pdf'),
               pages=len(PdfReader(paper/'combined.pdf').pages))
    source_checks = dict(unique_labels=len(labels), resolved_citation_keys=len(citations),
                         unresolved_references=0, unresolved_citations=0)
    inputs = [dict(file=rel, sha256=sha(paper/rel)) for rel in manuscript_inputs(paper)]
    graphics = [dict(file=g, sha256=sha(paper/g))
                for g in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', text)]
    identity = dict(scope='Accompanying manuscript.', path_base='manuscript_directory', source=source, pdf=pdf,
                    inputs=inputs, graphics=graphics, source_checks=source_checks)
    environments = {name:text.count(r'\begin{'+name+'}') for name in ENVIRONMENTS}
    displays = {'bracketed':text.count(r'\[')}
    displays.update({name:text.count(r'\begin{'+name+'}') for name in DISPLAY_ENVIRONMENTS
                     if r'\begin{'+name+'}' in text})
    numbers = result_numbers(text)
    layout = dict(
        status='passed',
        scope='Source and PDF identities, labels, citations, environment counts and result numbering.',
        source_sha256=source['sha256'], pdf_sha256=pdf['sha256'],
        pdf_pages=pdf['pages'], source_checks=source_checks,
        environment_counts=environments, display_math_blocks_by_type=displays,
        result_numbers=numbers, figure_labels=re.findall(r'\\label\{(fig:[^}]+)\}', text),
        table_labels=re.findall(r'\\label\{(tab:[^}]+)\}', text),
        appendix_labels=[l for l in re.findall(r'\\label\{(app:[^}]+)\}', text)])
    return text, identity, layout

def main(manuscript_dir=None,data_only=False):
    results = F/'results'
    paper=None if data_only else locate_manuscript(F,manuscript_dir)
    if paper is not None:
        _, identity, layout = manuscript_records(F,paper)
        (results/'manuscript_record.json').write_text(json.dumps(identity, indent=2)+'\n')
        (results/'manuscript_structure.json').write_text(json.dumps(layout, indent=2)+'\n')
        completion_path = results/'comparisons/completion.json'
        if completion_path.exists():
            completion = json.loads(completion_path.read_text())
            completion['manuscript_source_sha256'] = identity['source']['sha256']
            completion['manuscript_pdf_sha256'] = identity['pdf']['sha256']
            completion_path.write_text(json.dumps(completion, indent=2)+'\n')
    import csv
    def row_count(name):
        with (results/name).open() as stream:return sum(1 for _ in csv.DictReader(stream))
    generation=dict(
        scope='Complete release: evidence for the paper plus the additional studies. This recorder reads outputs and does not generate measurements.',
        retention_table_cells=row_count('rank_table.csv'),
        retention_table_method_records=sum(int(row['samples']) for row in csv.DictReader((results/'rank_table.csv').open())),
        practical_choice_cells=row_count('practical_choices/choices.csv'),
        peanut_neck_method_records=row_count('peanut_neck/trials.csv'),
        small_examples_report='results/small_examples_verification.json',
        practical_claims_report='results/practical_claims_verification.json',
        extra_layers_peanut_report='results/extra_layers_peanut_verification.json',
        extra_layer_cells=row_count('retention_extra_layers/retention.csv'),
        peanut_trough_clouds=len(list((results/'peanut_trough/checkpoints').glob('*.npz'))),
        peanut_trough_method_records=row_count('peanut_trough/closure.csv'),
        manuscript_evidence=[
            'Tables 1-2: established ingredients and early notation; written in manuscript.',
            'Table 3: retention_compact.csv; 24 primary cells from saved ranks and bridge-free records.',
            'Figures 1-3: closure_schematic.png, arc_support_fan.png, retention_extra_layers.png; the manuscript uses only these PNGs.',
            'Three theorems: absence before closure, arc-supported module, and circle events.',
            'Corollaries 1-2: circle error ceiling and small-gap closure optimality.',
            'Propositions 1-2: uniform-circle radius enlargement and smooth-curve locality.',
            'Examples 1-2 and Remarks 1-2: exact or rational finite certificates; proofs in paper.',
            'Sections 7-8: clean event diagrams, fixed-k and bridge-free counts, extra layers, diagram accuracy and the explicit independent-key heuristic.',
            'Appendix A: two-curve numerical locality check; Appendix B: ties; Appendix C: seven-point counterexample.',
            'Additional studies, outside the focused paper: noise, selectors, CkNN, ellipse sweeps, peanut/trough and varying-waist neck studies.'
        ],
        protocol_sha256=sha(F/'code/config.json'),
        clean_samples=len(list((results/'clean_closure/checkpoints').glob('*.json'))),
        clean_budget_evaluations=row_count('clean_closure/cells.csv'),
        clean_closed_prefix_diagrams=row_count('useful_diagram/saved_prefixes.csv'),
        rank_method_records=row_count('extra_k/rank_sweep/samples.csv'),
        bridge_method_records=row_count('comparisons/bridgeless/sweep_samples.csv'),
        matched_radial_samples=len(list((results/'matched_noise/checkpoints').glob('*.json'))),
        ellipse_coordinate_sets=len(list((results/'ellipse_checks/checkpoints').glob('*.json'))),
        ellipse_sweep_samples=len(list((results/'ellipse/checkpoints').glob('*.json'))),
        distinct_noise_clouds=len(list((results/'noise/checkpoints').glob('*.npz'))),
        distinct_euclidean_diagrams=len(list((results/'comparisons/euclidean/checkpoints').glob('*.npz'))),
        noise_budget_evaluations=row_count('noise/cells.csv'),
        noise_selector_records=row_count('noise/selectors.csv'),
        manuscript_png_figures=len(list((F/'manuscript_figures').glob('*.png'))),
        figure_pdf_png_pairs=sum(1 for folder in ('results/figure_exports','results/additional_figures') for path in (F/folder).glob('*.pdf') if path.with_suffix('.png').is_file() or (folder=='results/figure_exports' and (F/'manuscript_figures'/path.with_suffix('.png').name).is_file())),
        generation_commands='code/reproduce.py; individual study generators are listed in code/README.md',
        clock_times='Measured again on each run; not fixed numerical targets.')
    (results/'generation_summary.json').write_text(json.dumps(generation,indent=2)+'\n')
    excluded = {results/'file_manifest.json', results/'verification_report.json'}
    # Logs can still grow while the recorder and verifier run. They are
    # diagnostic output, outside the hash-checked scientific payload.
    logs = results/'run_logs'
    paths = sorted(p for p in F.rglob('*') if p.is_file() and p not in excluded
                   and logs not in p.parents and '__pycache__' not in p.parts and p.suffix!='.pyc')
    manifest = {str(p.relative_to(F)):sha(p) for p in paths}
    (results/'file_manifest.json').write_text(json.dumps({
        'algorithm':'SHA-256',
        'excluded':['results/file_manifest.json', 'results/verification_report.json', 'results/run_logs/'],
        'sha256':manifest}, indent=2)+'\n')
    print(f"Recorded {len(manifest)} release files; external manuscript comparison: {paper is not None}.")

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--manuscript-dir',type=Path)
    group.add_argument('--data-only',action='store_true',help='Record the release without reading a manuscript; the saved manuscript records are kept.')
    args=parser.parse_args()
    main(args.manuscript_dir,args.data_only)
