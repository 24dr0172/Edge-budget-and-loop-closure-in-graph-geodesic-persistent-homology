"""Record the current manuscript, structure checks and release file hashes."""
from datetime import date
from pathlib import Path
import argparse, hashlib, json, os, re, sys
from external_manuscript import locate_manuscript
sys.dont_write_bytecode = True
from pypdf import PdfReader

F = Path(__file__).resolve().parents[1]
ENVIRONMENTS = ('definition', 'theorem', 'proposition', 'corollary',
                'structurallemma', 'appendixlemma', 'remark', 'example', 'proof')
DISPLAY_ENVIRONMENTS = ('equation', 'equation*', 'align', 'align*',
                        'gather', 'gather*', 'multline', 'multline*', 'displaymath')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def manuscript_records(root,manuscript_dir=None):
    paper = locate_manuscript(root,manuscript_dir)
    if paper is None:
        raise FileNotFoundError('Supply --manuscript-dir to compare an external manuscript.')
    for name in ('combined.tex', 'combined.pdf'):
        if not (paper/name).is_file():
            raise FileNotFoundError(f'Place {name} in the external manuscript directory before recording it.')
    text = (paper/'combined.tex').read_text()
    labels = re.findall(r'\\label\{([^}]+)\}', text)
    references = re.findall(r'\\(?:ref|eqref)\{([^}]+)\}', text)
    citations = {k for group in re.findall(r'\\cite\w*\{([^}]+)\}', text) for k in group.split(',')}
    bibliography = set(re.findall(r'\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}', text))
    if len(labels) != len(set(labels)):
        raise ValueError('The manuscript has duplicate labels.')
    if set(references) - set(labels):
        raise ValueError('The manuscript has unresolved source references.')
    if citations - bibliography:
        raise ValueError('The manuscript has unresolved source citations.')
    source = dict(file=os.path.relpath(paper/'combined.tex',root/'results'), sha256=sha(paper/'combined.tex'))
    pdf = dict(file=os.path.relpath(paper/'combined.pdf',root/'results'), sha256=sha(paper/'combined.pdf'),
               pages=len(PdfReader(paper/'combined.pdf').pages))
    source_checks = dict(unique_labels=len(labels), resolved_citation_keys=len(citations),
                         unresolved_references=0, unresolved_citations=0)
    identity = dict(scope='Current accompanying manuscript.', source=source, pdf=pdf,
                    source_checks=source_checks)
    environments = {name:text.count(r'\begin{'+name+'}') for name in ENVIRONMENTS}
    displays = {'bracketed':text.count(r'\[')}
    displays.update({name:text.count(r'\begin{'+name+'}') for name in DISPLAY_ENVIRONMENTS
                     if r'\begin{'+name+'}' in text})
    result_numbers = {}
    for name in ('theorem', 'proposition', 'corollary', 'example', 'remark'):
        blocks = re.findall(r'\\begin\{'+name+r'\}([\s\S]*?)\\end\{'+name+r'\}', text)
        result_numbers[name] = {
            label:i for i, block in enumerate(blocks, 1)
            for label in re.findall(r'\\label\{([^}]+)\}', block)}
    layout = dict(
        status='passed',
        scope='Current source/PDF identities, labels, citations, environment counts and result numbering.',
        current_source_sha256=source['sha256'], current_pdf_sha256=pdf['sha256'],
        pdf_pages=pdf['pages'], source_checks=source_checks,
        environment_counts=environments, display_math_blocks_by_type=displays,
        result_numbers=result_numbers)
    return text, identity, layout

def main(manuscript_dir=None):
    results = F/'results'
    paper=locate_manuscript(F,manuscript_dir)
    text=''
    if paper is not None:
        text, identity, layout = manuscript_records(F,paper)
        (results/'manuscript_revision.json').write_text(json.dumps(identity, indent=2)+'\n')
        (results/'manuscript_restructure.json').write_text(json.dumps(layout, indent=2)+'\n')
        completion_path = results/'phase1/completion.json'
        if completion_path.exists():
            completion = json.loads(completion_path.read_text())
            completion['manuscript_source_sha256'] = identity['source']['sha256']
            completion['manuscript_pdf_sha256'] = identity['pdf']['sha256']
            completion_path.write_text(json.dumps(completion, indent=2)+'\n')
    inventory_path = results/'release_inventory.json'
    if inventory_path.exists():
        inventory = json.loads(inventory_path.read_text())
        for key in ('payload_files_excluding_inventory_and_manifest',
                    'payload_bytes_excluding_inventory_and_manifest', 'phase1_completion'):
            inventory.pop(key, None)
        inventory['date'] = date.today().isoformat()
        if text:inventory['appendix_figures'] = text.split(r'\appendix', 1)[1].count(r'\begin{figure}')
        inventory['figure_pdf_png_pairs'] = len(list((F/'manuscript_figures').glob('*.pdf')))
        inventory['expanded_noise_distinct_clouds'] = 220
        inventory['comparison_record'] = 'results/phase1/completion.json'
        excluded = ['results/release_inventory.json', 'results/file_manifest.json',
                    'results/verification_report.json']
        inventory['payload_exclusions'] = excluded
        payload = [p for p in F.rglob('*') if p.is_file() and str(p.relative_to(F)) not in excluded]
        inventory['payload_files'] = len(payload)
        inventory['payload_bytes'] = sum(p.stat().st_size for p in payload)
        inventory_path.write_text(json.dumps(inventory, indent=2)+'\n')
    import csv
    def row_count(name):
        with (results/name).open() as stream:return sum(1 for _ in csv.DictReader(stream))
    generation=dict(
        scope='Current seed protocol and generated scientific outputs; run verify_release.py to check their agreement with the paper.',
        protocol_sha256=sha(F/'code/config.json'),
        clean_samples=len(list((results/'clean_closure/checkpoints').glob('*.json'))),
        clean_budget_evaluations=row_count('clean_closure/cells.csv'),
        clean_closed_prefix_diagrams=row_count('useful_diagram/saved_prefixes.csv'),
        rank_method_records=row_count('extra_k/candidate_audit/samples.csv'),
        bridge_method_records=row_count('phase1/bridgeless/sweep_samples.csv'),
        matched_radial_samples=len(list((results/'matched_noise/checkpoints').glob('*.json'))),
        ellipse_coordinate_sets=len(list((results/'ellipse_checks/checkpoints').glob('*.json'))),
        ellipse_sweep_samples=len(list((results/'ellipse/checkpoints').glob('*.json'))),
        distinct_noise_clouds=len(list((results/'noise/checkpoints').glob('*.npz'))),
        distinct_euclidean_diagrams=len(list((results/'phase1/euclidean/checkpoints').glob('*.npz'))),
        noise_budget_evaluations=row_count('noise/cells.csv'),
        noise_selector_records=row_count('noise/selectors.csv'),
        figure_pdf_png_pairs=sum(1 for path in (F/'manuscript_figures').glob('*.pdf') if path.with_suffix('.png').is_file()),
        generation_commands='code/reproduce.py; individual study generators are listed in code/README.md',
        clock_times='Measured again on each run; not fixed numerical targets.')
    (results/'generation_summary.json').write_text(json.dumps(generation,indent=2)+'\n')
    excluded = {results/'file_manifest.json', results/'verification_report.json'}
    paths = sorted(p for p in F.rglob('*') if p.is_file() and p not in excluded and '__pycache__' not in p.parts and p.suffix!='.pyc')
    manifest = {str(p.relative_to(F)):sha(p) for p in paths}
    (results/'file_manifest.json').write_text(json.dumps({
        'algorithm':'SHA-256',
        'excluded':['results/file_manifest.json', 'results/verification_report.json'],
        'sha256':manifest}, indent=2)+'\n')
    print(f"Recorded {len(manifest)} release files; external manuscript comparison: {paper is not None}.")

if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--manuscript-dir',type=Path)
    main(parser.parse_args().manuscript_dir)
