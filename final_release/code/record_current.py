"""Record the current manuscript, structure checks and release file hashes."""
from datetime import date
from pathlib import Path
import hashlib, json, re, sys
sys.dont_write_bytecode = True
from pypdf import PdfReader

F = Path(__file__).resolve().parents[1]
ENVIRONMENTS = ('definition', 'theorem', 'proposition', 'corollary',
                'structurallemma', 'appendixlemma', 'remark', 'example', 'proof')
DISPLAY_ENVIRONMENTS = ('equation', 'equation*', 'align', 'align*',
                        'gather', 'gather*', 'multline', 'multline*', 'displaymath')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def manuscript_records(root):
    paper = root.parent
    for name in ('combined.tex', 'combined.pdf'):
        if not (paper/name).is_file():
            raise FileNotFoundError(f'Place {name} beside final_release before recording the manuscript.')
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
    source = dict(file='../../combined.tex', sha256=sha(paper/'combined.tex'))
    pdf = dict(file='../../combined.pdf', sha256=sha(paper/'combined.pdf'),
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

def main():
    results = F/'results'
    text, identity, layout = manuscript_records(F)
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
        inventory['appendix_figures'] = text.split(r'\appendix', 1)[1].count(r'\begin{figure}')
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
    excluded = {results/'file_manifest.json', results/'verification_report.json'}
    paths = sorted(p for p in F.rglob('*') if p.is_file() and p not in excluded)
    manifest = {str(p.relative_to(F)):sha(p) for p in paths}
    (results/'file_manifest.json').write_text(json.dumps({
        'algorithm':'SHA-256',
        'excluded':['results/file_manifest.json', 'results/verification_report.json'],
        'sha256':manifest}, indent=2)+'\n')
    print(f"Recorded current manuscript ({identity['pdf']['pages']} pages) and {len(manifest)} release files.")

if __name__ == '__main__':
    main()
