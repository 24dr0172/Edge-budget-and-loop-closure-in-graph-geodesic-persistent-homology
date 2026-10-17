"""Record the accompanying manuscript identity and current release file hashes."""
from pathlib import Path
import hashlib,json,re,sys
sys.dont_write_bytecode=True
from pypdf import PdfReader
F=Path(__file__).resolve().parents[1];PAPER=F.parent;R=F/'results'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    text=(PAPER/'combined.tex').read_text()
    labels=re.findall(r'\\label\{([^}]+)\}',text)
    references=re.findall(r'\\(?:ref|eqref)\{([^}]+)\}',text)
    citations={k for group in re.findall(r'\\cite\w*\{([^}]+)\}',text) for k in group.split(',')}
    bibliography=set(re.findall(r'\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}',text))
    assert len(labels)==len(set(labels)) and not(set(references)-set(labels))
    assert not(citations-bibliography)
    source=dict(file='../../combined.tex',sha256=sha(PAPER/'combined.tex'))
    pdf=dict(file='../../combined.pdf',sha256=sha(PAPER/'combined.pdf'),pages=len(PdfReader(PAPER/'combined.pdf').pages))
    identity=dict(scope='Edge budget and loop closure; accompanying manuscript and unified local release.',source=source,pdf=pdf,
        source_checks=dict(unique_labels=len(labels),resolved_citation_keys=len(citations),unresolved_references=0,unresolved_citations=0))
    (R/'manuscript_revision.json').write_text(json.dumps(identity,indent=2)+'\n')
    inventory_path=R/'release_inventory.json'
    if inventory_path.exists():
        inventory=json.loads(inventory_path.read_text())
        for key in ['payload_files_excluding_inventory_and_manifest','payload_bytes_excluding_inventory_and_manifest']:
            inventory.pop(key,None)
        inventory['date']='2026-10-16'
        inventory['appendix_figures']=text.split(r'\appendix',1)[1].count(r'\begin{figure}')
        inventory['figure_pdf_png_pairs']=len(list((F/'manuscript_figures').glob('*.pdf')))
        inventory['expanded_noise_distinct_clouds']=220
        inventory['phase1_completion']='results/phase1/completion.json'
        inventory['payload_exclusions']=['results/release_inventory.json','results/file_manifest.json','results/verification_report.json']
        payload=[p for p in F.rglob('*') if p.is_file() and str(p.relative_to(F)) not in inventory['payload_exclusions']]
        inventory['payload_files']=len(payload)
        inventory['payload_bytes']=sum(p.stat().st_size for p in payload)
        inventory_path.write_text(json.dumps(inventory,indent=2)+chr(10))
    excluded={R/'file_manifest.json',R/'verification_report.json'}
    paths=sorted(p for p in F.rglob('*') if p.is_file() and p not in excluded)
    manifest={str(p.relative_to(F)):sha(p) for p in paths}
    (R/'file_manifest.json').write_text(json.dumps({'algorithm':'SHA-256','excluded':['results/file_manifest.json','results/verification_report.json'],'sha256':manifest},indent=2)+'\n')
    print(f"Recorded manuscript identity ({pdf['pages']} pages) and {len(manifest)} release files.")
if __name__=='__main__':main()
