"""Prepare an external manuscript and generated figures for optional local LaTeX use."""
from pathlib import Path
import argparse,shutil,subprocess,tempfile
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True,help='External combined.tex source.')
    p.add_argument('--output',type=Path,help='New empty compilation folder; default is a temporary folder.')
    p.add_argument('--compile',action='store_true',help='Run an already installed pdflatex twice.')
    a=p.parse_args();out=a.output.resolve() if a.output else Path(tempfile.mkdtemp(prefix='edge-budget-paper-'))
    if out==ROOT or ROOT in out.parents:p.error('Keep manuscript builds outside the release.')
    if a.source.resolve()==ROOT or ROOT in a.source.resolve().parents:p.error('Supply an external source file.')
    if out.exists() and any(out.iterdir()):p.error('Choose an empty output folder.')
    out.mkdir(parents=True,exist_ok=True)
    shutil.copy2(a.source,out/'combined.tex')
    shutil.copytree(ROOT/'manuscript_figures',out/'github release/manuscript_figures')
    if a.compile:
        for _ in range(2):subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error','combined.tex'],cwd=out,check=True)
    print(out)
if __name__=='__main__':main()
