"""Locate an optional manuscript kept outside the numerical release."""
from pathlib import Path

def locate_manuscript(root,directory=None):
    root=Path(root).resolve()
    paper=Path(directory).expanduser().resolve() if directory is not None else root.parent
    if paper==root or root in paper.parents:
        raise ValueError('Keep the manuscript outside the release; pass its directory with --manuscript-dir.')
    if all((paper/name).is_file() for name in ('combined.tex','combined.pdf')):
        return paper
    if directory is not None:
        raise FileNotFoundError('The manuscript directory must contain combined.tex and combined.pdf.')
    return None
