"""Locate an optional manuscript kept outside the numerical release."""
from pathlib import Path
import re

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


_INPUT = re.compile(r'\\input\{([^}]+)\}')


def _code_part(line):
    """Split a LaTeX line into its code and comment parts (unescaped %)."""
    match = re.search(r'(?<!\\)%', line)
    return (line, '') if match is None else (line[:match.start()], line[match.start():])


def manuscript_inputs(paper, name='combined.tex'):
    """Relative paths of every active \\input file, recursively, in source order."""
    paper = Path(paper)
    found = []
    def walk(path):
        for line in path.read_text().splitlines():
            code, _ = _code_part(line)
            for target in _INPUT.findall(code):
                rel = target if target.endswith('.tex') else target + '.tex'
                if rel not in found:
                    found.append(rel)
                    walk(paper/rel)
    walk(paper/name)
    return found


def read_manuscript(paper, name='combined.tex'):
    """Return the manuscript source with every active \\input expanded.

    The paper's own text, labels and references may live in \\input files; the
    structural and numerical checks read this expanded text. Commented-out
    \\input lines are left unexpanded."""
    paper = Path(paper)
    def expand(path, depth=0):
        if depth > 10:
            raise ValueError('Nested \\input deeper than 10 levels.')
        out = []
        for line in path.read_text().splitlines(keepends=True):
            code, comment = _code_part(line)
            if _INPUT.search(code) is None:
                out.append(line)
                continue
            def repl(match):
                rel = match.group(1) if match.group(1).endswith('.tex') else match.group(1) + '.tex'
                return expand(paper/rel, depth + 1)
            out.append(_INPUT.sub(repl, code) + comment)
        return ''.join(out)
    return expand(paper/name)


def strip_comments(text):
    """Remove LaTeX comments (unescaped % to end of line), keeping line breaks."""
    return ''.join(_code_part(line)[0] + ('\n' if line.endswith('\n') else '')
                   for line in text.splitlines(keepends=True))


def manuscript_text(paper, name='combined.tex'):
    """Expanded source without comments: the text the compiled paper is built from."""
    return strip_comments(read_manuscript(paper, name))
