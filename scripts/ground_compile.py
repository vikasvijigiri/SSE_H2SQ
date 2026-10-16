# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Compile the local Overleaf validation report; keep auxiliary files ignored."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('paper',type=Path);args=parser.parse_args()
    paper=args.paper.resolve();assert paper.is_relative_to(ROOT/'docs/manuscript/benchmark')
    assert (paper/'main.tex').is_file()
    output=ROOT/'.cache/latex-ground'/paper.name;output.mkdir(parents=True,exist_ok=True)
    commands=[['pdflatex.exe','--disable-installer','-no-shell-escape','-interaction=nonstopmode','-halt-on-error',f'-output-directory={output}','main.tex'],
        ['bibtex.exe',str(output/'main')]]
    logs=[]
    for command in [commands[0],commands[1],commands[0],commands[0]]:
        run=subprocess.run(command,cwd=paper,capture_output=True,text=True,errors='replace',timeout=180)
        logs.append(run.stdout+run.stderr)
        if run.returncode:
            (output/'build-output.log').write_text('\n'.join(logs),encoding='utf-8')
            raise SystemExit(f'LaTeX build failed; inspect {output}.\n'+logs[-1][-2000:])
    log=(output/'main.log').read_text(errors='replace')
    if re.search(r'undefined references|Citation .* undefined|Reference .* undefined',log):raise SystemExit('Unresolved citations or references remain.')
    if 'Overfull \\hbox' in log:raise SystemExit('An overfull horizontal box needs layout review.')
    shutil.copyfile(output/'main.pdf',paper/'validation-report.pdf')
    result=dict(compiled=True,engine='pdfLaTeX + BibTeX + two resolving passes; automatic installation and shell escape disabled',
        pages=int(re.search(r'Output written.*?\((\d+) pages?',log,re.S).group(1)),
        pdf_sha256=hashlib.sha256((paper/'validation-report.pdf').read_bytes()).hexdigest(),auxiliary_directory=str(output))
    (paper/'metadata/latex-build.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    manifest={p.relative_to(paper).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paper.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    (paper/'metadata/manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
