# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Static asset/citation check for the local PRB/Overleaf source tree."""
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[2]
PAPER=ROOT/"manuscript_latest"/"review"/"overleaf"


def main():
    seen=set();sources=[]
    def visit(path):
        path=path.resolve()
        if path in seen:return
        seen.add(path);sources.append(path)
        text=re.sub(r"(?<!\\)%[^\n]*","",path.read_text(encoding="utf-8"))
        for name in re.findall(r"\\input\{([^}]+)\}|\\include\{([^}]+)\}",text):
            value=next(part for part in name if part)
            target=PAPER/(value if Path(value).suffix else value+".tex")
            if not target.is_file():raise SystemExit(f"Missing TeX input: {target.relative_to(PAPER)}")
            visit(target)
        for name in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}",text):
            target=PAPER/name
            if target.suffix:
                candidates=[target]
            else:
                candidates=[target.with_suffix(ext) for ext in (".pdf",".png",".jpg",".jpeg",".eps")]
            if not any(p.is_file() for p in candidates):
                raise SystemExit(f"Missing figure asset: {name} (from {path.relative_to(PAPER)})")
    visit(PAPER/"main.tex")
    joined="\n".join(re.sub(r"(?<!\\)%[^\n]*","",p.read_text(encoding="utf-8")) for p in sources)
    bibmatch=re.search(r"\\bibliography\{([^}]+)\}",joined)
    if not bibmatch:raise SystemExit("No bibliography command found.")
    bib=PAPER/(bibmatch.group(1)+".bib")
    if not bib.is_file():raise SystemExit(f"Missing bibliography: {bib.name}")
    keys=set(re.findall(r"@\w+\s*\{\s*([^,\s]+)",bib.read_text(encoding="utf-8")))
    cites={key.strip() for group in re.findall(r"\\cite\w*\*?(?:\[[^]]*\])*\{([^}]+)\}",joined) for key in group.split(",")}
    missing=sorted(cites-keys)
    if missing:raise SystemExit("Unresolved citation keys: "+", ".join(missing))
    env=[]
    for match in re.finditer(r"\\(begin|end)\{([^}]+)\}",joined):
        op,name=match.groups()
        if op=="begin":env.append(name)
        elif not env or env.pop()!=name:raise SystemExit(f"Mismatched \\end{{{name}}}")
    if env:raise SystemExit("Unclosed TeX environments: "+", ".join(env))
    print(f"Overleaf static check passed: {len(sources)} TeX sources, {len(cites)} citations, all inputs and figures present.")
    print("This check does not compile TeX or assess scientific correctness.")


if __name__=="__main__":main()
