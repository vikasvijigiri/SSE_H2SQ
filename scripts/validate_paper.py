"""Validate the organized manuscript folder without creating an archive."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "docs" / "manuscript"


def main():
    sources = [PAPER / "main.tex", PAPER / "preamble.tex", PAPER / "frontmatter.tex"]
    sources.extend(sorted((PAPER / "sections").glob("*.tex")))
    active = "\n".join(re.sub(r"(?<!\\)%[^\n]*", "", p.read_text(encoding="utf-8")) for p in sources)
    missing = []
    for pattern, suffix in [(r"\\input\{([^}]+)\}", ".tex"),
                            (r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", ""),
                            (r"\\bibliography\{([^}]+)\}", ".bib")]:
        for name in re.findall(pattern, active):
            target = PAPER / (name if Path(name).suffix else name + suffix)
            if not target.is_file():
                missing.append(str(target.relative_to(PAPER)))
    if missing:
        raise SystemExit("Missing manuscript assets: " + ", ".join(missing))
    bibliography = PAPER / "bibliography" / "references.bib"
    keys = set(re.findall(r"@\w+\s*\{\s*([^,\s]+)", bibliography.read_text(encoding="utf-8")))
    citations = {key.strip() for group in re.findall(r"\\cite\w*\*?(?:\[[^]]*\])*\{([^}]+)\}", active) for key in group.split(",")}
    unresolved = sorted(citations - keys)
    if unresolved:
        print("WARNING: existing unresolved bibliography keys: " + ", ".join(unresolved))
    print(f"Checked {len(sources)} TeX sources and {len(citations)} citation references.")


if __name__ == "__main__":
    main()
