# Paper folder

This directory is the manuscript project. Open `main.tex` to edit or compile it.
For Overleaf, add these source files and folders to a project, select `main.tex`
as the main document, and use pdfLaTeX with BibTeX.

```text
docs/manuscript/
|-- main.tex                     # Entry point and section order
|-- preamble.tex                 # REVTeX class, packages, macros
|-- frontmatter.tex              # Title, authors, abstract
|-- sections/
|   |-- 01-introduction.tex
|   |-- 02-model-and-method.tex
|   |-- 03-physical-quantities.tex
|   |-- 04-results.tex
|   |-- 05-discussion.tex
|   |-- 06-acknowledgments.tex
|   `-- 07-appendix.tex
|-- figures/                     # Paper figure assets
|-- bibliography/
|   `-- references.bib            # Citation database
|-- references/                  # Background reading PDFs
|-- latexmkrc                    # Build configuration
|-- compile.ps1                  # Windows PowerShell build
|-- compile.bat                  # Windows batch build
`-- README.md
```

Build from this folder with `latexmk -pdf main.tex`; clean auxiliary files with
`latexmk -c`. From the repository root, `python scripts/validate_paper.py` checks
source and figure paths and reports unresolved bibliography keys.

The existing REVTeX 4.1 / PRD formatting and scientific text are preserved.
Before submission, resolve citation `Ref`, review visible editorial text, and
reconcile the Suzuki-Trotter methods description with the forthcoming SLQMC/SSE
benchmark. PDF compilation remains unverified because local MiKTeX setup is incomplete.
