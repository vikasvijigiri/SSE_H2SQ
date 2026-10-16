$ErrorActionPreference = "Stop"
Push-Location $PSScriptRoot
try {
    & latexmk -pdf main.tex
    if ($LASTEXITCODE -ne 0) { throw "Manuscript compilation failed ($LASTEXITCODE)" }
} finally { Pop-Location }
