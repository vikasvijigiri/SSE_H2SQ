@echo off
setlocal
pushd "%~dp0"
latexmk -pdf main.tex
set "BUILD_EXIT=%ERRORLEVEL%"
popd
exit /b %BUILD_EXIT%
