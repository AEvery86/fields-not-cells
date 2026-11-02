@echo off
setlocal
cd /d "%~dp0"
set "RS="
for /f "delims=" %%i in ('where Rscript 2^>nul') do if not defined RS set "RS=%%i"
if not defined RS for /d %%D in ("%ProgramFiles%\R\R-*") do set "RS=%%D\bin\Rscript.exe"
if not defined RS (
  echo Rscript was not found. Please install R first: R-4.6.2-win.exe in your Downloads folder.
  pause
  exit /b 1
)
echo Using "%RS%"
"%RS%" 00_install_packages.R
if errorlevel 1 (
  echo Package installation failed. Please send the messages above to Claude.
  pause
  exit /b 1
)
"%RS%" 01_validate_A7.R full
echo.
echo Finished. The report is output\R_vs_Python_report_full.txt
pause
