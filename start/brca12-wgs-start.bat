@echo off
rem brca12-wgs Start for Windows: runs the start file in WSL 2 (Windows Subsystem for Linux).
rem The pipeline tools (Nextflow, samtools, GATK, ...) are Linux programs and need WSL 2.
setlocal
chcp 65001 >nul
cd /d "%~dp0"
where wsl >nul 2>nul
if errorlevel 1 goto nowsl
wsl -e bash -c "exit 0" >nul 2>nul
if errorlevel 1 goto nowsl
wsl -e bash ./brca12-wgs-start.sh
exit /b %errorlevel%
:nowsl
echo.
echo   WSL 2 with a Linux distribution is required. Install it once in an administrator
echo   PowerShell with:  wsl --install   then restart Windows and start this file again.
echo   Details: README.txt
echo.
echo   Deutsch: WSL 2 mit einer Linux-Distribution wird benoetigt. Einmalig in einer
echo   PowerShell als Administrator:  wsl --install   dann Windows neu starten.
echo.
pause
exit /b 1
