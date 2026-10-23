@echo off
cd /d "%~dp0"
if not exist .venv\Scripts\python.exe (
  py -3 -m venv .venv
  if errorlevel 1 goto erreur
)
.venv\Scripts\python.exe -m pip install -r requirements.txt
if errorlevel 1 goto erreur
.venv\Scripts\python.exe -m streamlit run app.py
goto :eof
:erreur
echo Installation impossible. Verifiez Python 3.11 ou 3.12 et la connexion Internet.
pause
