set DIR_OF_THIS_SCRIPT=%~dp0

cd "%DIR_OF_THIS_SCRIPT%..\"

python.exe .\cwlab.py up --config "%DIR_OF_THIS_SCRIPT%\config_windows.yml"


