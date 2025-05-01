@echo off
REM Build FASTER GUI as a standalone Windows executable
%~dp0.venv\Scripts\pyinstaller.exe --onefile --noconsole --name FASTER-GUI --add-data "src/faster/config/marker_info.json;config" src/faster/gui.py

echo.
echo EXE created: dist\FASTER-GUI.exe 