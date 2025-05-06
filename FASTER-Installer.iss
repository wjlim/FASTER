[Setup]
AppName=FASTER STR Analysis
AppVersion=1.0.0
DefaultDirName={pf}\FASTER
DefaultGroupName=FASTER
OutputDir=dist
OutputBaseFilename=FASTER-Setup
Compression=lzma
SolidCompression=yes

[Files]
Source: "dist\FASTER-GUI.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "src\faster\config\marker_info.json"; DestDir: "{app}\config"; Flags: ignoreversion

[Icons]
Name: "{group}\FASTER STR Analysis"; Filename: "{app}\FASTER-GUI.exe"
Name: "{group}\Uninstall FASTER"; Filename: "{uninstallexe}"