[Setup]
AppName=FASTER STR Analysis
AppVersion=1.0.0
DefaultDirName={pf}\FASTER
DefaultGroupName=FASTER
OutputDir=dist
OutputBaseFilename=FASTER-Setup
Compression=lzma
SolidCompression=yes
WizardStyle=modern
DisableProgramGroupPage=yes
DisableDirPage=no
DisableReadyPage=no
DisableFinishedPage=no
SetupIconFile=src\faster\assets\icon.ico

[Files]
Source: "dist\FASTER-GUI.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "src\faster\config\marker_info.json"; DestDir: "{app}\config"; Flags: ignoreversion
Source: "src\faster\config\variant_catalog.thermofisher_24markers.json"; DestDir: "{app}\config"; Flags: ignoreversion
Source: "src\faster\bin\ExpansionHunter"; DestDir: "{app}\bin"; Flags: ignoreversion
Source: "README.md"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\FASTER STR Analysis"; Filename: "{app}\FASTER-GUI.exe"; IconFilename: "{app}\icon.ico"
Name: "{group}\Uninstall FASTER"; Filename: "{uninstallexe}"
Name: "{commondesktop}\FASTER STR Analysis"; Filename: "{app}\FASTER-GUI.exe"; IconFilename: "{app}\icon.ico"

[Run]
Filename: "{app}\FASTER-GUI.exe"; Description: "Launch FASTER STR Analysis"; Flags: postinstall nowait

[Code]
function InitializeSetup(): Boolean;
begin
  Result := True;
end;