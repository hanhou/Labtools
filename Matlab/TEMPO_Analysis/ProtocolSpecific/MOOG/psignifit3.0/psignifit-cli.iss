; This setup script requires the script modpath.iss written by Jared Breland
; currently, the script is not in the psignifit source distribution due to unresolved
; license issues. If you want to build your own windows installer, you have to
; download the script yourself from
; http://www.legroom.net/software


[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{22600284-5AE4-44F6-BF95-7CD80B34E86E}
AppName=psignifit
AppVersion=3.0 beta
;AppVerName=psignifit-cli 3.0 beta
AppPublisher=psignifit development team
AppPublisherURL=http://psignifit.sourceforge.net
AppSupportURL=http://psignifit.sourceforge.net
AppUpdatesURL=http://psignifit.sourceforge.net
DefaultDirName={pf}\psignifit3
DefaultGroupName=psignifit
LicenseFile=COPYING
InfoBeforeFile=doc-src\WELCOME.rst
OutputBaseFilename=psignifit-cli_3_beta_installer
OutputDir=WindowsInstaller
Compression=lzma
SolidCompression=yes
ChangesEnvironment=true

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Files]
Source: "cli\*.exe"; DestDir: "{app}"; Flags: ignoreversion; Components: installcli;
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{group}\{cm:UninstallProgram,psignifit}"; Filename: "{uninstallexe}"

[Types]
Name: full; Description: "Full installation";
Name: cli;  Description: "Only install the command line interface";

[Components]
Name: installcli; Description: "Command line interface (required for the matlab interface)"; Types: full cli;

[Tasks]
Name: modifypath; Description: Add psignifit-cli to your environment path (Highly recommended); Flags: checkedonce

[Code]
const
  ModPathName = 'modifypath';
  ModPathType = 'user';

function ModPathDir(): TArrayOfString;
begin
  setArrayLength(Result, 1)
  Result[0] := ExpandConstant('{app}');
end;
#include "modpath.iss"
