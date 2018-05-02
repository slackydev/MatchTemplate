unit SimbaPlugin;
{==============================================================================]
  Copyright © 2018, Jarl Krister Holta
  
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  -----

  Standalone unit for exporting methods to Simba using version 2 of Simba ABI
[==============================================================================}
interface
{$mode objfpc}{$H+}
{$macro on}

{$DEFINE callconv :=
  {$IFDEF WINDOWS}{$IFDEF CPU32}cdecl;{$ELSE}{$ENDIF}{$ENDIF}
  {$IFDEF LINUX}{$IFDEF CPU32}cdecl;{$ELSE}{$ENDIF}{$ENDIF}
}

uses SysUtils;

var
  Methods: array of record ProcAddr: Pointer; ProcDef: PChar; end;
  TypeDefs: array of record TypeName, TypeDef: PChar; end;
  OldMemoryManager: TMemoryManager;
  MemIsset: Boolean = False;

procedure AddGlobalMethod(ProcAddr: Pointer; ProcDef: PChar);
procedure AddGlobalType(TypeName, TypeDef: PChar);

implementation


procedure AddGlobalMethod(ProcAddr: Pointer; ProcDef: PChar);
var L: Integer;
begin
  L := Length(Methods);
  SetLength(Methods, L + 1);
  Methods[l].ProcAddr := ProcAddr;
  Methods[l].ProcDef := ProcDef;
end;

procedure AddGlobalType(TypeName, TypeDef: PChar);
var L: Integer;
begin
  L := Length(TypeDefs);
  SetLength(TypeDefs, L + 1);
  TypeDefs[l].TypeName := TypeName;
  TypeDefs[l].TypeDef := TypeDef;
end;


// ----------------------------------------------------------------------------
// Methods exported to simba for getting plugin's methods & types

function GetPluginABIVersion: Integer; callconv export;
begin
  Result := 2;
end;

procedure SetPluginMemManager(MemMgr: TMemoryManager); callconv export;
begin
  if MemIsset then Exit;
  GetMemoryManager(OldMemoryManager);
  SetMemoryManager(MemMgr);
  memisset := True;
end;

function GetFunctionCount: Integer; callconv export;
begin
  Result := Length(Methods);
end;

function GetFunctionInfo(x: Integer; var ProcAddr: Pointer; var ProcDef: PChar): Integer; callconv export;
begin
  Result := x;
  if (x > -1) and (x < Length(Methods)) then
  begin
    ProcAddr := Methods[x].ProcAddr;
    StrPCopy(ProcDef, Methods[x].ProcDef);
  end;
end;

function GetTypeCount: Integer; callconv export;
begin
  Result := Length(TypeDefs);
end;

function GetTypeInfo(x: Integer; var TypeName, TypeDef: PChar): Integer; callconv export;
begin
  Result := x;
  if (x > -1) and (x < Length(TypeDefs)) then
  begin
    StrPCopy(TypeName, TypeDefs[x].TypeName);
    StrPCopy(TypeDef,  TypeDefs[x].TypeDef);
  end;
end;

procedure OnDetach; callconv export;
begin
  SetMemoryManager(OldMemoryManager);
end;

exports GetPluginABIVersion;
exports SetPluginMemManager;
exports GetFunctionCount;
exports GetFunctionInfo;
exports GetTypeCount;
exports GetTypeInfo;
exports OnDetach; 

end.
