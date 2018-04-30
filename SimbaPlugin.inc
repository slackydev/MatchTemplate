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
[==============================================================================}
var
  Methods: array of record procAddr: Pointer; ProcDef:PChar; end;
  TypeDefs: array of record TypeName, TypeDef:PChar; end;
  TypesLoaded: Boolean = False;
  MethodsLoaded: Boolean = False;
  
  OldMemoryManager: TMemoryManager;
  MemIsset: Boolean = False;


procedure AddMethod(ProcAddr: Pointer; ProcDef: PChar);
var L: Integer;
begin
  L := Length(Methods);
  SetLength(Methods, L + 1);
  Methods[l].ProcAddr := ProcAddr;
  Methods[l].ProcDef := ProcDef;
end;

procedure AddType(TypeName, TypeDef: PChar);
var L: Integer;
begin
  L := Length(TypeDefs);
  SetLength(TypeDefs, L + 1);
  TypeDefs[l].TypeName := TypeName;
  TypeDefs[l].TypeDef := TypeDef;
end;


procedure LoadExports;
begin
  AddType('TReal', {$IF SizeOf(TReal)=4}'Single'{$ELSE}'Double'{$ENDIF});
  AddType('TRealArray',   'array of TReal');
  AddType('T2DRealArray', 'array of TRealArray');
  AddType('ETMFormula',   '(TM_CCORR, TM_CCORR_NORMED, TM_CCOEFF, TM_CCOEFF_NORMED, TM_SQDIFF, TM_SQDIFF_NORMED);');


  AddMethod(@MatchTemplate_Wrap, 'function MatchTemplate(constref Img, Sub: T2DIntArray; Formula:ETMFormula=TM_CCOEFF_NORMED): T2DRealArray;');
  AddMethod(@LoadFFTWFrom,       'function  LoadFFTWFrom(constref Path: String): LongBool;');
  AddMethod(@DisableFFTW,        'procedure DisableFFTW();');
  AddMethod(@EnableFFTW,         'function  EnableFFTW(): LongBool');
  AddMethod(@SetMaxFFTThreads,   'procedure SetMaxFFTThreads(MaxThreads: Int32);');
  AddMethod(@EnableFFTCache,     'procedure EnableFFTCache(Enabled: LongBool);');
  AddMethod(@ClearFFTCache,      'procedure ClearFFTCache();');

  // stats helper methods
  AddMethod(@expMin,        'function T2DRealArray.Min(): TReal; constref;');
  AddMethod(@expMax,        'function T2DRealArray.Max(): TReal; constref;');
  AddMethod(@expArgMin,     'function T2DRealArray.ArgMin(): TPoint; constref;');
  AddMethod(@expArgMax,     'function T2DRealArray.ArgMax(): TPoint; constref;');
  AddMethod(@expNormMinMax, 'function T2DRealArray.NormMinMax(A,B: TReal): T2DRealArray; constref;');

  // image helpers
  AddMethod(@expCompareImageAt, 'function T2DIntArray.CompareImageAt(Templ: T2DIntArray; Pt: TPoint; Tol: Int32): Single; constref;');
  AddMethod(@expDownscaleImage, 'function T2DIntArray.DownscaleImage(Scale: Int32): T2DIntArray; constref;');
  AddMethod(@expRotateImage,    'function T2DIntArray.RotateImage(Angle: Single; Expand, Smooth: LongBool): T2DIntArray; constref;');
  AddMethod(@expCrop,           'function T2DIntArray.Crop(B: TBox): T2DIntArray; constref;');

  MethodsLoaded := True;
  TypesLoaded := True;
end;

procedure FreeMethods;
begin
  SetLength(Methods, 0);
  MethodsLoaded := False;
end;

procedure FreeTypes;
begin
  SetLength(TypeDefs, 0);
  TypesLoaded := False;
end;



function GetPluginABIVersion: Integer; callconv export;
begin
  Result := 2;
end;

procedure SetPluginMemManager(MemMgr : TMemoryManager); callconv export;
begin
  if memisset then
    exit;
  GetMemoryManager(OldMemoryManager);
  SetMemoryManager(MemMgr);
  memisset := True;
end;


function GetFunctionCount: Integer; callconv export;
begin
  if not MethodsLoaded then LoadExports;
  Result := Length(Methods);
end;

function GetFunctionInfo(x: Integer; var ProcAddr: Pointer; var ProcDef: PChar): Integer; callconv export;
begin
  Result := x;
  if (x > -1) and InRange(x, 0, High(Methods)) then
  begin
    ProcAddr := Methods[x].procAddr;
    StrPCopy(ProcDef, Methods[x].ProcDef);
    if (x = High(Methods)) then FreeMethods;
  end;
end;

function GetTypeCount: Integer; callconv export;
begin
  if not TypesLoaded then LoadExports;
  Result := Length(TypeDefs);
end;

function GetTypeInfo(x: Integer; var TypeName, TypeDef: PChar): Integer; callconv export;
begin
  Result := x;
  if (x > -1) and InRange(x, 0, High(TypeDefs)) then
  begin
    StrPCopy(TypeName, TypeDefs[x].TypeName);
    StrPCopy(TypeDef,  TypeDefs[x].TypeDef);
    if (x = High(TypeDefs)) then FreeTypes;
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
