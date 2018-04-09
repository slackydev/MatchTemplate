{==============================================================================]
  Copyright (c) 2018, Jarl `slacky` Holta
  Project: libfft
  License: GNU Lesser GPL (http://www.gnu.org/licenses/lgpl.html)
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
  AddType('TReal',   {$IF SizeOf(TReal)=4}'Single'{$ELSE}'Double'{$ENDIF});
  AddType('TRealArray',   'array of TReal');
  AddType('T2DRealArray', 'array of TRealArray');

  AddType('Complex',         'packed record re, im: TReal; end;');
  AddType('TComplexArray',   'array of Complex');
  AddType('T2DComplexArray', 'array of TComplexArray');

  AddMethod(@MatchTemplate, 'function MatchTemplate(const img, sub: T2DIntArray): T2DRealArray;');

  AddMethod(@fftpack.FFTInit, 'function FFTInit(n: Int32): TComplexArray;');
  AddMethod(@fftpack.FFT,     'function FFT(a,wsave: TComplexArray; Inplace: LongBool = False): TComplexArray;');
  AddMethod(@fftpack.IFFT,    'function IFFT(a,wsave: TComplexArray; Inplace: LongBool = False): TComplexArray;');

  AddMethod(@fftpack.RFFTInit, 'function RFFTInit(n: Int32): TRealArray;');
  AddMethod(@fftpack.RFFT,     'function RFFT(a,wsave: TRealArray; Inplace: LongBool = False): TRealArray;');
  AddMethod(@fftpack.IRFFT,    'function IRFFT(a,wsave: TRealArray; Inplace: LongBool = False): TRealArray;');

  AddMethod(@fftpack.FFT2,  'function FFT2(m: T2DComplexArray; Inplace: LongBool = False): T2DComplexArray;');
  AddMethod(@fftpack.IFFT2, 'function IFFT2(m: T2DComplexArray; Inplace: LongBool = False): T2DComplexArray;');

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
