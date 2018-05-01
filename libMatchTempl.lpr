library libMatchTempl;
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
{$I header.inc}

uses
  SysUtils, Math, core, matchTempl, matrix, FFTW3, FFTPACK4, cpuinfo;

// -----------------------------------------------------------------------------
// Complex data, FFI is not reliable here
procedure LoadFFTWFrom(const Params: PParamArray; const Result:Pointer); callconv export;
begin
  FFTW.Free();
  LongBool(Result^) := FFTW.Init([AnsiString(Params^[0]^)], Min(4, GetSystemThreadCount()));
end;

// ----------------------------------------------------
// Simple data, FFI is fine
procedure SetMaxFFTThreads(MaxThreads: Int32); callconv
begin
  FFTW.MaxThreads := MaxThreads;
  FFTPACK.MaxThreads := MaxThreads;
end;

procedure DisableFFTW(); callconv
begin
  FFTW.IsLoaded := False;
end;

function EnableFFTW(): LongBool; callconv
begin
  Result := False;
  if FFTW.Handle <> 0 then
  begin
    FFTW.IsLoaded := True;
    Result := True;
  end;
end;

procedure EnableFFTCache(Enabled: LongBool); callconv
begin
  TM_USE_CACHE := Enabled;
end;

procedure ClearFFTCache(); callconv
begin
  ClearCache();
end;


// -----------------------------------------------------------------------------
// export a couple of stats related method to work with template matching result
// complex data, FFI is not fine.

procedure expMatchTemplate(const Params: PParamArray; const Result:Pointer); callconv export;
begin
  T2DRealArray(Result^) := MatchTemplate(T2DIntArray(Params^[0]^), T2DIntArray(Params^[1]^), ETMFormula(Params^[2]^));
end;

procedure expMin(const Params: PParamArray; const Result:Pointer); callconv export;
var _:TReal;
begin
  MinMax(T2DRealArray(Params^[0]^), TReal(Result^), _);
end;

procedure expMax(const Params: PParamArray; const Result:Pointer); callconv export;
var _:TReal;
begin
  MinMax(T2DRealArray(Params^[0]^), _, TReal(Result^));
end;

procedure expArgMin(const Params: PParamArray; const Result:Pointer); callconv export;
begin
  TPoint(Result^) := ArgMin(T2DRealArray(Params^[0]^));
end;

procedure expArgMax(const Params: PParamArray; const Result:Pointer); callconv export;
begin
  TPoint(Result^) := ArgMax(T2DRealArray(Params^[0]^));
end;

procedure expNormMinMax(const Params: PParamArray; const Result:Pointer); callconv export;
begin
  T2DRealArray(Result^) := NormMinMax(T2DRealArray(Params^[0]^), TReal(Params^[1]^), TReal(Params^[2]^));
end;

procedure expCompareImageAt(const Params: PParamArray; const Result:Pointer); callconv export;
begin
  Single(Result^) := CompareImageAt(T2DIntArray(Params^[0]^), T2DIntArray(Params^[1]^), TPoint(Params^[2]^), Int32(Params^[2]^));
end;

procedure expDownscaleImage(const Params: PParamArray; const Result:Pointer); callconv export;
begin
  T2DIntArray(Result^) := DownscaleImage(T2DIntArray(Params^[0]^), Int32(Params^[1]^));
end;

procedure expRotateImage(const Params: PParamArray; const Result:Pointer); callconv export;
begin
  T2DIntArray(Result^) := RotateImage(T2DIntArray(Params^[0]^), Single(Params^[1]^), LongBool(Params^[2]^), LongBool(Params^[3]^));
end;

procedure expCrop(const Params: PParamArray; const Result:Pointer); callconv export;
begin
  T2DIntArray(Result^) := Crop(T2DIntArray(Params^[0]^), TBox(Params^[1]^));
end;

{$I SimbaPlugin.inc}


end.
