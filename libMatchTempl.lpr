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
  SysUtils, Math, core, matchTempl, matrix, FFTW3, FFTPACK4, threading;

function MatchTemplate_Wrap(constref Image, Templ: T2DIntArray; TMFormula: ETMFormula): T2DRealArray; callconv
begin
  Result := MatchTemplate(Image, Templ, TMFormula);
end;

function LoadFFTWFrom(constref Path: String): LongBool; callconv
begin
  FFTW.Free();
  Result := FFTW.Init([Path], Min(4, GetSystemThreadCount()));
end;

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


// export a couple of stats related method to work with template matching result
function expMin(constref a: T2DRealArray): TReal; callconv
         var _:TReal; begin MinMax(a, Result, _); end;

function expMax(constref a: T2DRealArray): TReal; callconv
         var _:TReal; begin MinMax(a, _, Result); end;

function expArgMin(constref a: T2DRealArray): TPoint; callconv
         begin Result := ArgMin(a); end;

function expArgMax(constref a: T2DRealArray): TPoint; callconv
         begin Result := ArgMax(a); end;

function expNormMinMax(constref a: T2DRealArray; Alpha, Beta: TReal): T2DRealArray; callconv
         begin Result := NormMinMax(a, Alpha, Beta); end;

function expCompareImageAt(constref Image: T2DIntArray; Templ: T2DIntArray; Pt: TPoint; Tol: Int32): Single; callconv
         begin Result := CompareImageAt(Image, Templ, Pt, Tol); end;

function expDownscaleImage(constref Image: T2DIntArray; Scale: Int32): T2DIntArray; callconv
         begin Result := DownscaleImage(Image, Scale); end;

function expRotateImage(constref Image: T2DIntArray; Angle: Single; Expand, Smooth: LongBool): T2DIntArray; callconv
         begin Result := RotateImage(Image, Angle, Expand, Smooth); end;

function expCrop(constref a: T2DIntArray; B: TBox): T2DIntArray; callconv
         begin Result := Crop(a, B); end;

{$I SimbaPlugin.inc}


end.
