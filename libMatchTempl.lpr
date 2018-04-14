library libMatchTempl;
{==============================================================================]
  Copyright (c) 2018, Jarl `slacky` Holta
  Project: libfft
  License: GNU Lesser GPL (http://www.gnu.org/licenses/lgpl.html)
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

function expCompareImageAt(constref Image, Templ: T2DIntArray; Pt: TPoint; Tol: Int32): Single; callconv
         begin Result := CompareImageAt(Image, Templ, Pt, Tol); end;

function expDownscaleImage(constref Image: T2DIntArray; Scale: Int32): T2DIntArray; callconv
         begin Result := DownscaleImage(Image, Scale); end;

function expRotateImage(constref Image: T2DIntArray; Angle: Single; Expand, Smooth: LongBool): T2DIntArray; callconv
         begin Result := RotateImage(Image, Angle, Expand, Smooth); end;

function expCrop(constref a: T2DIntArray; B: TBox): T2DIntArray; callconv
         begin Result := Crop(a, B); end;

{$I SimbaPlugin.inc}


end.
