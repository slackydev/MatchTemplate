library libMatchTempl;
{==============================================================================]
  Copyright (c) 2018, Jarl `slacky` Holta
  Project: libfft
  License: GNU Lesser GPL (http://www.gnu.org/licenses/lgpl.html)
[==============================================================================}
{$I header.inc}

uses
  SysUtils, Math, core, matchTempl, FFTW3, FFTPACK4, threading;

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
    FFTW.IsLoaded := False;
    Result := True;
  end;
end;


{$I SimbaPlugin.inc}


end.
