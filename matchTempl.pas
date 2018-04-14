unit matchTempl;
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

interface

uses
  SysUtils, core;

type
  ETMFormula = (TM_CCORR, TM_CCORR_NORMED, TM_CCOEFF, TM_CCOEFF_NORMED, TM_SQDIFF, TM_SQDIFF_NORMED);

function MatchTemplate(constref Image, Templ: T2DIntArray; TMFormula: ETMFormula): T2DRealArray;


implementation

uses
  math, matrix, threading, FFTPACK4, FFTW3;


procedure InitMatrix(out a: T2DRealArray; H,W: Int32; InitValue:Int32);
var
  x,y: Int32;
begin
  SetLength(a, H,W);
  for y:=0 to H-1 do
    for x:=0 to W-1 do a[y,x] := InitValue;
end;


// --------------------------------------------------------------------------------
// a * conj(b)

procedure __mulSpectrumConj_thread(params: PParamArray);
var
  x,y: Int32;
  a,b: T2DComplexArray;
  box: TBox;
  re,im: TReal;
begin
  a := T2DComplexArray(Params^[0]^);
  b := T2DComplexArray(Params^[1]^);
  box := PBox(Params^[2])^;

  for y:=box.y1 to box.y2 do
    for x:=box.x1 to box.x2 do
      begin
        re := (a[y,x].re *  b[y,x].re) - (a[y,x].im * -b[y,x].im);
        im := (a[y,x].re * -b[y,x].im) + (a[y,x].im *  b[y,x].re);
        b[y,x].re := re;
        b[y,x].im := im;
      end;
end;

function MulSpectrumConj(a,b: T2DComplexArray): T2DComplexArray;
var
  W,H,tc: Int32;
  ThreadPool: TThreadPool;
begin
  tc := GetSystemThreadCount div 2;
  ThreadPool := TThreadPool.Create(tc);
  H := Length(a);
  if H = 0 then Exit;
  W := Length(a[0]);
  ThreadPool.MatrixFunc(@__mulSpectrumConj_thread, [@a,@b], W,H, tc, 300*300);
  ThreadPool.Free();
  Result := b;
end;


// --------------------------------------------------------------------------------
// cross correlate

function CCORR(Image, Templ: T2DRealArray): T2DRealArray;
var
  x,y,aw,ah,tw,th: Int32;
  a,b: T2DComplexArray;
begin
  Size(Image, aw,ah);
  Size(Templ, tw,th);

  if FFTW.IsLoaded then
  begin
    SetLength(Image, FFTW.OptimalDFTSize(aH), FFTW.OptimalDFTSize(aW));
    SetLength(Templ, FFTW.OptimalDFTSize(aH), FFTW.OptimalDFTSize(aW));
    a := FFTW.FFT2_R2C(Image);
    b := FFTW.FFT2_R2C(Templ);
  end else
  begin
    SetLength(a, FFTPACK.OptimalDFTSize(aH), FFTPACK.OptimalDFTSize(aW));
    SetLength(b, FFTPACK.OptimalDFTSize(aH), FFTPACK.OptimalDFTSize(aW));
    for y:=0 to ah-1 do for x:=0 to aw-1 do a[y,x].re := Image[y,x];
    for y:=0 to th-1 do for x:=0 to tw-1 do b[y,x].re := Templ[y,x];
    a := FFTPACK.FFT2(a);
    b := FFTPACK.FFT2(b);
  end;

  b := MulSpectrumConj(a,b);

  if FFTW.IsLoaded then
  begin
    Result := FFTW.FFT2_C2R(b);
    SetLength(Result, ah-th+1, aw-tw+1);
  end else
  begin
    b := FFTPACK.IFFT2(b);
    SetLength(Result, ah-th+1, aw-tw+1);
    for y:=0 to ah-th do for x:=0 to aw-tw do Result[y,x] := b[y,x].re;
  end;
end;


// ----------------------------------------------------------------------------

function CCORR_RGB_HELPER(Image, Templ: T2DIntArray; out aR,aG,aB, tR,tG,tB: T2DRealArray): T2DRealArray;
var
  x,y,W,H: Int32;
  xR,xG,xB: T2DRealArray;
begin
  SplitRGB(Image, aR,aG,aB);
  SplitRGB(Templ, tR,tG,tB);

  xR := CCORR(aR,tR);
  xG := CCORR(aG,tG);
  xB := CCORR(aB,tB);

  Result := xR;
  Size(Result, W,H);
  for y:=0 to H-1 do
    for x:=0 to W-1 do
      Result[y,x] := xR[y,x] + xG[y,x] + xB[y,x];
end;



function CCORR_RGB(Image, Templ: T2DIntArray; Normed: Boolean): T2DRealArray;
var
  x,y,tw,th,aw,ah: Int32;
  invSize, numer, denom, tplSdv, tplMean, tplSigma, mR,sR, mG,sG, mB,sB, wndSum2: Double;
  sum2r, sum2g, sum2b: T2DDoubleArray;
  xcorr, aR,aG,aB, tR,tG,tB: T2DRealArray;
begin
  xcorr := CCORR_RGB_HELPER(Image, Templ, aR,aG,aB, tR,tG,tB);

  if not Normed then
    Exit(xcorr);

  Size(Templ, tw,th);
  invSize := Double(1.0) / Double(tw*th);

  MeanStdev(tR, mR, sR); tR := nil;
  MeanStdev(tG, mG, sG); tG := nil;
  MeanStdev(tB, mB, sB); tB := nil;

  tplMean := Sqr(mR) + Sqr(mG) + Sqr(mB);
  tplSdv  := Sqr(sR) + Sqr(sG) + Sqr(sB);

  tplSigma := Sqrt(tplSdv + tplMean) / Sqrt(invSize);

  SumsPd(aR, sum2r); aR := nil;
  SumsPd(aG, sum2g); aG := nil;
  SumsPd(aB, sum2b); aB := nil;

  Size(sum2r, aw,ah);
  SetLength(Result, ah-th, aw-tw);
  for y:=0 to ah-th-1 do
    for x:=0 to aw-tw-1 do
    begin
      wndSum2 := sum2r[Y,X] - sum2r[Y,X+tw] - sum2r[Y+th,X] + sum2r[Y+th,X+tw];
      wndSum2 += sum2g[Y,X] - sum2g[Y,X+tw] - sum2g[Y+th,X] + sum2g[Y+th,X+tw];
      wndSum2 += sum2b[Y,X] - sum2b[Y,X+tw] - sum2b[Y+th,X] + sum2b[Y+th,X+tw];

      numer := xcorr[y,x];
      denom := tplSigma * Sqrt(wndSum2);

      if abs(numer) < denom then
        Result[y,x] := numer / denom
      else if abs(numer) < denom*1.25 then
        if numer > 0 then Result[y,x] := 1 else Result[y,x] := -1;
    end;
end;

function CCOEFF_RGB(Image, Templ: T2DIntArray; Normed: Boolean): T2DRealArray;
var
  x,y,tw,th,aw,ah: Int32;
  invSize, numer, denom, tplSdv, tplSigma, wndSum2, wndMean2: Double;
  wndSumR, wndSumG, wndSumB: Double;
  mR,sR, mG,sG, mB,sB: Double;
  sumR, sumG, sumB, sum2r, sum2g, sum2b: T2DDoubleArray;
  xcorr, aR,aG,aB, tR,tG,tB: T2DRealArray;
begin
  xcorr := CCORR_RGB_HELPER(Image, Templ, aR,aG,aB, tR,tG,tB);

  Size(Templ, tw,th);
  invSize := Double(1.0) / Double(tw*th);

  if not Normed then
  begin
    mR := Mean(tR); tR := nil;
    mG := Mean(tG); tG := nil;
    mB := Mean(tB); tB := nil;
  end else
  begin
    MeanStdev(tR, mR, sR); tR := nil;
    MeanStdev(tG, mG, sG); tG := nil;
    MeanStdev(tB, mB, sB); tB := nil;

    tplSdv  := Sqr(sR) + Sqr(sG) + Sqr(sB);

    if tplSdv < 0.00001 then
    begin
      InitMatrix(Result, Length(xcorr), Length(xcorr[0]), 1);
      Exit;
    end;

    tplSigma := Sqrt(tplSdv) / Sqrt(invSize);
  end;

  sumR := SumsPd(aR, sum2r); aR := nil;
  sumG := SumsPd(aG, sum2g); aG := nil;
  sumB := SumsPd(aB, sum2b); aB := nil;

  Size(sumR, aw,ah);
  SetLength(Result, ah-th, aw-tw);
  for y:=0 to ah-th-1 do
    for x:=0 to aw-tw-1 do
    begin
      wndSumR  := sumR[Y,X] - sumR[Y,X+tw] - sumR[Y+th,X] + sumR[Y+th,X+tw];
      wndSumG  := sumG[Y,X] - sumG[Y,X+tw] - sumG[Y+th,X] + sumG[Y+th,X+tw];
      wndSumB  := sumB[Y,X] - sumB[Y,X+tw] - sumB[Y+th,X] + sumB[Y+th,X+tw];

      wndMean2 := Sqr(wndSumR) + Sqr(wndSumG) + Sqr(wndSumB);
      numer    := xcorr[y,x] - ((wndSumR * mR) + (wndSumG * mG) + (wndSumB * mB));
      if Normed then
      begin
        wndSum2  := sum2r[Y,X] - sum2r[Y,X+tw] - sum2r[Y+th,X] + sum2r[Y+th,X+tw];
        wndSum2  += sum2g[Y,X] - sum2g[Y,X+tw] - sum2g[Y+th,X] + sum2g[Y+th,X+tw];
        wndSum2  += sum2b[Y,X] - sum2b[Y,X+tw] - sum2b[Y+th,X] + sum2b[Y+th,X+tw];
        wndMean2 := wndMean2 * invSize;

        denom := tplSigma * Sqrt(Max(0, wndSum2 - wndMean2));
        if abs(numer) < denom then
          Result[y,x] := numer / denom
        else if abs(numer) < denom*1.25 then
          if numer > 0 then Result[y,x] := 1 else Result[y,x] := -1;
      end else
        Result[y,x] := numer;
    end;
end;

function SQDIFF_RGB(Image, Templ: T2DIntArray; Normed: Boolean): T2DRealArray;
var
  x,y,tw,th,aw,ah: Int32;
  invSize, numer, denom, tplSigma, tplSum2, wndSum2: Double;
  tplMean, tplSdv, mR,sR, mG,sG, mB,sB:Double;
  sum2r, sum2g, sum2b: T2DDoubleArray;
  xcorr, aR,aG,aB, tR,tG,tB: T2DRealArray;
begin
  xcorr := CCORR_RGB_HELPER(Image, Templ, aR,aG,aB, tR,tG,tB);

  Size(Templ, tw,th);
  invSize := Double(1.0) / Double(tw*th);

  MeanStdev(tR, mR, sR); tR := nil;
  MeanStdev(tG, mG, sG); tG := nil;
  MeanStdev(tB, mB, sB); tB := nil;

  tplMean := Sqr(mR) + Sqr(mG) + Sqr(mB);
  tplSdv  := Sqr(sR) + Sqr(sG) + Sqr(sB);

  tplSigma := Sqrt(tplSdv + tplMean) / Sqrt(invSize);
  tplSum2  := (tplSdv + tplMean) / invSize;

  SumsPd(aR, sum2r); aR := nil;
  SumsPd(aG, sum2g); aG := nil;
  SumsPd(aB, sum2b); aB := nil;

  Size(sum2r, aw,ah);
  SetLength(Result, ah-th, aw-tw);
  for y:=0 to ah-th-1 do
    for x:=0 to aw-tw-1 do
    begin
      wndSum2 := sum2r[Y,X] - sum2r[Y,X+tw] - sum2r[Y+th,X] + sum2r[Y+th,X+tw];
      wndSum2 += sum2g[Y,X] - sum2g[Y,X+tw] - sum2g[Y+th,X] + sum2g[Y+th,X+tw];
      wndSum2 += sum2b[Y,X] - sum2b[Y,X+tw] - sum2b[Y+th,X] + sum2b[Y+th,X+tw];

      numer   := Max(0, wndSum2 - Double(2.0)*xcorr[y,x] + tplSum2);
      if Normed then begin
        denom := tplSigma * Sqrt(wndSum2);
        if abs(numer) < denom then
          Result[y,x] := numer / denom
        else
          Result[y,x] := 1;
      end else
        Result[y,x] := numer;
    end;
end;


// ----------------------------------------------------------------------------


function MatchTemplate(constref Image, Templ: T2DIntArray; TMFormula: ETMFormula): T2DRealArray;
begin
  if FFTW.IsLoaded then
    FFTW.PrepareThreads(Length(Image)*Length(Image[0]));

  case TMFormula of
    TM_CCORR:         Result := CCORR_RGB(Image, Templ, False);
    TM_CCORR_NORMED:  Result := CCORR_RGB(Image, Templ, True);
    TM_CCOEFF:        Result := CCOEFF_RGB(Image, Templ, False);
    TM_CCOEFF_NORMED: Result := CCOEFF_RGB(Image, Templ, True);
    TM_SQDIFF:        Result := SQDIFF_RGB(Image, Templ, False);
    TM_SQDIFF_NORMED: Result := SQDIFF_RGB(Image, Templ, True);
    else
      raise Exception.Create('Not implemented');
  end;
end;


end.
