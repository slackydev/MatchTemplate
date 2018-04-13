unit matchTempl;
{==============================================================================]
  Copyright (c) 2018, Jarl `slacky` Holta
  Project: libfft
  License: GNU Lesser GPL (http://www.gnu.org/licenses/lgpl.html)
[==============================================================================}
{$I header.inc}

interface

uses
  SysUtils, core;

type
  ETMFormula = (TM_CCORR, TM_CCORR_NORMED, TM_CCOEFF, TM_CCOEFF_NORMED, TM_SQDIFF, TM_SQDIFF_NORMED);

function MatchTemplate(constref img, sub: T2DIntArray; TMFormula: ETMFormula): T2DRealArray;


implementation

uses
  math, matrix, threading, FFTPACK4, FFTW3;


procedure InitMatrix(var a: T2DRealArray; H,W: Int32; InitValue:Int32);
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

function CCORR(big, sub: T2DRealArray): T2DRealArray;
var
  x,y,w,h,tw,th: Int32;
  a,b: T2DComplexArray;
begin
  Size(big, w,h);
  Size(sub, tw,th);

  if FFTW.IsLoaded then
  begin
    SetLength(big, FFTW.OptimalDFTSize(H), FFTW.OptimalDFTSize(W));
    SetLength(sub, FFTW.OptimalDFTSize(H), FFTW.OptimalDFTSize(W));
    a := FFTW.FFT2_R2C(big);
    b := FFTW.FFT2_R2C(sub);
  end else
  begin
    SetLength(a, FFTPACK.OptimalDFTSize(H), FFTPACK.OptimalDFTSize(W));
    SetLength(b, FFTPACK.OptimalDFTSize(H), FFTPACK.OptimalDFTSize(W));
    for y:=0 to w-1  do for x:=0 to w-1  do a[y,x].re := big[y,x];
    for y:=0 to th-1 do for x:=0 to tw-1 do b[y,x].re := sub[y,x];
    a := FFTPACK.FFT2(a);
    b := FFTPACK.FFT2(b);
  end;

  b := MulSpectrumConj(a,b);

  if FFTW.IsLoaded then
  begin
    Result := FFTW.FFT2_C2R(b);
    SetLength(Result, H,W);
    Exit;
  end else
  begin
    b := FFTPACK.IFFT2(b);
    SetLength(Result, H,W);
    for y:=0 to H-1 do for x:=0 to W-1 do Result[y,x] := b[y,x].re;
  end;

  SetLength(Result, h-th+1,w-tw+1);
end;

function CCORR_NORMED(a,t: T2DRealArray): T2DRealArray;
var
  x,y,tw,th,aw,ah: Int32;
  invSize, numer, denom, tplSdv, tplMean, tplSigma, wndSum2: Double;
  xcorr: T2DRealArray;
  sum, sum2: T2DDoubleArray;
begin
  xcorr := CCORR(a,t);
  Size(t, tw,th);
  invSize := Double(1.0) / Double(tw*th);
  MeanStdev(t, tplMean, tplSdv);
  tplSigma := Sqrt(Sqr(tplSdv) + Sqr(tplMean)) / Sqrt(invSize);

  sum := SumsPd(a, sum2);
  Size(sum, aw,ah);
  SetLength(Result, ah-th, aw-tw);
  for y:=0 to ah-th-1 do
    for x:=0 to aw-tw-1 do
    begin
      wndSum2 := (sum2[Y,X] - sum2[Y,X+tw] - sum2[Y+th,X] + sum2[Y+th,X+tw]);
      numer := xcorr[y,x];
      denom := tplSigma * Sqrt(wndSum2);

      if abs(numer) < denom then
        Result[y,x] := numer / denom
      else if abs(numer) < denom*1.25 then
        if numer > 0 then Result[y,x] := 1 else Result[y,x] := -1;
    end;
end;

function CCOEFF(a,t: T2DRealArray; Normed: Boolean): T2DRealArray;
var
  x,y,tw,th,aw,ah: Int32;
  invSize, numer, denom, tplSdv, tplMean, tplSigma, wndDiff, wndSum: Double;
  xcorr: T2DRealArray;
  sum, sum2: T2DDoubleArray;
begin
  xcorr := CCORR(a,t);
  Size(t, tw,th);

  invSize := Double(1.0) / Double(tw*th);
  MeanStdev(t, tplMean, tplSdv);
  tplSigma := tplSdv / Sqrt(invSize);

  if tplSdv < 0.00001 then
  begin
    InitMatrix(Result, Length(xcorr), Length(xcorr[0]), 1);
    Exit;
  end;

  sum := SumsPd(a, sum2);
  Size(sum, aw,ah);
  SetLength(Result, ah-th, aw-tw);
  for y:=0 to ah-th-1 do
    for x:=0 to aw-tw-1 do
    begin
      wndSum := sum[Y,X] - sum[Y,X+tw] - sum[Y+th,X] + sum[Y+th,X+tw];
      numer  := xcorr[y,x] - (wndSum * tplMean);

      if Normed then
      begin
        wndDiff := (sum2[Y,X] - sum2[Y,X+tw] - sum2[Y+th,X] + sum2[Y+th,X+tw]) - (Sqr(wndSum) * invSize);
        denom   := tplSigma * Sqrt(Max(0,wndDiff));

        if abs(numer) < denom then
          Result[y,x] := numer / denom
        else if abs(numer) < denom*1.25 then
          if numer > 0 then Result[y,x] := 1 else Result[y,x] := -1;
      end else
        Result[y,x] := numer;
    end;
end;

function SQDIFF(a,t: T2DRealArray; Normed: Boolean): T2DRealArray;
var
  x,y,tw,th,aw,ah: Int32;
  invSize, numer, denom, tplSdv, tplMean, tplSigma, tplSum2, wndSum2: Double;
  xcorr: T2DRealArray;
  sum2: T2DDoubleArray;
begin
  xcorr := CCORR(a,t);
  Size(t, tw,th);

  invSize := Double(1.0) / Double(tw*th);
  MeanStdev(t, tplMean, tplSdv);
  tplSigma := Sqrt(Sqr(tplSdv) + Sqr(tplMean)) / Sqrt(invSize);
  tplSum2  := (Sqr(tplSdv) + Sqr(tplMean)) / invSize;

  SumsPd(a, sum2);
  Size(sum2, aw,ah);
  SetLength(Result, ah-th, aw-tw);
  for y:=0 to ah-th-1 do
    for x:=0 to aw-tw-1 do
    begin
      wndSum2 := sum2[Y,X] - sum2[Y,X+tw] - sum2[Y+th,X] + sum2[Y+th,X+tw];
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

function CCORR_RGB(constref img, sub: T2DIntArray; Normed: Boolean): T2DRealArray;
var
  R1,G1,B1, R2,G2,B2: T2DRealArray;
  x,y,W,H: Int32;
begin
  Size(img, W,H);
  SplitRGB(img, R1,G1,B1);
  SplitRGB(sub, R2,G2,B2);

  if FFTW.IsLoaded then
    FFTW.PrepareThreads(W*H);

  if Normed then
  begin
    R1 := CCORR_NORMED(R1,R2);  R2 := nil;
    G1 := CCORR_NORMED(G1,G2);  G2 := nil;
    B1 := CCORR_NORMED(B1,B2);  B2 := nil;
  end else
  begin
    R1 := CCORR(R1,R2);  R2 := nil;
    G1 := CCORR(G1,G2);  G2 := nil;
    B1 := CCORR(B1,B2);  B2 := nil;
  end;

  Result := R1;
  Size(R1, w,h);
  for y:=0 to h-1 do
    for x:=0 to w-1 do
      if Normed then
        Result[y,x] := (R1[y,x] + B1[y,x] + G1[y,x]) * 0.333333333334
      else
        Result[y,x] := (R1[y,x] + B1[y,x] + G1[y,x]);
end;

function CCOEFF_RGB(constref img, sub: T2DIntArray; Normed:Boolean): T2DRealArray;
var
  R1,G1,B1, R2,G2,B2: T2DRealArray;
  x,y,W,H: Int32;
begin
  Size(img, W,H);
  SplitRGB(img, R1,G1,B1);
  SplitRGB(sub, R2,G2,B2);

  if FFTW.IsLoaded then
    FFTW.PrepareThreads(W*H);

  R1 := CCOEFF(R1,R2, Normed);  R2 := nil;
  G1 := CCOEFF(G1,G2, Normed);  G2 := nil;
  B1 := CCOEFF(B1,B2, Normed);  B2 := nil;

  Result := R1;
  Size(Result, W,H);
  for y:=0 to H-1 do
    for x:=0 to W-1 do
      if Normed then
        Result[y,x] := (R1[y,x] + G1[y,x] + B1[y,x]) * 0.333333333334
      else
        Result[y,x] := (R1[y,x] + G1[y,x] + B1[y,x]);
end;

function SQDIFF_RGB(constref img, sub: T2DIntArray; Normed:Boolean): T2DRealArray;
var
  R1,G1,B1, R2,G2,B2: T2DRealArray;
  x,y,W,H: Int32;
begin
  Size(img, W,H);
  SplitRGB(img, R1,G1,B1);
  SplitRGB(sub, R2,G2,B2);

  if FFTW.IsLoaded then
    FFTW.PrepareThreads(W*H);

  R1 := SQDIFF(R1,R2, Normed);  R2 := nil;
  G1 := SQDIFF(G1,G2, Normed);  G2 := nil;
  B1 := SQDIFF(B1,B2, Normed);  B2 := nil;

  Result := R1;
  Size(Result, W,H);
  for y:=0 to H-1 do
    for x:=0 to W-1 do
      if Normed then
        Result[y,x] := (R1[y,x] + G1[y,x] + B1[y,x]) * 0.333333333334
      else
        Result[y,x] := (R1[y,x] + G1[y,x] + B1[y,x]);
end;


// ----------------------------------------------------------------------------


function MatchTemplate(constref img, sub: T2DIntArray; TMFormula: ETMFormula): T2DRealArray;
begin
  case TMFormula of
    TM_CCORR:         Result := CCORR_RGB(img,sub, False);
    TM_CCORR_NORMED:  Result := CCORR_RGB(img,sub, True);
    TM_CCOEFF:        Result := CCOEFF_RGB(img,sub, False);
    TM_CCOEFF_NORMED: Result := CCOEFF_RGB(img,sub, True);
    TM_SQDIFF:        Result := SQDIFF_RGB(img,sub, False);
    TM_SQDIFF_NORMED: Result := SQDIFF_RGB(img,sub, True);
    else
      raise Exception.Create('Not implemented');
  end;
end;


end.
