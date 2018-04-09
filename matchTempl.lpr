library matchTempl;
{==============================================================================]
  Copyright (c) 2018, Jarl `slacky` Holta
  Project: libfft
  License: GNU Lesser GPL (http://www.gnu.org/licenses/lgpl.html)
[==============================================================================}
{$I header.inc}

uses
  SysUtils, Math, 
  core, matrix, threading, fftpack;

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

function CrossCorr_1c(const img, sub: T2DRealArray): T2DRealArray; callconv
var
  x,y,H,W,kH,kW: Int32;
  a,b: T2DComplexArray;
begin
  Size(img, W,H);
  Size(sub, kW,kH);

  SetLength(a, OptimalDFTSize(H), OptimalDFTSize(W));
  SetLength(b, OptimalDFTSize(H), OptimalDFTSize(W));

  for y:=0 to H-1  do for x:=0 to W-1  do a[y,x].re := img[y,x];
  for y:=0 to kH-1 do for x:=0 to kW-1 do b[y,x].re := sub[y,x];

  a := FFT2MT(a, False);
  b := FFT2MT(b, False);
  b := MulSpectrumConj(a,b);
  b := FFT2MT(b, True);

  SetLength(Result, H,W);
  for y:=0 to H-1 do for x:=0 to W-1 do Result[y,x] := b[y,x].re;
end;

function MatchTemplate_1c(const a,t: T2DRealArray): T2DRealArray; callconv
var
  x,y,tw,th,aw,ah: Int32;
  sigma_t, denom, numer, invSize, std_t, mean_t: Single;
  ls_a, xcorr: T2DRealArray;

  ls_diff: Double;
  ls2_a: T2DDoubleArray;
begin
  // [ShiftResize] Result dimensions must match LocalSum results
  xcorr := CrossCorr_1c(ShiftResize(a, High(t), High(t[0])), t);

  Size(t, tw,th);
  invSize := 1.0 / (tw*th);
  MeanStdev(t, mean_t, std_t);
  sigma_t := Sqrt(tw*th-1) * std_t;

  ls_a  := LocalSum(a, th,tw, 0);
  ls2_a := LocalSum(SqrToDouble(a), th,tw, 0);

  Size(ls_a, aw,ah);
  SetLength(Result, ah,aw);
  Dec(th); Dec(tw);
  for y:=th to ah-th+1 do
    for x:=tw to aw-tw+1 do
    begin
      ls_diff := ls2_a[y,x] - (Sqr(ls_a[y,x]) * invSize);
      if ls_diff > 0 then
      begin
        denom := sigma_t * Sqrt(ls_diff);
        numer := xcorr[y,x] - ls_a[y,x] * mean_t;
        if denom <> 0 then Result[y-th,x-tw] := numer / denom;
      end;
    end;

  Size(a, aw,ah);
  SetLength(Result, ah-th,aw-tw);
end;

function MatchTemplate(const img, sub: T2DIntArray): T2DRealArray; callconv
var
  R1,G1,B1, R2,G2,B2: T2DRealArray;
  G,B: T2DRealArray;
  x,y,W,H: Int32;
begin
  SplitRGB(img, R1,G1,B1);
  SplitRGB(sub, R2,G2,B2);
  Result := MatchTemplate_1c(R1,R2);       R1 := nil; R2 := nil;
  G      := MatchTemplate_1c(G1,G2);       G1 := nil; G2 := nil;
  B      := MatchTemplate_1c(B1,B2);       B1 := nil; B2 := nil;

  // This is imperfect, but my guess is for this usage, it's OK
  Size(Result, W,H);
  for y:=0 to H-1 do
    for x:=0 to W-1 do
      Result[y,x] := (Result[y,x] + B[y,x] + G[y,x]) * 0.33333333;
end;


{$I SimbaPlugin.inc}


end.
