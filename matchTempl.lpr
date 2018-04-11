library matchTempl;
{==============================================================================]
  Copyright (c) 2018, Jarl `slacky` Holta
  Project: libfft
  License: GNU Lesser GPL (http://www.gnu.org/licenses/lgpl.html)
[==============================================================================}
{$I header.inc}

uses
  SysUtils, Math, 
  core, matrix, threading, fftpack, FFTW3;

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

function CrossCorrChannel(img, sub: T2DRealArray): T2DRealArray; callconv
var
  x,y,H,W,kH,kW: Int32;
  a,b: T2DComplexArray;
begin
  Size(img, W,H);
  Size(sub, kW,kH);

  if FFTW3.HAS_FFTW then
  begin
    SetLength(img, FFTW3.OptimalDFTSize(H), FFTW3.OptimalEvenDFTSize(W));
    SetLength(sub, FFTW3.OptimalDFTSize(H), FFTW3.OptimalEvenDFTSize(W));
    a := FFTW3.FFT2_R2C(img);
    b := FFTW3.FFT2_R2C(sub);
  end else
  begin
    SetLength(a, FFTPACK.OptimalDFTSize(H), FFTPACK.OptimalDFTSize(W));
    SetLength(b, FFTPACK.OptimalDFTSize(H), FFTPACK.OptimalDFTSize(W));
    for y:=0 to H-1  do for x:=0 to W-1  do a[y,x].re := img[y,x];
    for y:=0 to kH-1 do for x:=0 to kW-1 do b[y,x].re := sub[y,x];
    a := FFT2MT(a, False);
    b := FFT2MT(b, False);
  end;

  b := MulSpectrumConj(a,b);

  if FFTW3.HAS_FFTW then
  begin
    Result := FFTW3.FFT2_C2R(b);
    SetLength(Result, H,W);
    Exit;
  end else
  begin
    b := FFT2MT(b, True);
    SetLength(Result, H,W);
    for y:=0 to H-1 do for x:=0 to W-1 do Result[y,x] := b[y,x].re;
  end;
end;


function MatchTemplateFS(constref a,t: T2DRealArray): T2DRealArray; callconv
var
  x,y,tw,th,aw,ah: Int32;
  invSize, std_t, mean_t, sigma_t, denom, numer: Single;
  xcorr: T2DRealArray;
  asum, ls_diff: Double;
  sum, sum2: T2DDoubleArray;

  time: Double;

  function RegionSum(Y,X, M,N: Int32; Sums: T2DDoubleArray): Double; inline;
  begin
    Result := Sums[Y,X] - Sums[Y,X+N] - Sums[Y+M,X] + Sums[Y+M,X+N];
  end;
begin
  time := MarkTime();
  xcorr := CrossCorrChannel(a,t);
  WriteLn(Format('XCORR:  %.3f ms', [MarkTime() - time]));

  Size(t, tw,th);
  invSize := 1.0 / (tw*th);
  MeanStdev(t, mean_t, std_t);
  sigma_t := Sqrt(tw*th-1) * std_t;

  if sigma_t <> 0 then
  begin
    time := MarkTime();
    sum := SumsPd(a, sum2);
    WriteLn(Format('Sums:  %.3f ms', [MarkTime() - time]));

    time := MarkTime();
    Size(sum, aw,ah);
    SetLength(Result, ah,aw);
    for y:=0 to ah-th-1 do
      for x:=0 to aw-tw-1 do
      begin
        asum    := RegionSum(Y,X, th,tw, sum);
        ls_diff := RegionSum(Y,X, th,tw, sum2) - (Sqr(asum) * invSize);
        if ls_diff > 0 then
        begin
          denom := sigma_t * Sqrt(ls_diff);
          numer := xcorr[y,x] - asum * mean_t;
          Result[y,x] := numer / denom;
        end;
      end;
    WriteLn(Format('Final: %.3f ms', [MarkTime() - time]));
  end;

  Size(a, aw,ah);
  SetLength(Result, ah-th+1,aw-tw+1);
end;


function MatchTemplateSS(constref a,t: T2DRealArray): T2DRealArray; callconv
var
  x,y,tw,th,aw,ah: Int32;
  invSize, std_t, mean_t, sigma_t, denom, numer: Single;
  sum, xcorr: T2DRealArray;
  asum, ls_diff: Double;
  sum2: T2DDoubleArray;

  time: Double;
begin
  Size(t, tw,th);

  time := MarkTime();
  xcorr := CrossCorrChannel(ShiftResize(a,th-1,tw-1), t);
  WriteLn(Format('XCORR:  %.3f ms', [MarkTime() - time]));

  invSize := 1.0 / (tw*th);
  MeanStdev(t, mean_t, std_t);
  sigma_t := Sqrt(tw*th-1) * std_t;

  if sigma_t <> 0 then
  begin
    time := MarkTime();
    sum  := LocalSum(a,th,tw);
    sum2 := LocalSum(SqrToDouble(a),th,tw);
    WriteLn(Format('Sums:  %.3f ms', [MarkTime() - time]));

    time := MarkTime();
    Size(sum, aw,ah);
    SetLength(Result, ah,aw);
    Dec(th); Dec(tw);
    for y:=th to ah-th-1 do
      for x:=tw to aw-tw-1 do
      begin
        ls_diff := sum2[y,x] - (Sqr(sum[y,x]) * invSize);
        if ls_diff > 0 then
        begin
          denom := sigma_t * Sqrt(ls_diff);
          numer := xcorr[y,x] - sum[y,x] * mean_t;
          Result[y-th,x-tw] := numer / denom;
        end;
      end;
    WriteLn(Format('Final: %.3f ms', [MarkTime() - time]));
  end;

  Size(a, aw,ah);
  SetLength(Result, ah-th,aw-tw);
end;


function MatchTemplate(constref img, sub: T2DIntArray; maxThreads:Int32=-1): T2DRealArray; callconv
var
  R1,G1,B1, R2,G2,B2: T2DRealArray;
  x,y,W,H: Int32;
begin
  Size(img, W,H);
  SplitRGB(img, R1,G1,B1);
  SplitRGB(sub, R2,G2,B2);

  if FFTW3.HAS_FFTW then
  begin
    if (maxThreads <= -1) and (W*H < FFTW3.MIN_THREDING_SZ) then
      fftw_plan_with_nthreads(1)
    else if (maxThreads <= -1) then
      fftw_plan_with_nthreads(GetSystemThreadCount() div 2)
    else
      fftw_plan_with_nthreads(Min(maxThreads,GetSystemThreadCount()));
  end;

  if W*H <= 700*700 then //unsure if this is a truely safe number.
  begin
    // This can only handle input sizes up to.. As it will get inaccurate due to
    // floatingpoint roundoff when it sums up some very large numbers.
    R1 := MatchTemplateFS(R1,R2);  R2 := nil;
    G1 := MatchTemplateFS(G1,G2);  G2 := nil;
    B1 := MatchTemplateFS(B1,B2);  B2 := nil;
  end else
  begin
    R1 := MatchTemplateSS(R1,R2);  R2 := nil;
    G1 := MatchTemplateSS(G1,G2);  G2 := nil;
    B1 := MatchTemplateSS(B1,B2);  B2 := nil;
  end;

  // This is imperfect, but my guess is for this usage, it's OK
  Result := R1;
  Size(Result, W,H);
  for y:=0 to H-1 do
    for x:=0 to W-1 do
      Result[y,x] := (R1[y,x] + B1[y,x] + G1[y,x]) * 0.33333333;
end;


{$I SimbaPlugin.inc}


end.
