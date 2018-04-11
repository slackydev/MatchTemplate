unit fftpack;
{==============================================================================]
  fftpack.pas: A set of FFT routines in FPC
  Algorithmically based on Fortran-77 FFTPACK by Paul N. Swarztrauber (V4, 1985)

  Ported to Free Pascal by Jarl `slacky` Holta
[==============================================================================}
{$I header.inc}
interface

uses
  SysUtils, fftpack_core, core;

type
  CFFT_FUNC = function(data, plan: TComplexArray; Inplace:LongBool=False): TComplexArray; callconv

const
  __OptimalDFT: array[0..168] of Int32 = (
    8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48,
    50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128,
    135, 144, 150, 160, 162, 180, 192, 200, 216, 225, 240, 243, 250,
    256, 270, 288, 300, 320, 324, 360, 375, 384, 400, 405, 432, 450,
    480, 486, 500, 512, 540, 576, 600, 625, 640, 648, 675, 720, 729,
    750, 768, 800, 810, 864, 900, 960, 972, 1000, 1024, 1080, 1125,
    1152, 1200, 1215, 1250, 1280, 1296, 1350, 1440, 1458, 1500, 1536,
    1600, 1620, 1728, 1800, 1875, 1920, 1944, 2000, 2025, 2048, 2160,
    2187, 2250, 2304, 2400, 2430, 2500, 2560, 2592, 2700, 2880, 2916,
    3000, 3072, 3125, 3200, 3240, 3375, 3456, 3600, 3645, 3750, 3840,
    3888, 4000, 4050, 4096, 4320, 4374, 4500, 4608, 4800, 4860, 5000,
    5120, 5184, 5400, 5625, 5760, 5832, 6000, 6075, 6144, 6250, 6400,
    6480, 6561, 6750, 6912, 7200, 7290, 7500, 7680, 7776, 8000, 8100,
    8192, 8640, 8748, 9000, 9216, 9375, 9600, 9720, 10000
  );

const
  MIN_THREDING_SZ = 333*333;


function OptimalDFTSize(target: Int32): Int32;

function FFTInit(n: Int32): TComplexArray; callconv
function FFT(a, wsave: TComplexArray; Inplace: LongBool = False): TComplexArray; callconv
function IFFT(a, wsave: TComplexArray; Inplace: LongBool = False): TComplexArray; callconv
function BFFT(a, wsave: TComplexArray; Inplace: LongBool = False): TComplexArray; callconv

function RFFTInit(n: Int32): TRealArray; callconv
function RFFT(a, wsave: TRealArray; Inplace: LongBool = False): TRealArray; callconv
function IRFFT(a, wsave: TRealArray; Inplace: LongBool = False): TRealArray; callconv
function BRFFT(a, wsave: TRealArray; Inplace: LongBool = False): TRealArray; callconv

function FFT2(m: T2DComplexArray; Inplace: LongBool = False): T2DComplexArray; callconv
function IFFT2(m: T2DComplexArray; Inplace: LongBool = False): T2DComplexArray; callconv
function FFT2MT(m: T2DComplexArray; Inverse: LongBool): T2DComplexArray;


implementation

uses
  Math, matrix, threading;

function OptimalDFTSize(target: Int32): Int32;
var
  n,match,quotient,p2,p5,p35: Int32;
begin
  if (target <= 6) then
    Exit(target);
  
  if NextPow2(target) = target then
    Exit(target);
  
  n := 0;
  if target <= __OptimalDFT[High(__OptimalDFT)] then
  begin
    while __OptimalDFT[n] < target do Inc(n);
    Exit(__OptimalDFT[n]);
  end;

  match := $7FFFFFFF;
  p5 := 1;
  while p5 < target do
  begin
    p35 := p5;
    while p35 < target do
    begin
      quotient := Ceil(target / p35);
      p2 := NextPow2(quotient);
      N := p2 * p35;

      if N = target then Exit(N);
      if N < match  then match := N;

      p35 *= 3;
      if p35 = target then Exit(p35)
    end;
    if p35 < match then
      match := p35;
    
    p5 *= 5;
    if p5 = target then
      Exit(p5);
  end;
  Result := Min(p5, match);
end;

// --------------------------------------------------------------------------------
// complex 2 complex FFT

{$DEFINE CPLX_BUFFSZ := 2*n + 15}
function FFTInit(n: Int32): TComplexArray; callconv
begin
  SetLength(Result, CPLX_BUFFSZ); //4*n + 15
  cffti(n, @Result[0]);
end;

function FFT(a, wsave: TComplexArray; Inplace: LongBool = False): TComplexArray; callconv
var n: Int32;
begin
  if Inplace then Result := a
  else            Result := Copy(a);
  n := Length(a); 
  Assert(Length(wsave) = CPLX_BUFFSZ, Format('Invalid work array for fft size (a: %d, w: %d)',[CPLX_BUFFSZ, Length(wsave)]));
  cfftf(n, @Result[0], @wsave[0]);
end;

function IFFT(a, wsave: TComplexArray; Inplace: LongBool = False): TComplexArray; callconv
var
  n: Int32;
  f: Single;
begin
  if Inplace then Result := a
  else            Result := Copy(a);
  n := Length(a);
  Assert(Length(wsave) = CPLX_BUFFSZ, Format('Invalid work array for fft size (a: %d, w: %d)',[CPLX_BUFFSZ, Length(wsave)]));
  cfftb(n, @Result[0], @wsave[0]);

  f := 1.0 / n;
  for n:=0 to High(a) do
  begin
    Result[n].re *= f;
    Result[n].im *= f;
  end;
end;

function BFFT(a, wsave: TComplexArray; Inplace: LongBool = False): TComplexArray; callconv
var
  n: Int32;
begin
  if Inplace then Result := a
  else            Result := Copy(a);
  n := Length(a);
  Assert(Length(wsave) = CPLX_BUFFSZ, Format('Invalid work array for fft size (a: %d, w: %d)',[CPLX_BUFFSZ, Length(wsave)]));
  cfftb(n, @Result[0], @wsave[0]);
end;

// --------------------------------------------------------------------------------
// real 2 real FFT

function RFFTInit(n: Int32): TRealArray; callconv
begin
  SetLength(Result, 2*n + 15);
  rffti(n, @Result[0]);
end;

function RFFT(a, wsave: TRealArray; Inplace: LongBool = False): TRealArray; callconv
var n: Int32;
begin
  if Inplace then Result := a
  else            Result := Copy(a);
  n := Length(a);
  Assert(Length(wsave) = (2*n + 15), Format('Invalid work array for fft size (a: %d, w: %d)',[2*n + 8, Length(wsave)]));
  rfftf(n, @Result[0], @wsave[0]);
end;

function IRFFT(a, wsave: TRealArray; Inplace: LongBool = False): TRealArray; callconv
var
  n: Int32;
  f: Single;
begin
  if Inplace then Result := a
  else            Result := Copy(a);
  n := Length(a);
  Assert(Length(wsave) = (2*n + 15), Format('Invalid work array for fft size (a: %d, w: %d)',[2*n + 8, Length(wsave)]));
  rfftb(n, @Result[0], @wsave[0]);
  
  f := 1.0 / n;
  for n:=0 to High(a) do Result[n] *= f;
end;

function BRFFT(a, wsave: TRealArray; Inplace: LongBool = False): TRealArray; callconv
var
  n: Int32;
begin
  if Inplace then Result := a
  else            Result := Copy(a);
  n := Length(a);
  Assert(Length(wsave) = (2*n + 15), Format('Invalid work array for fft size (a: %d, w: %d)',[2*n + 8, Length(wsave)]));
  rfftb(n, @Result[0], @wsave[0]);
end;


// --------------------------------------------------------------------------------
// 2d fft

// complex 2 complex
function FFT2(m: T2DComplexArray; Inplace: LongBool = False): T2DComplexArray; callconv
var
  x,y,w,h: Int32;
  buffer: T2DComplexArray;
  t,plan: TComplexArray;
begin
  H := Length(m);
  if H = 0 then Exit;
  W := Length(m[0]);

  if Inplace then Result := m
  else            SetLength(Result, H,W);

  SetLength(buffer, W,H);
  plan := FFTInit(W);
  for y:=0 to h-1 do
  begin
    t := fft(m[y],plan, Inplace);
    for x:=0 to w-1 do buffer[x,y] := t[x];
  end;

  plan := FFTInit(H);
  for x:=0 to w-1 do
  begin
    t := fft(buffer[x], plan, True);
    for y:=0 to h-1 do Result[y,x] := t[y];
  end;
end;

function IFFT2(m: T2DComplexArray; Inplace: LongBool = False): T2DComplexArray; callconv
var
  x,y,w,h: Int32;
  f: Single;
  buffer: T2DComplexArray;
  t,plan: TComplexArray;
begin
  H := Length(m);
  if H = 0 then Exit;
  W := Length(m[0]);

  if Inplace then Result := m
  else            SetLength(Result, H,W);

  SetLength(buffer, W,H);
  f := 1.0 / W;
  plan := FFTInit(W);
  for y:=0 to h-1 do
  begin
    t := bfft(m[y], plan, Inplace);
    for x:=0 to w-1 do
    begin
      buffer[x,y].re := t[x].re * f;
      buffer[x,y].im := t[x].im * f;
    end;
  end;

  f := 1.0 / H;
  plan := FFTInit(H);
  for x:=0 to w-1 do
  begin
    t := bfft(buffer[x], plan, True);
    for y:=0 to h-1 do
    begin
      Result[y,x].re := t[y].re * f;
      Result[y,x].im := t[y].im * f;
    end;
  end;
end;


// --------------------------------------------------------------------------------
// 2d complex fft - threaded

procedure __fft2_thread(params: PParamArray);
var
  y: Int32;
  data: T2DComplexArray;
  plan: TComplexArray;
  func: CFFT_FUNC;
begin
  data := T2DComplexArray(Params^[0]^);
  plan := Copy(TComplexArray(Params^[1]^)); //copy plan/workbase as it's also a buffer
  func := CFFT_FUNC(Params^[2]^);

  for y:=PBox(Params^[3])^.y1 to PBox(Params^[3])^.y2 do
    func(data[y], plan, True);
end;

function FFT2MT(m: T2DComplexArray; Inverse: LongBool): T2DComplexArray;
var
  W,H,tc: Int32;
  plan: TComplexArray;
  ThreadPool: TThreadPool;
  func: CFFT_FUNC;
begin
  func := @fft;
  if Inverse then func := @ifft;

  tc := GetSystemThreadCount div 2;
  ThreadPool := TThreadPool.Create(tc);

  Size(m, W,H);
  plan := FFTInit(W);
  ThreadPool.MatrixFunc(@__fft2_thread, [@m, @plan, @func], W,H, tc, MIN_THREDING_SZ);

  m := Rot90(m);
  Size(m, W,H);
  plan := FFTInit(W);
  ThreadPool.MatrixFunc(@__fft2_thread, [@m, @plan, @func], W,H, tc, MIN_THREDING_SZ);

  Result := Rot90(m);
  ThreadPool.Free();
end;


end.
