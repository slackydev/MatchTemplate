unit FFTPACK4;
{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 For more info see: Copyright.txt
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
{$I header.inc}
interface

uses
  SysUtils, FFTPACK4_core, core;

type
  CFFT_FUNC = function(data, plan: TComplexArray; Inplace:Boolean=False): TComplexArray;

  TFFTPACK = record
    MaxThreads: Int32;
    MinThreadingSize: Int32;

    procedure Init(AMaxThreads:Int32);
    procedure Free();

    function OptimalDFTSize(target: Int32): Int32;

    function InitFFT(n: Int32): TComplexArray;
    function FFT(a, wsave: TComplexArray; Inplace: Boolean=False): TComplexArray;
    function IFFT(a, wsave: TComplexArray; Inplace: Boolean=False): TComplexArray;

    function InitRFFT(n: Int32): TRealArray;
    function RFFT(a, wsave: TRealArray; Inplace: Boolean=False): TRealArray;
    function IRFFT(a, wsave: TRealArray; Inplace: Boolean=False): TRealArray;

    function FFT2MT(m: T2DComplexArray; Inverse: Boolean): T2DComplexArray;
    function FFT2(m: T2DComplexArray): T2DComplexArray;
    function IFFT2(m: T2DComplexArray): T2DComplexArray;
  end;

const
  MIN_THREDING_SZ = 333*333;

var
  FFTPACK: TFFTPACK;

implementation

uses
  Math, matrix, threading;


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


{$DEFINE CPLX_BUFFSZ := 2*n + 15}
{$DEFINE REAL_BUFFSZ := 2*n + 15}


procedure TFFTPACK.Init(AMaxThreads: Int32);
begin
  Self.MaxThreads := AMaxThreads;
  Self.MinThreadingSize := MIN_THREDING_SZ;
end;

procedure TFFTPACK.Free();
begin
  Self.MaxThreads := 0;
  Self.MinThreadingSize := 0;
end;


// --------------------------------------------------------------------------------
// Compute the optimal size for FFT

function TFFTPACK.OptimalDFTSize(target: Int32): Int32;
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

function TFFTPACK.InitFFT(n: Int32): TComplexArray;
begin
  SetLength(Result, CPLX_BUFFSZ);
  cffti(n, @Result[0]);
end;

function TFFTPACK.FFT(a, wsave: TComplexArray; Inplace:Boolean=False): TComplexArray;
var n: Int32;
begin
  if Inplace then Result := a
  else            Result := Copy(a);
  n := Length(a); 
  Assert(Length(wsave) = CPLX_BUFFSZ, Format('Invalid work array for fft size (a: %d, w: %d)',[CPLX_BUFFSZ, Length(wsave)]));
  cfftf(n, @Result[0], @wsave[0]);
end;

function TFFTPACK.IFFT(a, wsave: TComplexArray; Inplace:Boolean=False): TComplexArray;
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

// --------------------------------------------------------------------------------
// real 2 real FFT

function TFFTPACK.InitRFFT(n: Int32): TRealArray;
begin
  SetLength(Result, REAL_BUFFSZ);
  rffti(n, @Result[0]);
end;

function TFFTPACK.RFFT(a, wsave: TRealArray; Inplace:Boolean=False): TRealArray;
var n: Int32;
begin
  if Inplace then Result := a
  else            Result := Copy(a);
  n := Length(a);
  Assert(Length(wsave) = REAL_BUFFSZ, Format('Invalid work array for fft size (a: %d, w: %d)',[REAL_BUFFSZ, Length(wsave)]));
  rfftf(n, @Result[0], @wsave[0]);
end;

function TFFTPACK.IRFFT(a, wsave: TRealArray; Inplace:Boolean=False): TRealArray;
var
  n: Int32;
  f: Single;
begin
  if Inplace then Result := a
  else            Result := Copy(a);
  n := Length(a);
  Assert(Length(wsave) = REAL_BUFFSZ, Format('Invalid work array for fft size (a: %d, w: %d)',[REAL_BUFFSZ, Length(wsave)]));
  rfftb(n, @Result[0], @wsave[0]);

  f := 1.0 / n;
  for n:=0 to High(a) do Result[n] *= f;
end;


// --------------------------------------------------------------------------------
// 2d complex fft (supports threading)

procedure __fft2_thread(params: PParamArray);
var
  y: Int32;
  data: T2DComplexArray;
  plan: TComplexArray;
begin
  data := T2DComplexArray(Params^[0]^);
  plan := Copy(TComplexArray(Params^[1]^)); //copy plan/workbase as it's also a buffer

  {if not inverse}
  if not PBoolean(Params^[2])^ then
    for y:=PBox(Params^[3])^.y1 to PBox(Params^[3])^.y2 do
      FFTPACK.FFT(data[y], plan, True)
  {if inverse}
  else
    for y:=PBox(Params^[3])^.y1 to PBox(Params^[3])^.y2 do
      FFTPACK.IFFT(data[y], plan, True);
end;

function TFFTPACK.FFT2MT(m: T2DComplexArray; Inverse: Boolean): T2DComplexArray;
var
  W,H,tc: Int32;
  plan: TComplexArray;
  ThreadPool: TThreadPool;
begin
  tc := Self.MaxThreads;
  ThreadPool := TThreadPool.Create(tc);

  Size(m, W,H);
  plan := InitFFT(W);
  ThreadPool.MatrixFunc(@__fft2_thread, [@m, @plan, @inverse], W,H, tc, Self.MinThreadingSize);

  m := Rot90(m);
  Size(m, W,H);
  plan := InitFFT(W);
  ThreadPool.MatrixFunc(@__fft2_thread, [@m, @plan, @inverse], W,H, tc, Self.MinThreadingSize);

  Result := Rot90(m);
  ThreadPool.Free();
end;


function TFFTPACK.FFT2(m: T2DComplexArray): T2DComplexArray;
begin
  if Length(m) = 0 then Exit;
  Result := FFT2MT(m, False);
end;

function TFFTPACK.IFFT2(m: T2DComplexArray): T2DComplexArray;
begin
  if Length(m) = 0 then Exit;
  Result := FFT2MT(m, True);
end;


// ----------------------------------------------------------------------------
// Initialize unit

initialization
  FFTPACK.Init(GetSystemThreadCount() div 2);


end.
