unit FFTW3;
{==============================================================================]
  fftpack.pas: A set of FFT routines in FPC
  Algorithmically based on Fortran-77 FFTPACK by Paul N. Swarztrauber (V4, 1985)

  Ported to Free Pascal by Jarl `slacky` Holta
[==============================================================================}
{$I header.inc}
interface

uses
  SysUtils, dynlibs, core;

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
function OptimalEvenDFTSize(target: Int32): Int32;

function FFT2(m: T2DComplexArray): T2DComplexArray; callconv
function IFFT2(m: T2DComplexArray): T2DComplexArray; callconv
procedure FFT_RGB(var R,G,B: T2DComplexArray); callconv
function FFT2_R2C(m: T2DRealArray): T2DComplexArray; callconv
function FFT2_C2R(m: T2DComplexArray): T2DRealArray; callconv

const
  FFTW_FORWARD        = -1;
  FFTW_BACKWARD        = 1;

  FFTW_MEASURE         = 0;
  FFTW_DESTROY_INPUT   = 1;   {1U << 0}
  FFTW_UNALIGNED       = 2;   {1U << 1}
  FFTW_CONSERVE_MEMORY = 4;   {1U << 2}
  FFTW_EXHAUSTIVE      = 8;   {1U << 3} {NO_EXHAUSTIVE is default }
  FFTW_PRESERVE_INPUT  = 16;  {1U << 4} {cancels FFTW_DESTROY_INPUT}
  FFTW_PATIENT         = 32;  {1U << 5} {IMPATIENT is default }
  FFTW_ESTIMATE        = 64;  {1U << 6}   

type
  FFTW_PLAN  = type Pointer;

  FFTW_PLAN_REC = record
    Plan: FFTW_PLAN;
    Dispose: Boolean;
  end;

  FFTW_PLAN_LIST_2D = array of record
    Plan: FFTW_PLAN_REC;
    IsRunning: Boolean;
    W,H: Int32;
  end;

var
  FFTW: TLibHandle;
  HAS_FFTW: Boolean;

  fftw_init_threads:       function(): Int32; cdecl;
  fftw_plan_with_nthreads: procedure(nthreads: Int32); cdecl;
  fftw_cleanup_threads:    procedure(); cdecl;

  fftw_malloc:    function(n: Int32): Pointer; cdecl;
  fftw_free:      procedure(p: Pointer); cdecl;
  fftw_plan_dft1: function(n: Int32; inData, outData: PComplex; sign: Int32; flags: UInt32): FFTW_PLAN; cdecl;
  fftw_plan_dft2: function(n0, n1: Int32; inData, outData: PComplex; sign: Int32; flags: UInt32): FFTW_PLAN; cdecl;
  fftw_plan_dft3: function(n0, n1, n2: Int32; inData, outData: PComplex; sign: Int32; flags: UInt32): FFTW_PLAN; cdecl;
  fftw_plan_dft2_r2c: function(n0, n1: Int32; inData: PReal; outData: PComplex; flags: UInt32): FFTW_PLAN; cdecl;
  fftw_plan_dft2_c2r: function(n0, n1: Int32; inData: PComplex; outData: PReal; flags: UInt32): FFTW_PLAN; cdecl;
  fftw_exec:      procedure(plan: FFTW_PLAN); cdecl;
  fftw_exec_dft:  procedure(plan: FFTW_PLAN; inData, outData: PComplex); cdecl;

  fftw_free_plan: procedure(plan: FFTW_PLAN); cdecl;
  

implementation

uses
  Math, matrix, threading;

procedure LoadFFTW();
begin
  {$IF Defined(WIN32)}
  FFTW := LoadLibrary('libfftw3f-3_32.dll');
  {$ELSEIF Defined(WIN64)}
  FFTW := LoadLibrary('libfftw3f-3_64.dll');
  {$ELSEIF Defined(CPU386) and Defined(LINUX)}
  FFTW := LoadLibrary('libfftw3f-3_32.so');
  {$ELSEIF Defined(CPUX86_64) and Defined(LINUX)}
  FFTW := LoadLibrary('libfftw3f-3_64.so');
  {$ENDIF}
  HAS_FFTW := FFTW <> 0;

  if not HAS_FFTW then
  begin
    {$IF Defined(Windows)}
    FFTW := LoadLibrary('ibfftw3f-3.dll');
    {$ELSEIF Defined(Linux)}
    FFTW := LoadLibrary('libfftw3f-3.so');
    {$ENDIF}
    HAS_FFTW := FFTW <> 0;
  end;


  if HAS_FFTW then
  begin
    WriteLn('Found FFTW3.3 - Loaded successfully');

    Pointer(fftw_init_threads)       := GetProcAddress(FFTW, 'fftwf_init_threads');
    Pointer(fftw_plan_with_nthreads) := GetProcAddress(FFTW, 'fftwf_plan_with_nthreads');
    Pointer(fftw_cleanup_threads)    := GetProcAddress(FFTW, 'fftwf_cleanup_threads');

    Pointer(fftw_malloc)        := GetProcAddress(FFTW, 'fftwf_malloc');
    Pointer(fftw_free)          := GetProcAddress(FFTW, 'fftwf_free');

    Pointer(fftw_plan_dft1)     := GetProcAddress(FFTW, 'fftwf_plan_dft_1d');
    Pointer(fftw_plan_dft2)     := GetProcAddress(FFTW, 'fftwf_plan_dft_2d');
    Pointer(fftw_plan_dft3)     := GetProcAddress(FFTW, 'fftwf_plan_dft_3d');
    Pointer(fftw_plan_dft2_r2c) := GetProcAddress(FFTW, 'fftwf_plan_dft_r2c_2d');
    Pointer(fftw_plan_dft2_c2r) := GetProcAddress(FFTW, 'fftwf_plan_dft_c2r_2d');
    Pointer(fftw_free_plan)     := GetProcAddress(FFTW, 'fftwf_destroy_plan');

    Pointer(fftw_exec)          := GetProcAddress(FFTW, 'fftwf_execute');
    Pointer(fftw_exec_dft)      := GetProcAddress(FFTW, 'fftwf_execute_dft');

    fftw_init_threads();
    fftw_plan_with_nthreads(GetSystemThreadCount() div 2);
  end else
    WriteLn('FFTW3.3 was not found');
end; 
  

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

function OptimalEvenDFTSize(target: Int32): Int32;
var
  n: Int32;
begin
  n := Target;
  while True do
  begin
    n := OptimalDFTSize(n);
    if n and 1 = 0 then Exit(n);
    Inc(n);
  end;
end;

(*
function GetPlan2D(W,H: Int32): FFTW_PLAN_REC;
var
  i: Int32;
begin
  for i:=0 to High(FFTW_PLAN_2D) do
    if (FFTW_PLAN_2D[i].W = W) and (FFTW_PLAN_2D[i].H = H) then
    begin
      if FFTW_PLAN_2D[i].IsRunning then
      begin
        Result.Plan := fftw_plan_dft2(H,W, _data, _res, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        Exit;
      end;
    end;
end;
*)

// --------------------------------------------------------------------------------
// complex 2 complex 2d fft

function FFT2(m: T2DComplexArray): T2DComplexArray; callconv
var
  y,w,h,row,sz: Int32;
  _data, _res: PComplex;  
  plan: FFTW_PLAN;
  t: Double;
begin
  Size(m, W,H);
  if W = 0 then Exit;

  row := SizeOf(Complex) * W;
  sz  := row * H;

  _data := fftw_malloc(sz);
  _res  := fftw_malloc(sz);
  plan  := fftw_plan_dft2(H,W, _data, _res, FFTW_FORWARD, FFTW_ESTIMATE);

  for y:=0 to h-1 do Move(m[y,0], _data[y*w], row);
  fftw_exec(plan);
  for y:=0 to h-1 do Move(_res[y*w], m[y,0], row);

  fftw_free(_data);
  fftw_free(_res);
  fftw_free_plan(plan);

  Result := m;
end;

function IFFT2(m: T2DComplexArray): T2DComplexArray; callconv
var
  x,y,row,sz,w,h: Int32;
  f: Single;   
  _data, _res: PComplex;  
  plan: FFTW_PLAN;
begin
  Size(m, W,H);
  if W = 0 then Exit;
  
  row := SizeOf(Complex) * W;
  sz  := row * H;

  _data := fftw_malloc(sz);
  _res  := fftw_malloc(sz);
  plan  := fftw_plan_dft2(H,W, _data, _res, FFTW_BACKWARD, FFTW_ESTIMATE);
  
  for y:=0 to h-1 do Move(m[y,0], _data[y*w], row);
  fftw_exec(plan);
  for y:=0 to h-1 do Move(_res[y*w], m[y,0], row);

  fftw_free(_data);
  fftw_free(_res);
  fftw_free_plan(plan);

  Result := m;
  f := 1.0 / (W*H);
  for y:=0 to H-1 do
    for x:=0 to W-1 do
    begin
      Result[y,x].re := Result[y,x].re * f;
      Result[y,x].im := Result[y,x].im * f;
    end;
end;

// --------------------------------------------------------------------------------
// complex 2 complex 3d fft

procedure FFT_RGB(var R,G,B: T2DComplexArray); callconv
var
  y,w,h,row,sz: Int32;
  _data, _res: PComplex;
  plan: FFTW_PLAN;
  t: Double;
begin
  Size(R, W,H);
  if W = 0 then Exit;

  row := SizeOf(Complex) * W;
  sz  := row * H;

  _data := fftw_malloc(3*sz);
  _res  := fftw_malloc(3*sz);
  plan  := fftw_plan_dft3(3,H,W, _data, _res, FFTW_FORWARD, FFTW_ESTIMATE);

  for y:=0 to h-1 do Move(R[y,0], _data[1*y*w], row);
  for y:=0 to h-1 do Move(G[y,0], _data[2*y*w], row);
  for y:=0 to h-1 do Move(B[y,0], _data[3*y*w], row);
  fftw_exec(plan);
  for y:=0 to h-1 do Move(_res[1*y*w], R[y,0], row);
  for y:=0 to h-1 do Move(_res[2*y*w], G[y,0], row);
  for y:=0 to h-1 do Move(_res[3*y*w], B[y,0], row);

  fftw_free(_data);
  fftw_free(_res);
  fftw_free_plan(plan);
end;

// --------------------------------------------------------------------------------
// real 2 complex 2d fft     /// NOT WORKING | UNFINISHED

function FFT2_R2C(m: T2DRealArray): T2DComplexArray; callconv
var
  y,w,h,newW: Int32;
  input: PReal;
  output: PComplex;
  plan: FFTW_PLAN;
begin
  Size(m, W,H);
  if W = 0 then Exit;
  SetLength(Result, H,W);

  newW := W div 2+1;

  input  := fftw_malloc(SizeOf(TReal)   * H*W);
  output := fftw_malloc(SizeOf(Complex) * H*newW);
  plan := fftw_plan_dft2_r2c(H,W, input, output, FFTW_ESTIMATE);

  for y:=0 to H-1 do Move(m[y,0], input[y*w], SizeOf(TReal) * W);
  fftw_exec(plan);
  for y:=0 to H-1 do Move(output[y*newW], Result[y,0], SizeOf(Complex) * newW);

  fftw_free(input);
  fftw_free(output);
  fftw_free_plan(plan);
end;

function FFT2_C2R(m: T2DComplexArray): T2DRealArray; callconv
var
  x,y,w,h,newW: Int32;
  input: PComplex;
  output: PReal;
  plan: FFTW_PLAN;
  f: Single;
begin
  Size(m, W,H);
  if W = 0 then Exit;
  SetLength(Result, H,W);

  newW := W div 2+1;

  input  := fftw_malloc(SizeOf(Complex) * H*newW);
  output := fftw_malloc(SizeOf(TReal)   * H*W);
  plan  := fftw_plan_dft2_c2r(H,W, input, output, FFTW_ESTIMATE);

  for y:=0 to h-1 do Move(m[y,0], input[y*newW], SizeOf(Complex) * newW);
  fftw_exec(plan);
  for y:=0 to h-1 do Move(output[y*w], Result[y,0], SizeOf(TReal) * W);

  fftw_free(input);
  fftw_free(output);
  fftw_free_plan(plan);

  f := 1.0 / (W*H);
  for y:=0 to H-1 do
    for x:=0 to W-1 do Result[y,x] := Result[y,x] * f;
end;

initialization
  LoadFFTW();

end.
