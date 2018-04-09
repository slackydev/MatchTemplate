unit core;
{==============================================================================]
  Copyright (c) 2018, Jarl `slacky` Holta
  Project: libfft
  License: GNU Lesser GPL (http://www.gnu.org/licenses/lgpl.html)
[==============================================================================}
{$I header.inc}

interface

uses
  SysUtils;
  
type
  PParamArray = ^TParamArray;
  TParamArray = array[Word] of Pointer;

  TReal = Single; //dont change me!
  Complex = packed record
    re, im: TReal;
  end;
  TComplexArray   = array of Complex;
  T2DComplexArray = array of TComplexArray;

  TRealArray   = array of TReal;
  T2DRealArray = array of array of TReal;

  TSingleArray   = array of Single;
  T2DSingleArray = array of array of Single;

  TDoubleArray   = array of Double;
  T2DDoubleArray = array of array of Double;

  TIntArray   = array of Int32;
  T2DIntArray = array of array of Int32;

  PBox = ^TBox;
  TBox = record X1,Y1,X2,Y2: Int32; end;

function NextPow2(n: Int32): Int32; inline;
function ParamArray(arr: array of Pointer): TParamArray;
function Box(X1,Y1,X2,Y2: Int32): TBox; inline;
function MarkTime(): Double;

//-----------------------------------------------------------------------
implementation

uses
  {$IFDEF WINDOWS}Windows{$ELSE}BaseUnix, Unix{$ENDIF};

function NextPow2(n: Int32): Int32;
begin
  n := n - 1;
  n := n or (n shr 1);
  n := n or (n shr 2);
  n := n or (n shr 4);
  n := n or (n shr 8);
  n := n or (n shr 16);
  n := n or (n shr 32);
  Result := n + 1;
end; 

function MarkTime(): Double;
var 
  frequency,count: Int64;
  {$IFDEF UNIX} 
  TV:TTimeVal; TZ:PTimeZone;
  {$ENDIF}
begin
  {$IFDEF WINDOWS}
  QueryPerformanceFrequency(frequency);
  QueryPerformanceCounter(count);
  Result := count / frequency * 1000;
  {$ELSE}
  TZ := nil;
  fpGetTimeOfDay(@TV, TZ);
  count := Int64(TV.tv_sec) * 1000000 + Int64(TV.tv_usec);
  Result := count / 1000;
  {$ENDIF} 
end;



function ParamArray(arr:Array of Pointer): TParamArray;
var i:Int32;
begin
  for i:=0 to High(arr) do Result[i] := arr[i];
end;

function Box(X1,Y1,X2,Y2: Int32): TBox;
begin
  Result.X1 := X1;
  Result.Y1 := Y1;
  Result.X2 := X2;
  Result.Y2 := Y2;
end;


end.
