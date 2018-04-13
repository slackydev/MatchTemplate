unit matrix;
{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2013, Jarl K. <Slacky> Holta || http://github.com/WarPie
 All rights reserved.
 For more info see: Copyright.txt
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
{$I header.inc}

interface

uses
  SysUtils, core;

procedure Size(a: T2DIntArray; out W,H: Int32);
procedure Size(a: T2DI64Array; out W,H: Int32);
procedure Size(a: T2DSingleArray; out W,H: Int32);
procedure Size(a: T2DDoubleArray; out W,H: Int32);
procedure Size(a: T2DComplexArray; out W,H: Int32);

procedure SplitRGB(img: T2DIntArray; out R,G,B: T2DRealArray);

function IntToSingle(a: T2DIntArray): T2DSingleArray;
function IntToDouble(a: T2DIntArray): T2DDoubleArray;
function IntToReal(a: T2DIntArray): T2DRealArray;

function Pad2(a: T2DIntArray; M,N: Int32): T2DIntArray;
function Pad2(a: T2DI64Array; M,N: Int32): T2DI64Array;
function Pad2(a: T2DSingleArray; M,N: Int32): T2DSingleArray;
function Pad2(a: T2DDoubleArray; M,N: Int32): T2DDoubleArray;

function ShiftResize(a: T2DIntArray; M,N: Int32): T2DIntArray;
function ShiftResize(a: T2DSingleArray; M,N: Int32): T2DSingleArray;
function ShiftResize(a: T2DDoubleArray; M,N: Int32): T2DDoubleArray;

function LocalSum(a: T2DIntArray; th,tw: Int32): T2DIntArray;
function LocalSum(a: T2DI64Array; th,tw: Int32): T2DI64Array;
function LocalSum(a: T2DSingleArray; th,tw: Int32): T2DSingleArray;
function LocalSum(a: T2DDoubleArray; th,tw: Int32): T2DDoubleArray;

function Sums(a: T2DIntArray; out Square: T2DI64Array): T2DIntArray;
function Sums(a: T2DSingleArray; out Square: T2DDoubleArray): T2DDoubleArray;
function Sums(a: T2DDoubleArray; out Square: T2DDoubleArray): T2DDoubleArray;

function SumsPd(a: T2DIntArray; out Square: T2DI64Array): T2DI64Array;
function SumsPd(a: T2DSingleArray; out Square: T2DDoubleArray): T2DDoubleArray;
function SumsPd(a: T2DDoubleArray; out Square: T2DDoubleArray): T2DDoubleArray;

function Rot90(a: T2DIntArray): T2DIntArray;
function Rot90(a: T2DSingleArray): T2DSingleArray;
function Rot90(a: T2DDoubleArray): T2DDoubleArray;
function Rot90(a: T2DComplexArray): T2DComplexArray;

function Mean(a: T2DIntArray): Double;
function Mean(a: T2DSingleArray): Double;
function Mean(a: T2DDoubleArray): Double;
procedure MeanStdev(a: T2DIntArray; out Mean, Stdev: Double);
procedure MeanStdev(a: T2DSingleArray; out Mean, Stdev: Double);
procedure MeanStdev(a: T2DDoubleArray; out Mean, Stdev: Double);
procedure MinMax(a: T2DIntArray; out vMin,vMax: Int32);
procedure MinMax(a: T2DSingleArray; out vMin,vMax: Single);
procedure MinMax(a: T2DDoubleArray; out vMin,vMax: Double);

function Squared(a: T2DIntArray): T2DIntArray;
function Squared(a: T2DSingleArray): T2DSingleArray;
function Squared(a: T2DDoubleArray): T2DDoubleArray;
function SqrToDouble(a: T2DIntArray): T2DDoubleArray;
function SqrToDouble(a: T2DSingleArray): T2DDoubleArray;
function SqrToDouble(a: T2DDoubleArray): T2DDoubleArray;

//-----------------------------------------------------------------------
implementation

{$I size.inc}

procedure SplitRGB(img: T2DIntArray; out R,G,B: T2DRealArray);
var
  W,H,x,y:Int32;
begin
  Size(img, W,H);
  SetLength(R, H,W);
  SetLength(G, H,W);
  SetLength(B, H,W);
  for y:=0 to H-1 do
    for x:=0 to W-1 do begin
      R[y,x] := (img[y,x]{shr 00}and $FF);
      G[y,x] := (img[y,x] shr 08 and $FF);
      B[y,x] := (img[y,x] shr 16 and $FF);
    end;
end;


function IntToSingle(a: T2DIntArray): T2DSingleArray;
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H,W);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[i,j] := a[i,j];
end;

function IntToDouble(a: T2DIntArray): T2DDoubleArray;
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H,W);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[i,j] := a[i,j];
end;

function IntToReal(a: T2DIntArray): T2DRealArray;
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H,W);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[i,j] := a[i,j];
end;

{$I pad2.inc}
{$I shiftresize.inc}
{$I localsum.inc}
{$I sums.inc}
{$I rot90.inc}
{$I stats.inc}
{$I sqr.inc}


end.
