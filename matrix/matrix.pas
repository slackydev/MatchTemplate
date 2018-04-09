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

procedure Size(a: T2DIntArray; out W,H: Int32); callconv
procedure Size(a: T2DSingleArray; out W,H: Int32); callconv
procedure Size(a: T2DDoubleArray; out W,H: Int32); callconv
procedure Size(a: T2DComplexArray; out W,H: Int32); callconv

procedure SplitRGB(img: T2DIntArray; out R,G,B: T2DRealArray); callconv
function IntToSingle(a: T2DIntArray): T2DSingleArray; callconv
function IntToDouble(a: T2DIntArray): T2DDoubleArray; callconv

function Pad2(a: T2DIntArray; M,N: Int32; init:Int32=0): T2DIntArray; callconv
function Pad2(a: T2DSingleArray; M,N: Int32; init:Single=0): T2DSingleArray; callconv
function Pad2(a: T2DDoubleArray; M,N: Int32; init:Double=0): T2DDoubleArray; callconv

function ShiftResize(a: T2DIntArray; M,N: Int32): T2DIntArray; callconv
function ShiftResize(a: T2DSingleArray; M,N: Int32): T2DSingleArray; callconv
function ShiftResize(a: T2DDoubleArray; M,N: Int32): T2DDoubleArray; callconv

function LocalSum(a: T2DIntArray; th,tw: Int32; padval: Int32): T2DIntArray; callconv
function LocalSum(a: T2DSingleArray; th,tw: Int32; padval: Single): T2DSingleArray; callconv
function LocalSum(a: T2DDoubleArray; th,tw: Int32; padval: Double): T2DDoubleArray; callconv

function Rot90(a: T2DIntArray): T2DIntArray; callconv
function Rot90(a: T2DSingleArray): T2DSingleArray; callconv
function Rot90(a: T2DDoubleArray): T2DDoubleArray; callconv
function Rot90(a: T2DComplexArray): T2DComplexArray; callconv

function Mean(a: T2DIntArray): Double; callconv
function Mean(a: T2DSingleArray): Double; callconv
function Mean(a: T2DDoubleArray): Double; callconv
procedure MeanStdev(a: T2DIntArray; out Mean, Stdev: Double); callconv
procedure MeanStdev(a: T2DSingleArray; out Mean, Stdev: Single); callconv
procedure MeanStdev(a: T2DDoubleArray; out Mean, Stdev: Double); callconv
procedure MinMax(a: T2DIntArray; out vMin,vMax: Int32); callconv
procedure MinMax(a: T2DSingleArray; out vMin,vMax: Single); callconv
procedure MinMax(a: T2DDoubleArray; out vMin,vMax: Double); callconv

function Squared(a: T2DIntArray): T2DIntArray; callconv
function Squared(a: T2DSingleArray): T2DSingleArray; callconv
function Squared(a: T2DDoubleArray): T2DDoubleArray; callconv
function SqrToDouble(a: T2DIntArray): T2DDoubleArray; callconv
function SqrToDouble(a: T2DSingleArray): T2DDoubleArray; callconv
function SqrToDouble(a: T2DDoubleArray): T2DDoubleArray; callconv

//-----------------------------------------------------------------------
implementation

{$I size.inc}

procedure SplitRGB(img: T2DIntArray; out R,G,B: T2DRealArray); callconv
var W,H,x,y:Int32;
begin
  Size(img, W,H);
  SetLength(R, H,W);
  SetLength(G, H,W);
  SetLength(B, H,W);
  for y:=0 to H-1 do
    for x:=0 to W-1 do begin
      R[y,x] := img[y,x] and $FF;
      G[y,x] := img[y,x] shr 8 and $FF;
      B[y,x] := img[y,x] shr 16 and $FF;
    end;
end;


function IntToSingle(a: T2DIntArray): T2DSingleArray; callconv
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H,W);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[i,j] := a[i,j];
end;

function IntToDouble(a: T2DIntArray): T2DDoubleArray; callconv
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
{$I rot90.inc}
{$I stats.inc}
{$I sqr.inc}


end.
