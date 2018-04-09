{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 Resize to M,N, padding happens before the items, rather than after
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
function Squared(a: T2DIntArray): T2DIntArray; callconv
var W,H,i,j:Int32;
begin
  Size(a, W,H);
  SetLength(Result, H,W);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[i,j] := Sqr(a[i,j]);
end;  

function Squared(a: T2DSingleArray): T2DSingleArray; callconv
var W,H,i,j:Int32;
begin
  Size(a, W,H);
  SetLength(Result, H,W);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[i,j] := Sqr(a[i,j]);
end; 

function Squared(a: T2DDoubleArray): T2DDoubleArray; callconv
var W,H,i,j:Int32;
begin
  Size(a, W,H);
  SetLength(Result, H,W);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[i,j] := Sqr(a[i,j]);
end; 

// ----------------------------------------------------------
// Square with Result being Double

function SqrToDouble(a: T2DIntArray): T2DDoubleArray; callconv
var W,H,i,j:Int32;
begin
  Size(a, W,H);
  SetLength(Result, H,W);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[i,j] := Sqr(a[i,j]);
end; 

function SqrToDouble(a: T2DSingleArray): T2DDoubleArray; callconv
var W,H,i,j:Int32;
begin
  Size(a, W,H);
  SetLength(Result, H,W);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[i,j] := Sqr(a[i,j]);
end; 

function SqrToDouble(a: T2DDoubleArray): T2DDoubleArray; callconv
var W,H,i,j:Int32;
begin
  Size(a, W,H);
  SetLength(Result, H,W);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[i,j] := Sqr(a[i,j]);
end; 