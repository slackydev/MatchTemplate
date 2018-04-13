{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 For more info see: Copyright.txt
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
function Rot90(a: T2DIntArray): T2DIntArray;
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, W,H);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[j,i] := a[i,j];
end;

function Rot90(a: T2DSingleArray): T2DSingleArray;
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, W,H);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[j,i] := a[i,j];
end;

function Rot90(a: T2DDoubleArray): T2DDoubleArray;
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, W,H);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[j,i] := a[i,j];
end;

function Rot90(a: T2DComplexArray): T2DComplexArray;
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, W,H);
  for i:=0 to H-1 do
    for j:=0 to W-1 do Result[j,i] := a[i,j];
end;
