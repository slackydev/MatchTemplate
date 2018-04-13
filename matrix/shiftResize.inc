{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 Resize to M,N, padding happens before the items, rather than after
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
function ShiftResize(a: T2DIntArray; M,N: Int32): T2DIntArray;
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H+M, W+N);
  for i:=0 to H-1 do
    for j:=0 to W-1 do
      Result[i+M,j+N] := a[i,j];
end;   

function ShiftResize(a: T2DSingleArray; M,N: Int32): T2DSingleArray;
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H+M, W+N);
  for i:=0 to H-1 do
    for j:=0 to W-1 do
      Result[i+M,j+N] := a[i,j];
end;   

function ShiftResize(a: T2DDoubleArray; M,N: Int32): T2DDoubleArray;
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H+M, W+N);
  for i:=0 to H-1 do
    for j:=0 to W-1 do
      Result[i+M,j+N] := a[i,j];
end;   
