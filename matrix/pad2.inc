{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 For more info see: Copyright.txt
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
function Pad2(a: T2DIntArray; M,N: Int32; init:Int32=0): T2DIntArray; callconv
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H+M*2, W+N*2);
  if init <> 0 then
    for i:=0 to H+M*2-1 do
      for j:=0 to W+N*2-1 do
        Result[i,j] := init;

  for i:=0 to H-1 do
    Move(a[i,0], Result[i+M,N], W*SizeOf(Int32));
end;

function Pad2(a: T2DSingleArray; M,N: Int32; init:Single=0): T2DSingleArray; callconv
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H+M*2, W+N*2);
  if init <> 0 then
    for i:=0 to H+M*2-1 do
      for j:=0 to W+N*2-1 do
        Result[i,j] := init;

  for i:=0 to H-1 do
    Move(a[i,0], Result[i+M,N], W*SizeOf(Single));
end;

function Pad2(a: T2DDoubleArray; M,N: Int32; init:Double=0): T2DDoubleArray; callconv
var W,H,i,j:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H+M*2, W+N*2);
  if init <> 0 then
    for i:=0 to H+M*2-1 do
      for j:=0 to W+N*2-1 do
        Result[i,j] := init;

  for i:=0 to H-1 do
    Move(a[i,0], Result[i+M,N], W*SizeOf(Double));
end;
