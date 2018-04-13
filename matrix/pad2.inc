{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 For more info see: Copyright.txt
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
function Pad2(a: T2DIntArray; M,N: Int32): T2DIntArray;
var W,H,i:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H+M*2, W+N*2);
  for i:=0 to H-1 do
    Move(a[i,0], Result[i+M,N], W*SizeOf(Int32));
end;

function Pad2(a: T2DI64Array; M,N: Int32): T2DI64Array;
var W,H,i:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H+M*2, W+N*2);
  for i:=0 to H-1 do
    Move(a[i,0], Result[i+M,N], W*SizeOf(Int32));
end;

function Pad2(a: T2DSingleArray; M,N: Int32): T2DSingleArray;
var W,H,i:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H+M*2, W+N*2);
  for i:=0 to H-1 do
    Move(a[i,0], Result[i+M,N], W*SizeOf(Single));
end;

function Pad2(a: T2DDoubleArray; M,N: Int32): T2DDoubleArray;
var W,H,i:Int32;
begin
  Size(a, w,h);
  SetLength(Result, H+M*2, W+N*2);
  for i:=0 to H-1 do
    Move(a[i,0], Result[i+M,N], W*SizeOf(Double));
end;
