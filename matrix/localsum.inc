{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 For more info see: Copyright.txt
 Computes the local sum around the template size
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
function LocalSum(a: T2DIntArray; th,tw: Int32): T2DIntArray; callconv
var i,j,w,h:Int32;
begin
  a := Pad2(a, th,tw);
  h := Length(a);
  w := Length(a[0]);
  for j:=0 to w-1 do
  begin
    for i:=1 to h-1    do a[i,j] += a[i-1,j];
    for i:=0 to h-th-1 do a[i,j] := a[i+th,j] - a[i,j];
  end;

  for i:=0 to h-1 do
  begin
    for j:=1 to w-1    do a[i,j] += a[i,j-1];
    for j:=0 to w-tw-1 do a[i,j] := a[i,j+tw] - a[i,j];
  end;
  SetLength(a, h-th-1, w-tw-1);
  Result := a;
end;

function LocalSum(a: T2DSingleArray; th,tw: Int32): T2DSingleArray; callconv
var i,j,w,h:Int32;
begin
  a := Pad2(a, th,tw);
  h := Length(a);
  w := Length(a[0]);
  for j:=0 to w-1 do
  begin
    for i:=1 to h-1    do a[i,j] += a[i-1,j];
    for i:=0 to h-th-1 do a[i,j] := a[i+th,j] - a[i,j];
  end;

  for i:=0 to h-1 do
  begin
    for j:=1 to w-1    do a[i,j] += a[i,j-1];
    for j:=0 to w-tw-1 do a[i,j] := a[i,j+tw] - a[i,j];
  end;
  SetLength(a, h-th-1, w-tw-1);
  Result := a;
end;

function LocalSum(a: T2DDoubleArray; th,tw: Int32): T2DDoubleArray; callconv
var i,j,w,h:Int32;
begin
  a := Pad2(a, th,tw);
  h := Length(a);
  w := Length(a[0]);
  for j:=0 to w-1 do
  begin
    for i:=1 to h-1    do a[i,j] += a[i-1,j];
    for i:=0 to h-th-1 do a[i,j] := a[i+th,j] - a[i,j];
  end;

  for i:=0 to h-1 do
  begin
    for j:=1 to w-1    do a[i,j] += a[i,j-1];
    for j:=0 to w-tw-1 do a[i,j] := a[i,j+tw] - a[i,j];
  end;
  SetLength(a, h-th-1, w-tw-1);
  Result := a;
end;


