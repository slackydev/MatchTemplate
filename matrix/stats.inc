{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 For more info see: Copyright.txt
 Mean, MeanStdev, MinMax, 
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
// ----------------------------------------------------------------------------
// Mean
function Mean(a: T2DIntArray): Double; callconv
var
  W,H,i,j:Int32;
  sum: Double;
begin
  Size(a, W,H);
  sum := 0;
  for i:=0 to h-1 do
    for j:=0 to w-1 do sum += a[i,j];
  Result := sum / (W*H);
end; 

function Mean(a: T2DSingleArray): Double; callconv
var
  W,H,i,j:Int32;
  sum: Double;
begin
  Size(a, W,H);
  sum := 0;
  for i:=0 to h-1 do
    for j:=0 to w-1 do sum += a[i,j];
  Result := sum / (W*H);
end; 

function Mean(a: T2DDoubleArray): Double; callconv
var
  W,H,i,j:Int32;
  sum: Double;
begin
  Size(a, W,H);
  sum := 0;
  for i:=0 to h-1 do
    for j:=0 to w-1 do sum += a[i,j];
  Result := sum / (W*H);
end; 


// ----------------------------------------------------------------------------
// Mean and Stdev at once

procedure MeanStdev(a: T2DIntArray; out Mean, Stdev: Double); callconv
var
  W,H,i,j:Int32;
  sum, square: Extended;
begin
  Size(a, W,H);
  sum := 0;
  for i:=0 to h-1 do
    for j:=0 to w-1 do sum += a[i,j];
  Mean := sum / (W*H);

  square := 0;
  for i:=0 to h-1 do
    for j:=0 to w-1 do square += Sqr(a[i,j] - Mean);
  Stdev := Sqrt(square / (W*H));
end;

procedure MeanStdev(a: T2DSingleArray; out Mean, Stdev: Single); callconv
var
  W,H,i,j:Int32;
  sum, square: Extended;
begin
  Size(a, W,H);
  sum := 0;
  for i:=0 to h-1 do
    for j:=0 to w-1 do sum += a[i,j];
  Mean := sum / (W*H);

  square := 0;
  for i:=0 to h-1 do
    for j:=0 to w-1 do square += Sqr(a[i,j] - Mean);
  Stdev := Sqrt(square / (W*H));
end;

procedure MeanStdev(a: T2DDoubleArray; out Mean, Stdev: Double); callconv
var
  W,H,i,j:Int32;
  sum, square: Extended;
begin
  Size(a, W,H);
  sum := 0;
  for i:=0 to h-1 do
    for j:=0 to w-1 do sum += a[i,j];
  Mean := sum / (W*H);

  square := 0;
  for i:=0 to h-1 do
    for j:=0 to w-1 do square += Sqr(a[i,j] - Mean);
  Stdev := Sqrt(square / (W*H));
end;

// ----------------------------------------------------------------------------
// Mean

procedure MinMax(a: T2DIntArray; out vMin,vMax: Int32); callconv
var
  W,H,i,j:Int32;
begin
  Size(a, W,H);
  vMin := a[0,0]; vMax := a[0,0];
  for i:=0 to h-1 do
    for j:=0 to w-1 do
    begin
      if a[i,j] > vMax then vMax := a[i,j];
      if a[i,j] < vMin then vMin := a[i,j];
    end;
end; 

procedure MinMax(a: T2DSingleArray; out vMin,vMax: Single); callconv
var
  W,H,i,j:Int32;
begin
  Size(a, W,H);
  vMin := a[0,0]; vMax := a[0,0];
  for i:=0 to h-1 do
    for j:=0 to w-1 do
    begin
      if a[i,j] > vMax then vMax := a[i,j];
      if a[i,j] < vMin then vMin := a[i,j];
    end;
end; 

procedure MinMax(a: T2DDoubleArray; out vMin,vMax: Double); callconv
var
  W,H,i,j:Int32;
begin
  Size(a, W,H);
  vMin := a[0,0]; vMax := a[0,0];
  for i:=0 to h-1 do
    for j:=0 to w-1 do
    begin
      if a[i,j] > vMax then vMax := a[i,j];
      if a[i,j] < vMin then vMin := a[i,j];
    end;
end; 
