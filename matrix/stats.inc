{==============================================================================]
  Copyright © 2018, Jarl Krister Holta
  
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
[==============================================================================}
// ----------------------------------------------------------------------------
// Mean
function Mean(a: T2DIntArray): Double;
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

function Mean(a: T2DSingleArray): Double;
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

function Mean(a: T2DDoubleArray): Double;
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

procedure MeanStdev(a: T2DIntArray; out Mean, Stdev: Double);
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

procedure MeanStdev(a: T2DSingleArray; out Mean, Stdev: Double);
var
  W,H,i,j:Int32;
  sum, square: Double;
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

procedure MeanStdev(a: T2DDoubleArray; out Mean, Stdev: Double);
var
  W,H,i,j:Int32;
  sum, square: Double;
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
// MinMax

procedure MinMax(a: T2DIntArray; out vMin,vMax: Int32);
var
  W,H,i,j:Int32;
begin
  Size(a, W,H);
  vMin := 0; 
  vMax := 0;
  if H = 0 then Exit;
  
  vMin := a[0,0]; vMax := a[0,0];
  for i:=0 to h-1 do
    for j:=0 to w-1 do
    begin
      if a[i,j] > vMax then vMax := a[i,j];
      if a[i,j] < vMin then vMin := a[i,j];
    end;
end; 

procedure MinMax(a: T2DSingleArray; out vMin,vMax: Single);
var
  W,H,i,j:Int32;
begin
  Size(a, W,H);
  vMin := 0; 
  vMax := 0;
  if H = 0 then Exit;
  
  vMin := a[0,0]; vMax := a[0,0];
  for i:=0 to h-1 do
    for j:=0 to w-1 do
    begin
      if a[i,j] > vMax then vMax := a[i,j];
      if a[i,j] < vMin then vMin := a[i,j];
    end;
end; 

procedure MinMax(a: T2DDoubleArray; out vMin,vMax: Double);
var
  W,H,i,j:Int32;
begin
  Size(a, W,H);
  vMin := 0; 
  vMax := 0;
  if H = 0 then Exit;
  
  vMin := a[0,0]; vMax := a[0,0];
  for i:=0 to h-1 do
    for j:=0 to w-1 do
    begin
      if a[i,j] > vMax then vMax := a[i,j];
      if a[i,j] < vMin then vMin := a[i,j];
    end;
end; 


// ----------------------------------------------------------------------------
// Various (TReal)


function ArgMax(a: T2DRealArray): TPoint;
var X,Y,W,H:Int32;
begin
  Result := Point(0,0);
  Size(a, W,H);

  for Y:=0 to H-1 do
    for X:=0 to W-1 do
      if a[Y,X] > a[Result.y,Result.x] then
      begin
        Result.x := x;
        Result.y := y;
      end;
end;

function ArgMin(a: T2DRealArray): TPoint;
var X,Y,W,H:Int32;
begin
  Result := Point(0,0);
  Size(a, W,H);

  for Y:=0 to H-1 do
    for X:=0 to W-1 do
      if a[Y,X] < a[Result.y,Result.x] then
      begin
        Result.x := x;
        Result.y := y;
      end;
end;

function NormMinMax(a: T2DRealArray; Alpha, Beta: TReal): T2DRealArray;
var
  Lo,Hi,oldRange,newRange: TReal;
  X,Y,W,H: Int32;
begin
  Size(a, W,H);
  SetLength(Result, H,W);
  MinMax(a, Lo,Hi);

  oldRange := Hi-Lo;
  newRange := Beta-Alpha;
  if oldRange = 0 then Exit;
  for Y:=0 to H-1 do
    for X:=0 to W-1 do
      Result[Y,X] := (a[Y,X] - lo) / oldRange * newRange + Alpha;
end;
