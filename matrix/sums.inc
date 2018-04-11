{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 For more info see: Copyright.txt
 Computes the local sum around the template size
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
function Sums(a: T2DIntArray; out Square: T2DIntArray): T2DIntArray;
var
  x,y,W,H: Int32;
  sum,sqsum: Int64;
begin
  H := Length(a);
  W := Length(a[0]);
  SetLength(Result, H,W);
  SetLength(Square, H,W);

  Result[0,0] := a[0,0];
  Square[0,0] := Sqr(a[0,0]);
  for y:=1 to H-1 do
  begin
    Result[y,0] := Result[y-1,0] + a[y,0];
    Square[y,0] := Square[y-1,0] + Sqr(a[y,0]);
  end;

  for x:=1 to W-1 do
  begin
    Result[0,x] := Result[0,x-1] + a[0,x];
    Square[0,x] := Square[0,x-1] + Sqr(a[0,x]);
  end;

  for y:=1 to H-1 do
  begin
    sum   := a[y,0];
    sqsum := Sqr(sum);
    for x:=1 to W-1 do
    begin
      sum += a[y,x];
      Result[y,x] := Result[y-1,x] + sum;
      sqsum += Sqr(a[y,x]);
      Square[y,x] := Square[y-1,x] + sqsum;
    end;
  end;
end;


function Sums(a: T2DSingleArray; out Square: T2DDoubleArray): T2DDoubleArray;
var
  x,y,W,H: Int32;
  sum,sqsum: Double;
begin
  H := Length(a);
  W := Length(a[0]);
  SetLength(Result, H,W);
  SetLength(Square, H,W);

  Result[0,0] := a[0,0];
  Square[0,0] := Sqr(a[0,0]);
  for y:=1 to H-1 do
  begin
    Result[y,0] := Result[y-1,0] + a[y,0];
    Square[y,0] := Square[y-1,0] + Sqr(a[y,0]);
  end;

  for x:=1 to W-1 do
  begin
    Result[0,x] := Result[0,x-1] + a[0,x];
    Square[0,x] := Square[0,x-1] + Sqr(a[0,x]);
  end;

  for y:=1 to H-1 do
  begin
    sum   := a[y,0];
    sqsum := Sqr(sum);
    for x:=1 to W-1 do
    begin
      sum += a[y,x];
      Result[y,x] := Result[y-1,x] + sum;
      sqsum += Sqr(a[y,x]);
      Square[y,x] := Square[y-1,x] + sqsum;
    end;
  end;
end;

function Sums(a: T2DDoubleArray; out Square: T2DDoubleArray): T2DDoubleArray;
var
  x,y,W,H: Int32;
  sum,sqsum: Double;
begin
  H := Length(a);
  W := Length(a[0]);
  SetLength(Result, H,W);
  SetLength(Square, H,W);

  Result[0,0] := a[0,0];
  Square[0,0] := Sqr(a[0,0]);
  for y:=1 to H-1 do
  begin
    Result[y,0] := Result[y-1,0] + a[y,0];
    Square[y,0] := Square[y-1,0] + Sqr(a[y,0]);
  end;

  for x:=1 to W-1 do
  begin
    Result[0,x] := Result[0,x-1] + a[0,x];
    Square[0,x] := Square[0,x-1] + Sqr(a[0,x]);
  end;

  for y:=1 to H-1 do
  begin
    sum   := a[y,0];
    sqsum := Sqr(sum);
    for x:=1 to W-1 do
    begin
      sum += a[y,x];
      Result[y,x] := Result[y-1,x] + sum;
      sqsum += Sqr(a[y,x]);
      Square[y,x] := Square[y-1,x] + sqsum;
    end;
  end;
end;




function SumsPd(a: T2DIntArray; out Square: T2DIntArray): T2DIntArray;
var
  x,y,W,H: Int32;
  sum,sqsum: Int64;
begin
  H := Length(a);
  W := Length(a[0]);
  SetLength(Result, H+1,W+1);
  SetLength(Square, H+1,W+1);

  Result[1,1] := a[0,0];
  Square[1,1] := Sqr(a[0,0]);
  for y:=2 to H do
  begin
    Result[y,1] := Result[y-1,1] + a[y-1,0];
    Square[y,1] := Square[y-1,1] + Sqr(a[y-1,0]);
  end;

  for x:=2 to W do
  begin
    Result[1,x] := Result[1,x-1] + a[0,x-1];
    Square[1,x] := Square[1,x-1] + Sqr(a[0,x-1]);
  end;

  for y:=2 to H do
  begin
    sum   := a[y-1,0];
    sqsum := Sqr(sum);
    for x:=2 to W do
    begin
      sum += a[y-1,x-1];
      Result[y,x] := Result[y-1,x] + sum;
      sqsum += Sqr(a[y-1,x-1]);
      Square[y,x] := Square[y-1,x] + sqsum;
    end;
  end;
end;

function SumsPd(a: T2DSingleArray; out Square: T2DDoubleArray): T2DDoubleArray;
var
  x,y,W,H: Int32;
  sum,sqsum: Double;
begin
  H := Length(a);
  W := Length(a[0]);
  SetLength(Result, H+1,W+1);
  SetLength(Square, H+1,W+1);

  Result[1,1] := a[0,0];
  Square[1,1] := Sqr(a[0,0]);
  for y:=2 to H do
  begin
    Result[y,1] := Result[y-1,1] + a[y-1,0];
    Square[y,1] := Square[y-1,1] + Sqr(a[y-1,0]);
  end;

  for x:=2 to W do
  begin
    Result[1,x] := Result[1,x-1] + a[0,x-1];
    Square[1,x] := Square[1,x-1] + Sqr(a[0,x-1]);
  end;

  for y:=2 to H do
  begin
    sum   := a[y-1,0];
    sqsum := Sqr(sum);
    for x:=2 to W do
    begin
      sum += a[y-1,x-1];
      Result[y,x] := Result[y-1,x] + sum;
      sqsum += Sqr(a[y-1,x-1]);
      Square[y,x] := Square[y-1,x] + sqsum;
    end;
  end;
end;

function SumsPd(a: T2DDoubleArray; out Square: T2DDoubleArray): T2DDoubleArray;
var
  x,y,W,H: Int32;
  sum,sqsum: Double;
begin
  H := Length(a);
  W := Length(a[0]);
  SetLength(Result, H+1,W+1);
  SetLength(Square, H+1,W+1);

  Result[1,1] := a[0,0];
  Square[1,1] := Sqr(a[0,0]);
  for y:=2 to H do
  begin
    Result[y,1] := Result[y-1,1] + a[y-1,0];
    Square[y,1] := Square[y-1,1] + Sqr(a[y-1,0]);
  end;

  for x:=2 to W do
  begin
    Result[1,x] := Result[1,x-1] + a[0,x-1];
    Square[1,x] := Square[1,x-1] + Sqr(a[0,x-1]);
  end;

  for y:=2 to H do
  begin
    sum   := a[y-1,0];
    sqsum := Sqr(sum);
    for x:=2 to W do
    begin
      sum += a[y-1,x-1];
      Result[y,x] := Result[y-1,x] + sum;
      sqsum += Sqr(a[y-1,x-1]);
      Square[y,x] := Square[y-1,x] + sqsum;
    end;
  end;
end;

