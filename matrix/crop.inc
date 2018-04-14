{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 For more info see: Copyright.txt
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
function Crop(a: T2DIntArray; B: TBox): T2DIntArray;
var
  y,w,h:Int32;
begin
  Size(a, W,H);
  if (H = 0) then Exit;
  
  B.x1 := Min(B.x1, W-1);
  B.y1 := Min(B.y1, H-1);
  B.x2 := Min(B.x2, W-1);
  B.y2 := Min(B.y2, H-1);
  
  SetLength(Result, B.y2-B.y1+1, B.x2-B.x1+1);
  for Y:=B.y1 to B.y2 do
    Move(a[y][B.x1], Result[y-B.y1,0], (B.x2-B.x1+1)*SizeOf(a[0,0]));
end;

function Crop(a: T2DSingleArray; B: TBox): T2DSingleArray;
var y,w,h:Int32;
begin
  Size(a, W,H);
  if (H = 0) then Exit;
  
  B.x1 := Min(B.x1, W-1);
  B.y1 := Min(B.y1, H-1);
  B.x2 := Min(B.x2, W-1);
  B.y2 := Min(B.y2, H-1);
  
  SetLength(Result, B.y2-B.y1+1, B.x2-B.x1+1);
  for Y:=B.y1 to B.y2 do
    Move(a[y][B.x1], Result[y-B.y1,0], (B.x2-B.x1+1)*SizeOf(a[0,0]));
end;

function Crop(a: T2DDoubleArray; B: TBox): T2DDoubleArray;
var y,w,h:Int32;
begin
  Size(a, W,H);
  if (H = 0) then Exit;
  
  B.x1 := Min(B.x1, W-1);
  B.y1 := Min(B.y1, H-1);
  B.x2 := Min(B.x2, W-1);
  B.y2 := Min(B.y2, H-1);
  
  SetLength(Result, B.y2-B.y1+1, B.x2-B.x1+1);
  for Y:=B.y1 to B.y2 do
    Move(a[y][B.x1], Result[y-B.y1,0], (B.x2-B.x1+1)*SizeOf(a[0,0]));
end;
