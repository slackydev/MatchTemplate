{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 For more info see: Copyright.txt
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
{$DEFINE GET_VALUES_BODY :=
  L := High(Indices);
  Size(a, W,H);
  if H < 0 then Exit;
  
  SetLength(Result, L+1);
  c := 0;
  for i:=0 to L do
    if (Indices[i].x >= 0) and (Indices[i].y >= 0) and 
       (Indices[i].x <  W) and (Indices[i].y <  H) then
    begin
      Result[c] := a[Indices[i].y,Indices[i].x];
      Inc(c);
    end;
  SetLength(Result, c);
}

function GetValues(a: T2DIntArray; Indices: TPointArray): TIntArray;
var i,W,H,c,L:Int32;
begin 
  GET_VALUES_BODY
end;

function GetValues(a: T2DSingleArray; Indices: TPointArray): TSingleArray;
var i,W,H,c,L:Int32;
begin 
  GET_VALUES_BODY
end;

function GetValues(a: T2DDoubleArray; Indices: TPointArray): TDoubleArray;
var i,W,H,c,L:Int32;
begin 
  GET_VALUES_BODY
end;
