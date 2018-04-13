{=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=]
 Copyright (c) 2018, Jarl K. <Slacky> Holta || http://github.com/slackydev
 All rights reserved.
 For more info see: Copyright.txt
[=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=}
procedure Size(a: T2DIntArray; out W,H: Int32);
begin 
  H := Length(a); 
  if H > 0 then  W := Length(a[0])
  else           W := 0;
end;

procedure Size(a: T2DI64Array; out W,H: Int32);
begin
  H := Length(a);
  if H > 0 then  W := Length(a[0])
  else           W := 0;
end;

procedure Size(a: T2DSingleArray; out W,H: Int32);
begin 
  H := Length(a); 
  if H > 0 then  W := Length(a[0])
  else           W := 0;
end;

procedure Size(a: T2DDoubleArray; out W,H: Int32);
begin 
  H := Length(a); 
  if H > 0 then  W := Length(a[0])
  else           W := 0;
end;

procedure Size(a: T2DComplexArray; out W,H: Int32);
begin
  H := Length(a);
  if H > 0 then  W := Length(a[0])
  else           W := 0;
end;
