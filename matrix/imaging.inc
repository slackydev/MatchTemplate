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
procedure SplitRGB(Image: T2DIntArray; out R,G,B: T2DRealArray);
var
  W,H,x,y:Int32;
begin
  Size(Image, W,H);
  SetLength(R, H,W);
  SetLength(G, H,W);
  SetLength(B, H,W);
  for y:=0 to H-1 do
    for x:=0 to W-1 do begin
      R[y,x] := (Image[y,x]{shr 00}and $FF);
      G[y,x] := (Image[y,x] shr 08 and $FF);
      B[y,x] := (Image[y,x] shr 16 and $FF);
    end;
end;

function CompareImageAt(Image, Templ: T2DIntArray; Pt: TPoint; Tol: Int32): Single;
var
  x,y,w,h,SAD:Int32;
  c1,c2:TRGB32;
begin
  Size(Templ, W,H);
  if (W = 0) or (H = 0) then Exit(1);

  SAD := 0;
  for y:=0 to h-1 do
    for x:=0 to w-1 do
    begin
      c1 := TRGB32(Image[y+Pt.y,x+Pt.x]);
      c2 := TRGB32(Templ[y,x]);
      if (Abs(c1.R-c2.R) < Tol) and
         (Abs(c1.G-c2.G) < Tol) and
         (Abs(c1.B-c2.B) < Tol) then
        Inc(SAD);
    end;
  Result := SAD / (W*H);
end;

function DownscaleImage(Image: T2DIntArray; Scale: Int32): T2DIntArray;
var
  x,y,ys,w,h,nw,nh:Int32;
  invArea: Double;
  
  function BlendArea(x1,y1,x2,y2: Int32): Int32; inline;
  var
    x,y:Int32;
    R:Int32=0; G:Int32=0; B:Int32=0;
  begin
    for y:=y1 to y2 do
      for x:=x1 to x2 do
      begin
        R += Image[y,x] and $FF;
        G += Image[y,x] shr 08 and $FF;
        B += Image[y,x] shr 16 and $FF;
      end;
    R := Round(R * invArea);
    G := Round(G * invArea);
    B := Round(B * invArea);
    Result := R or G shl 8 or B shl 16;
  end;

begin
  Size(Image, w,h);
  if (w = 0) or (h = 0) or (Scale = 0) then 
    Exit;
  
  nw := w div Scale;
  nh := h div Scale;
  invArea := Double(1.0) / Sqr(Scale);
  SetLength(Result, nH, nW);
  for y:=0 to nh-1 do
  begin
    ys := y*scale;
    for x:=0 to nw-1 do
      Result[y,x] := BlendArea(x*scale, ys, x*scale+scale-1, y*scale+scale-1);
  end;
end;


function GetRotatedArea(W, H: Int32; Angle: Single): TBox;
  function Rotate(p:TPoint; angle:Single; mx,my:Int32): TPoint; inline;
  begin
    Result.X := Round(mx + Cos(angle) * (p.x - mx) - Sin(angle) * (p.y - my));
    Result.Y := Round(my + Sin(angle) * (p.x - mx) + Cos(angle) * (p.y - my));
  end;
var
  A,B,C,D: TPoint;
begin
  A := Rotate(Point(0,h), angle, W div 2, H div 2);
  B := Rotate(Point(w,h), angle, W div 2, H div 2);
  C := Rotate(Point(w,0), angle, W div 2, H div 2);
  D := Rotate(Point(0,0), angle, W div 2, H div 2);
  
  Result.X1 := Min(A.x, Min(B.x, Min(C.x, D.x)));
  Result.Y1 := Min(A.y, Min(B.y, Min(C.y, D.y)));
  Result.X2 := Max(A.x, Max(B.x, Max(C.x, D.x)));
  Result.Y2 := Max(A.y, Max(B.y, Max(C.y, D.y)));
end;


function Interpolate(Image: T2DIntArray; x1,y1,x2,y2: Int32; dx,dy: Single): Int32; inline;
var
  R,G,B: Int32;
  p0,p1,p2,p3: TRGB32;
  tR,tG,tB, bR, bG, bB: Single;
begin
  p0 := TRGB32(Image[y1, x1]);
  p1 := TRGB32(Image[y1, x2]);
  p2 := TRGB32(Image[y2, x1]);
  p3 := TRGB32(Image[y2, x2]);

  tR := (1 - dx) * p0.R + dx * p1.R;
  tG := (1 - dx) * p0.G + dx * p1.G;
  tB := (1 - dx) * p0.B + dx * p1.B;
  bR := (1 - dx) * p2.R + dx * p3.R;
  bG := (1 - dx) * p2.G + dx * p3.G;
  bB := (1 - dx) * p2.B + dx * p3.B;

  R := Round((1 - dy) * tR + dy * bR);
  G := Round((1 - dy) * tG + dy * bG);
  B := Round((1 - dy) * tB + dy * bB);

  if (R < 0) then R := 0 else if (R > 255) then R := 255;
  if (G < 0) then G := 0 else if (G > 255) then G := 255;
  if (B < 0) then B := 0 else if (B > 255) then B := 255;

  Result := B or G shl 8 or R shl 16;
end;


(*
 Rotates the bitmap using bilinear interpolation
*)
function __RotateBI(Image: T2DIntArray; Angle: Single; Expand: Boolean): T2DIntArray;
var
  i,j,mx,my,W,H,nW,nH,fX,fY,cX,cY: Int32;
  rX,rY,cosa,sina: Single;
  B: TBox;
begin
  Size(Image, W,H);
  
  if Expand then
  begin
    B := GetRotatedArea(W,H, Angle);
    nW := B.X2 - B.X1 + 1;
    nH := B.Y2 - B.Y1 + 1;
  end else
  begin
    B := Box(0,0,0,0);
    nW := W;
    nH := H;
  end;
  
  mX := nW div 2;
  mY := nH div 2;
  SetLength(Result,nH,nW);
  cosa := Cos(Angle);
  sina := Sin(Angle);
  for i := 0 to nH-1 do begin
    for j := 0 to nW-1 do begin
      rx := mx + cosa * (j - mx) - sina * (i - my);
      ry := my + sina * (j - mx) + cosa * (i - my);
      
      fX := Trunc(rX)+ B.x1;
      fY := Trunc(rY)+ B.y1;
      cX := Ceil(rX) + B.x1;
      cY := Ceil(rY) + B.y1;
      
      if (fX >= 0) and (cX >= 0) and (fX < W) and (cX < W) and
         (fY >= 0) and (cY >= 0) and (fY < H) and (cY < H) then
        Result[i,j] := Interpolate(Image, fx,fy,cx,cy,  rX - (fX - B.x1), rY - (fY - B.y1));
    end;
  end;
end;

(*
 Rotates the bitmap using nearest neighbor
*)
function __RotateNN(Image: T2DIntArray; Angle: Single; Expand: Boolean): T2DIntArray;
var
  nW,nH,W,H,x,y,mx,my,j,i: Int32;
  B: TBox;
  cosa,sina: Single;
begin
  Size(Image, W,H);
  mx := W div 2;
  my := H div 2;
  
  if Expand then
  begin
    B := GetRotatedArea(W,H, Angle);
    nW := B.X2 - B.X1 + 1;
    nH := B.Y2 - B.Y1 + 1;
  end else
  begin
    B := Box(0,0,0,0);
    nW := W;
    nH := H;
  end;
  
  
  SetLength(Result, nH,nW);
  cosa := Cos(angle);
  sina := Sin(angle);
  
  for i:=0 to nH-1 do
    for j:=0 to nW-1 do
    begin
      x := Round(mx + cosa * (B.x1+j - mx) - sina * (B.y1+i - my));
      y := Round(my + sina * (B.x1+j - mx) + cosa * (B.y1+i - my));
      if (x >= 0) and (x < W) and (y >= 0) and (y < H) then
        Result[i,j] := Image[y,x];
    end;
end;


function RotateImage(Image: T2DIntArray; Angle: Single; Expand: Boolean; Smooth: Boolean=True): T2DIntArray;
var W,H: Int32;
begin
  Size(Image, W,H);
  if (W=0) or (H=0) then Exit;
  
  if Smooth then
    Result := __RotateBI(Image, Angle, Expand)
  else
    Result := __RotateNN(Image, Angle, Expand);
end;
