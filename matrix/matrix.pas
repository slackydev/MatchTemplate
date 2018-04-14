unit matrix;
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
{$I header.inc}

interface

uses
  SysUtils, core;

procedure Size(a: T2DIntArray; out W,H: Int32);
procedure Size(a: T2DI64Array; out W,H: Int32);
procedure Size(a: T2DSingleArray; out W,H: Int32);
procedure Size(a: T2DDoubleArray; out W,H: Int32);
procedure Size(a: T2DComplexArray; out W,H: Int32);

procedure SplitRGB(Image: T2DIntArray; out R,G,B: T2DRealArray);
function CompareImageAt(Image, Templ: T2DIntArray; Pt: TPoint; Tol: Int32): Single;
function DownscaleImage(Image: T2DIntArray; Scale: Int32): T2DIntArray;
function RotateImage(Image: T2DIntArray; Angle: Single; Expand: Boolean; Smooth: Boolean=True): T2DIntArray;

function Crop(a: T2DIntArray; B: TBox): T2DIntArray;
function Crop(a: T2DSingleArray; B: TBox): T2DSingleArray;
function Crop(a: T2DDoubleArray; B: TBox): T2DDoubleArray;

function GetValues(a: T2DIntArray; Indices: TPointArray): TIntArray;
function GetValues(a: T2DSingleArray; Indices: TPointArray): TSingleArray;
function GetValues(a: T2DDoubleArray; Indices: TPointArray): TDoubleArray;

function Sums(a: T2DIntArray; out Square: T2DI64Array): T2DIntArray;
function Sums(a: T2DSingleArray; out Square: T2DDoubleArray): T2DDoubleArray;
function Sums(a: T2DDoubleArray; out Square: T2DDoubleArray): T2DDoubleArray;

function SumsPd(a: T2DIntArray; out Square: T2DI64Array): T2DIntArray;
function SumsPd(a: T2DSingleArray; out Square: T2DDoubleArray): T2DDoubleArray;
function SumsPd(a: T2DDoubleArray; out Square: T2DDoubleArray): T2DDoubleArray;

function Rot90(a: T2DIntArray): T2DIntArray;
function Rot90(a: T2DSingleArray): T2DSingleArray;
function Rot90(a: T2DDoubleArray): T2DDoubleArray;
function Rot90(a: T2DComplexArray): T2DComplexArray;

function Mean(a: T2DIntArray): Double;
function Mean(a: T2DSingleArray): Double;
function Mean(a: T2DDoubleArray): Double;

procedure MeanStdev(a: T2DIntArray; out Mean, Stdev: Double);
procedure MeanStdev(a: T2DSingleArray; out Mean, Stdev: Double);
procedure MeanStdev(a: T2DDoubleArray; out Mean, Stdev: Double);

procedure MinMax(a: T2DIntArray; out vMin,vMax: Int32);
procedure MinMax(a: T2DSingleArray; out vMin,vMax: Single);
procedure MinMax(a: T2DDoubleArray; out vMin,vMax: Double);

function ArgMin(a: T2DRealArray): TPoint;
function ArgMax(a: T2DRealArray): TPoint;
function NormMinMax(a: T2DRealArray; Alpha, Beta: TReal): T2DRealArray;

implementation

uses
  Math;

{$I size.inc}
{$I imaging.inc}
{$I crop.inc}
{$I getValues.inc}
{$I sums.inc}
{$I rot90.inc}
{$I stats.inc}

end.
