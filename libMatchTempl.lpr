library libMatchTempl;
{==============================================================================]
  Copyright (c) 2018, Jarl `slacky` Holta
  Project: libfft
  License: GNU Lesser GPL (http://www.gnu.org/licenses/lgpl.html)
[==============================================================================}
{$I header.inc}

uses
  SysUtils, Math, core, matchTempl;

function MatchTemplate_Wrap(constref Img, Sub: T2DIntArray; TMFormula: ETMFormula): T2DRealArray; callconv
begin
  Result := MatchTemplate(Img, Sub, TMFormula);
end;

{$I SimbaPlugin.inc}


end.
