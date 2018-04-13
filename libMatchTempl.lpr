library libMatchTempl;
{==============================================================================]
  Copyright (c) 2018, Jarl `slacky` Holta
  Project: libfft
  License: GNU Lesser GPL (http://www.gnu.org/licenses/lgpl.html)
[==============================================================================}
{$I header.inc}

uses
  SysUtils, Math, core, matchTempl;

function MatchTemplate_Wrap(constref Image, Templ: T2DIntArray; TMFormula: ETMFormula): T2DRealArray; callconv
begin
  Result := MatchTemplate(Image, Templ, TMFormula);
end;

{$I SimbaPlugin.inc}


end.
