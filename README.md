Pure FPC implementation of FFT based template matching.

Template matching is a technique for finding areas of an image that match (are similar) to a template image (smaller image).

Exported:
```
type ETMFormula = (TM_CCORR, TM_CCORR_NORMED, TM_CCOEFF, TM_CCOEFF_NORMED, TM_SQDIFF, TM_SQDIFF_NORMED);
function MatchTemplate(constref img, sub: T2DIntArray; Formula: ETMFormula=TM_CCOEFF_NORMED): T2DRealArray;
```

