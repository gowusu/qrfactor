## Resubmission of an archived package

`qrfactor` was archived on CRAN in 2018. This 1.5 release restores it to
working order for current R.

The archival cause was an S3 registration defect: the package exported its
methods with `exportPattern(".")` but never registered them with
`S3method()`, so `qrfactor()` failed with
`no applicable method for 'qrfactor' applied to an object of class
"data.frame"`. This is fixed, and the methods are now explicitly
registered.

In addition, the retired spatial packages the original depended on
(`rgdal`, `maptools`, `mgraph`) have been removed and replaced with `sf`
and `sp`, and a crash under R >= 4.2 (a condition of length > 1) in the
map-plotting code has been fixed.

## Test environments

* Local: Windows 10 x64, R 4.6.0 (2026-04-24 ucrt)

## R CMD check results

0 errors | 0 warnings | 0 notes

On CRAN's incoming checks a NOTE is expected for a new submission /
previously archived package.

## Downstream dependencies

There are no reverse dependencies (the package was archived).
