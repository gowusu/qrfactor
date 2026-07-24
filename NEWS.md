# qrfactor 1.5

Restoration release. `qrfactor` was archived from CRAN in 2018; this
version repairs it for current R and CRAN.

## Bug fixes

* Registered the S3 methods (`qrfactor.default`, `print.qrfactor`,
  `summary.qrfactor`, `plot.qrfactor`) in `NAMESPACE`. The archived
  version exported them with `exportPattern(".")` but never registered
  them, so `qrfactor()` failed with
  `no applicable method for 'qrfactor' applied to an object of class
  "data.frame"`. This was the cause of the CRAN archival.

* Fixed a crash under R >= 4.2 in the map plots, where a condition of
  length greater than one is now an error rather than a silent
  first-element take.

## Dependency changes

* Replaced retired spatial packages. Shapefiles are now read with
  `sf::st_read()` instead of `rgdal::readOGR()`, and the `maptools` and
  `mgraph` dependencies were removed. A small compatibility layer bridges
  the `sf`/`sp` object differences.

* `Depends:` on contributed packages moved to `Imports:`; base-R
  functions are now explicitly imported via `importFrom()`.

## Documentation

* Rewrote the `DESCRIPTION` (complete-sentence Description, modern
  `Authors@R`).
