# qrfactor

<!-- badges: start -->
<!-- badges: end -->

**Simultaneous Q-mode and R-mode factor analysis for spatial data.**

`qrfactor` performs Q-mode and R-mode factor analysis *simultaneously* on
spatial and non-spatial data. A single function, `qrfactor()`, carries out

- principal component analysis (PCA),
- R-mode factor analysis,
- Q-mode factor analysis,
- simultaneous R- and Q-mode factor analysis,
- principal coordinate analysis (PCoA), and
- multidimensional scaling (MDS).

Loadings and scores are returned from the fitted object, and the `plot()`
method provides annotated biplots for combinations of eigenvectors,
loadings and scores. Input may be supplied as an ESRI shapefile, a
delimited text file or a data frame.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("gowusu/qrfactor")
```

For the spatial (shapefile) workflow you will also need the `sf` and `sp`
packages:

```r
install.packages(c("sf", "sp"))
```

## Quick start

Non-spatial (data frame) example:

```r
library(qrfactor)
data(UScereal, package = "MASS")
v <- c("calories", "protein", "sodium", "carbo", "sugars", "potassium")
m <- qrfactor(UScereal[v], scale = "n")
plot(m, plot = "r", type = "loading")
```

Spatial (shapefile) example, using the bundled African freshwater data:

```r
p <- system.file("external", package = "qrfactor")
m <- qrfactor(p, layer = "Africanfreshwater",
              var = c("Domestic", "Industry", "Agricultur"))
plot(m, plot = "map")
```

## Status

This is the 1.5 restoration of the package (originally archived from CRAN
in 2018). The archival cause — an unregistered S3 method — has been fixed,
and the retired `rgdal` / `maptools` / `mgraph` dependencies have been
replaced with an `sf`/`sp` compatibility layer.

## Author

George Owusu (<gowusu@gmail.com>)

## License

GPL-2.
