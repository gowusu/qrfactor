## =====================================================================
## qrfactor manual — reproduction script
##
## Runs EVERY qrfactor function/option and captures the results so the
## manual can be written from real output (not guesses).
##
##   * all console output  -> manual/repro_output.txt   (via sink)
##   * all figures         -> manual/figures/*.pdf      (multi-page)
##
## Run once from RGui:
##   source("C:/Users/PC/Documents/GitHub/qrfactor/manual/reproduce_manual.R")
##
## Needs: qrfactor (installed), sf, sp, MASS, cluster, mvoutlier, pvclust.
## =====================================================================

OUT  <- "C:/Users/PC/Documents/GitHub/qrfactor/manual"
FIG  <- file.path(OUT, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

library(qrfactor)

## Small helpers -------------------------------------------------------
sec  <- function(title) cat("\n\n========== ", title, " ==========\n")
show <- function(label, value) { cat("\n---", label, "---\n"); print(value) }
resetpar <- function() { try(par(mfrow = c(1,1)), silent = TRUE) }
## save a plotting expression to a PNG (named), trapping errors so one
## failure never halts the run. multi=TRUE -> one PNG per page (%02d),
## for calls that draw many pages (maps, paginated diagnostics).
figpng <- function(name, expr, w = 900, h = 600, multi = FALSE) {
  file <- file.path(FIG, if (multi) paste0(name, "_%02d.png") else paste0(name, ".png"))
  cat("\n[FIGURE] ", name, ": ")
  png(file, width = w, height = h, res = 110)
  res <- try(force(expr), silent = TRUE)
  try(dev.off(), silent = TRUE)
  if (inherits(res, "try-error"))
    cat("FAIL -> ", conditionMessage(attr(res, "condition")), "\n")
  else cat("ok\n")
  resetpar()
  invisible(res)
}

## =====================================================================
## Begin capture
## =====================================================================
sink(file.path(OUT, "repro_output.txt"), split = TRUE)
cat("qrfactor manual reproduction — run at capture time\n")
cat("qrfactor version:", as.character(utils::packageVersion("qrfactor")), "\n")
cat("R version:", R.version.string, "\n")

## ---------------------------------------------------------------------
sec("0. Data sources")
csv    <- system.file("external", "Africanfreshwater.csv", package = "qrfactor")
source <- system.file("external", package = "qrfactor")
layer  <- "Africanfreshwater"
var    <- c("Domestic","Industry","Agricultur","Resources","withdrawal","perCapitaW")
cat("csv path:   ", csv, "\n")
cat("shape dir:  ", source, "\n")
raw <- read.csv(csv)
show("head(raw)", head(raw))
show("dim(raw)",  dim(raw))
show("names(raw)", names(raw))
show("variables used", var)

## For the simplest intro (matches the package tests): UScereal
data(UScereal, package = "MASS")
uvars <- c("calories","protein","sodium","carbo","sugars","potassium")

## =====================================================================
## PART 1 — NON-SPATIAL (data-frame / CSV)
## =====================================================================

## ---- 1. Quick start (UScereal), the object -------------------------
sec("1. Quick start + the fitted object (UScereal)")
mUS <- qrfactor(UScereal[uvars], scale = "n")
show("class(mUS)", class(mUS))
show("names(mUS)  [the whole 'workbench']", names(mUS))
show("str of key parts", utils::str(mUS[c("eigen.value","variance","cumvariance",
                                          "r.loading","q.loading")]))

## ---- 2. Build the running model from the CSV -----------------------
sec("2. qrfactor() on the freshwater CSV")
mod <- qrfactor(csv, var = var)
show("mod$call", mod$call)
show("dim(mod$data)", dim(mod$data))
show("head(mod$data)", head(mod$data))

## ---- 3. print() and summary() --------------------------------------
sec("3. print(mod)")
print(mod)
sec("3b. summary(mod)")
summary(mod)

## ---- 4. PCA --------------------------------------------------------
sec("4. PCA")
show("eigen.value", mod$eigen.value)
show("standard deviations sqrt(eigen.value)", sqrt(mod$eigen.value))
show("variance (% explained)", mod$variance)
show("cumulative variance", mod$cumvariance)
show("pca loadings (mod$pca)", mod$pca)
show("pca scores head", head(mod$pca.scores))

## ---- 5. R-mode -----------------------------------------------------
sec("5. R-mode factor analysis")
show("r.loading", mod$r.loading)

## ---- 6. Q-mode -----------------------------------------------------
sec("6. Q-mode factor analysis")
show("q.loading (head)", head(mod$q.loading))
show("q.scores (head)",  head(mod$q.scores))

## ---- 7. Simultaneous Q + R ----------------------------------------
sec("7. Simultaneous Q and R mode")
show("combined loadings (head)", head(mod$loadings))
show("combined scores (head)",   head(mod$scores))
show("eigen.vector", mod$eigen.vector)
show("correlation matrix", round(mod$correlation, 3))

## ---- 8. scale = variants (show the effect on eigenvalues) ----------
sec("8. Effect of scale= on the decomposition")
for (sc in c("sd","normal","n","centre","data")) {
  m <- try(qrfactor(csv, var = var, scale = sc), silent = TRUE)
  if (inherits(m, "try-error")) { cat("\nscale=", sc, "FAILED\n"); next }
  cat("\nscale =", sc, "\n")
  cat("  eigen.value :", paste(round(m$eigen.value,3), collapse=" "), "\n")
  cat("  cumvariance :", paste(round(m$cumvariance,1), collapse=" "), "\n")
}

## ---- 9. transform argument (t=) ------------------------------------
sec("9. Variable transform (t=) and scale log/sqrt")
mlog <- try(qrfactor(csv, var = var, scale = "log"), silent = TRUE)
show("scale='log' eigen.value", if (inherits(mlog,"try-error")) mlog else mlog$eigen.value)
mtr  <- try(qrfactor(csv, var = var, t = c("withdrawal","Resources")), silent = TRUE)
show("t=c('withdrawal','Resources') eigen.value",
     if (inherits(mtr,"try-error")) mtr else mtr$eigen.value)

## ---- 10. All the non-spatial plots (named PNGs) --------------------
sec("10. Non-spatial plots -> figures/*.png")
figpng("fig_pca_all",     plot(mod, type = "PCA", plot = "all"))
figpng("fig_pca_loadings",plot(mod, type = "pca"))   # lowercase: true PCA-loadings plot
figpng("fig_ca_view",     plot(mod, type = "ca"))    # correspondence-labelled view
figpng("fig_rmode",       plot(mod, type = "loadings", plot = "r"))
figpng("fig_qmode",       plot(mod, type = "loadings", plot = "q"))
figpng("fig_qmode_scores",plot(mod, type = "scores",  plot = "q"))
figpng("fig_simultaneous",plot(mod, type = "loadings"))
figpng("fig_scores",      plot(mod, type = "scores"))
figpng("fig_plot_all",    plot(mod, type = "loadings", plot = "all"))
figpng("fig_scores_qr",   plot(mod, type = "scores", plot = "qr"))
figpng("fig_pcoa",        plot(qrfactor(csv,var=var), plot="all", type="coord", factors=c(1,3)))
figpng("fig_mds",         plot(qrfactor(csv,var=var,scale="n"), plot="r", type="mds", factors=c(1,2)))
figpng("fig_labelled",    plot(mod, rowname = "COUNTRY"))
figpng("fig_cex_values",  plot(mod, cex = c("withdrawal"), values = c("Domestic"), pch = 23))
figpng("fig_factors13",   plot(mod, factors = c(1,3)))
figpng("fig_xlim_ylim",   plot(mod, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5), main = "Freshwater Q/R biplot"))
## documented-argument stress cases (each exercises a scalar-condition path)
figpng("fig_values_true", plot(mod, cex = c("withdrawal"), values = TRUE, pch = 23))
figpng("fig_abline_num",  plot(mod, abline = c(-0.5, 0.5)))
figpng("fig_xlim_qr",     plot(mod, xlim = "qr"))
figpng("fig_rd_style",    plot(mod, rowname = "COUNTRY", cex = c("withdrawal"),
                                legend = "topleft", values = TRUE, pch = 23))
mtr_ok <- try(qrfactor(csv, var = var, t = c("withdrawal","Resources")), silent = TRUE)
if (!inherits(mtr_ok, "try-error"))
  figpng("fig_transformed", plot(mtr_ok, main = "sqrt-transformed withdrawal & Resources"))

## =====================================================================
## PART 2 — SPATIAL (shapefile)
## =====================================================================
sec("11. Spatial model from the shapefile")
smod <- qrfactor(source, layer = layer, var = var)
show("class(smod)", class(smod))
show("class(smod$gisdata)", class(smod$gisdata))
show("feature count nrow(smod$gisdata)", nrow(smod$gisdata))
## the write-back: score/cluster/index columns added to the attribute table
gd <- as.data.frame(smod$gisdata)
show("attribute-table columns after fit", names(gd))
show("a known row (Algeria) with its scores/clusters",
     gd[gd$COUNTRY == "Algeria",
        intersect(c("COUNTRY","Domestic","Factor1","Factor2","index","index2",
                    "means","cluster","cluster1","cluster2"), names(gd))])

## ---- 12. Maps -> figures/fig_map_*.png (one per spplot page) -------
sec("12. Spatial maps -> figures/fig_map_*.png")
figpng("fig_map", plot(smod, plot = "map"), multi = TRUE)

## ---- 13. Diagnostics -> figures/*.png ------------------------------
sec("13. Diagnostic plots -> figures/*.png")
## sink() the huge anova/kruskal console dumps to a side file so they do
## not swamp the main log; we only need to know they run.
figpng("fig_diagnose", plot(smod, type = "diagnose"), multi = TRUE)
figpng("fig_cluster",  plot(smod, plot = "cluster"),  multi = TRUE)
figpng("fig_compare",  plot(smod, plot = "compare"),  multi = TRUE)
sink_anova <- file(file.path(OUT, "repro_anova_dump.txt"), open = "wt")
sink(sink_anova)
figpng("fig_anova",         plot(smod, plot = "anova"), multi = TRUE)
figpng("fig_nonparametric", plot(smod, plot = "nonparametric"), multi = TRUE)
sink()                       # close the anova side-sink (main sink resumes)
close(sink_anova)

## ---- 14. In-memory spatial object + CSV<->shape join ---------------
sec("14. In-memory spatial object, and the m=/f= join")
gisobj <- smod$gisdata                      # a Spatial*DataFrame
mmem <- try(qrfactor(gisobj, layer = "gisobject", var = var), silent = TRUE)
show("qrfactor on in-memory object: class", if (inherits(mmem,"try-error")) mmem else class(mmem))
mjoin <- try(qrfactor(source, layer = layer, var = var, m = "COUNTRY", f = csv), silent = TRUE)
show("join via m='COUNTRY', f=csv: attribute columns",
     if (inherits(mjoin,"try-error")) mjoin else names(as.data.frame(mjoin$gisdata)))

## ---- 15. Accessor catalogue ----------------------------------------
sec("15. Every element the object exposes")
for (nm in names(mod)) {
  v <- mod[[nm]]
  cat(sprintf("  %-18s %-12s dim/len: %s\n", nm, class(v)[1],
              if (!is.null(dim(v))) paste(dim(v), collapse="x") else length(v)))
}

sec("16. sessionInfo()")
print(sessionInfo())

cat("\n\nDONE. Figures in:", FIG, "\n")
sink()

## Console breadcrumb (outside the sink)
cat("\nReproduction complete.\n",
    "Text  ->", file.path(OUT, "repro_output.txt"), "\n",
    "Figs  ->", FIG, "(fig_nonspatial.pdf, fig_maps.pdf, fig_diagnostics.pdf)\n")
