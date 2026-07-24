## =====================================================================
## paper3_diagnostics_refreshed.R
## Runs the full Paper 3 diagnostics on the REFRESHED AQUASTAT extract
## (Africanfreshwater_latest.csv, 1999-2017, with Resources + CODE).
##
## Non-spatial: correlation, KMO, Bartlett, parallel analysis, varimax,
##   eigenvalues, silhouette, k=2 cluster means + ANOVA.
## Spatial: refreshed Factor-1 scores joined to the shapefile geometry by
##   CODE, then global Moran's I + local Moran (LISA).
##
## Output -> paper/paper3_diagnostics_refreshed_output.txt  (I read this).
## Run:  source(".../paper/paper3_diagnostics_refreshed.R")
## =====================================================================

OUT <- "C:/Users/PC/Documents/GitHub/qrfactor/paper"
for (p in c("psych","spdep","sf","cluster"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(qrfactor)

sink(file.path(OUT, "paper3_diagnostics_refreshed_output.txt"), split = TRUE)
cat("Paper 3 diagnostics — REFRESHED AQUASTAT data\n")

dat <- read.csv(file.path(OUT, "Africanfreshwater_latest.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
cat("rows:", nrow(dat), " | year range:",
    paste(range(dat$year, na.rm = TRUE), collapse = "-"), "\n")
cat("columns:", paste(names(dat), collapse = ", "), "\n")

## ---- inspect ALL candidate variables (watch total-vs-sum collinearity) ----
allv <- c("Domestic","Industry","Agricultur","Resources","withdrawal","perCapitaW")
allv <- allv[allv %in% names(dat)]
cat("\n== correlation of ALL candidate variables ==\n")
print(round(cor(dat[allv], use = "pairwise.complete.obs"), 3))
cat("\n(check: is `withdrawal` ~ Domestic+Industry+Agricultur? if r>0.95 drop it)\n")
sums <- rowSums(dat[c("Domestic","Industry","Agricultur")], na.rm = TRUE)
if ("withdrawal" %in% names(dat))
  cat("cor(withdrawal, sum of 3 sectors) =",
      round(cor(dat$withdrawal, sums, use = "complete.obs"), 3), "\n")

## ---- analysis variable set (edit if collinearity dictates) ----------
vars <- c("Domestic","Industry","Agricultur","Resources","withdrawal")
X    <- dat[vars]
cc   <- complete.cases(X)
Xc   <- X[cc, ]
CODEc<- dat$CODE[cc]
cat("\ncomplete cases used:", nrow(Xc), "of", nrow(X), "\n")

## ---- M2 factorability + rotation ------------------------------------
cat("\n== correlation (analysis set) ==\n"); print(round(cor(Xc), 3))
cat("\n== KMO ==\n");     print(try(psych::KMO(Xc)))
cat("\n== Bartlett ==\n");print(try(psych::cortest.bartlett(cor(Xc), n = nrow(Xc))))
cat("\n== parallel analysis ==\n")
pa <- try(psych::fa.parallel(Xc, fa = "pc", n.iter = 100, plot = FALSE), silent = TRUE)
if (!inherits(pa,"try-error")) cat("suggested components:", pa$ncomp, "\n")
cat("\n== varimax 2-factor loadings + communalities ==\n")
pc2 <- try(psych::principal(Xc, nfactors = 2, rotate = "varimax"), silent = TRUE)
if (!inherits(pc2,"try-error")) {
  print(pc2$loadings); cat("\ncommunalities:\n"); print(round(pc2$communality,3))
}
m <- qrfactor(Xc)
cat("\n== qrfactor eigenvalues / variance% / cum% ==\n")
print(round(m$eigen.value,4)); print(round(m$variance,2)); print(round(m$cumvariance,2))

## ---- M7 clusters ----------------------------------------------------
cat("\n== silhouette k=2..6 ==\n")
set.seed(1); sc <- scale(Xc); d <- dist(sc)
for (k in 2:6) {
  km <- kmeans(sc, k, nstart = 25)
  cat("k=",k," sil=",round(mean(cluster::silhouette(km$cluster,d)[,3]),3),"\n")
}
set.seed(1); km2 <- kmeans(sc, 2, nstart = 25)
cat("\n== k=2 cluster mean profiles ==\n")
print(round(aggregate(Xc, list(cluster = km2$cluster), mean), 2))
cat("\n== one-way ANOVA per variable (cluster separation) ==\n")
for (v in vars) {
  a <- summary(aov(Xc[[v]] ~ factor(km2$cluster)))[[1]]
  cat(sprintf("%-12s F=%.2f  p=%s\n", v, a[["F value"]][1],
              format.pval(a[["Pr(>F)"]][1], digits = 3)))
}

## ---- M4 spatial autocorrelation of refreshed Factor 1 ---------------
cat("\n== Moran's I / LISA on refreshed Factor 1 ==\n")
spatial <- try({
  gd <- sf::st_read(system.file("external", package = "qrfactor"),
                    layer = "Africanfreshwater", quiet = TRUE)
  gd$CODE <- trimws(as.character(gd$CODE))
  F1 <- data.frame(CODE = trimws(CODEc), F1 = as.numeric(m$q.loading[,1]))
  gd$F1 <- F1$F1[match(gd$CODE, F1$CODE)]
  nb <- spdep::poly2nb(gd, queen = TRUE)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  mt <- spdep::moran.test(gd$F1, lw, zero.policy = TRUE, na.action = na.omit)
  cat(sprintf("Global Moran's I (Factor 1): I=%.3f  E[I]=%.3f  p=%s\n",
              mt$estimate[[1]], mt$estimate[[2]], format.pval(mt$p.value, digits=3)))
  li <- spdep::localmoran(gd$F1, lw, zero.policy = TRUE, na.action = na.omit)
  tab <- data.frame(COUNTRY = gd$COUNTRY, Ii = round(li[,1],3),
                    p = round(li[,ncol(li)],4))
  cat("\n-- LISA, most significant 12 --\n")
  print(head(tab[order(tab$p),], 12), row.names = FALSE)
}, silent = TRUE)
if (inherits(spatial,"try-error"))
  cat("spatial step failed:", conditionMessage(attr(spatial,"condition")), "\n")

cat("\n\nDONE\n"); sink()
cat("\nWrote", file.path(OUT,"paper3_diagnostics_refreshed_output.txt"), "\n")
