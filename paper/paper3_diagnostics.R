## =====================================================================
## paper3_diagnostics.R
## Reviewer-requested statistics for Paper 3 (addresses referee M2/M3/M4/M7).
##
## Produces, on the CURRENT bundled data (re-run later on the refreshed
## AQUASTAT extract, unchanged):
##   M3  cleaned variable set (drops perCapitaW, r=0.99 with Agricultur)
##   M2  KMO sampling adequacy, Bartlett sphericity, communalities,
##       VARIMAX-rotated loadings, Horn's parallel analysis (retention)
##   M7  cluster-number diagnostics (silhouette by k), seeded kmeans
##   M4  Global Moran's I + local Moran (LISA) on the factor indices
##
## Output -> paper/paper3_diagnostics_output.txt  (I read this file).
## Run:  source("C:/Users/PC/Documents/GitHub/qrfactor/paper/paper3_diagnostics.R")
## =====================================================================

OUT <- "C:/Users/PC/Documents/GitHub/qrfactor/paper"

## --- dependencies ----------------------------------------------------
need <- c("psych","spdep","sf","cluster","MASS")
for (p in need)
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(qrfactor)

sink(file.path(OUT, "paper3_diagnostics_output.txt"), split = TRUE)
cat("Paper 3 diagnostics — run at capture time\n")
cat("qrfactor", as.character(utils::packageVersion("qrfactor")),
    "| psych", as.character(utils::packageVersion("psych")),
    "| spdep", as.character(utils::packageVersion("spdep")), "\n")

## --- data ------------------------------------------------------------
csv <- system.file("external", "Africanfreshwater.csv", package = "qrfactor")
src <- system.file("external", package = "qrfactor")
raw <- read.csv(csv)

## M3: cleaned variable set — drop perCapitaW (collinear with Agricultur)
vars <- c("Domestic","Industry","Agricultur","Resources","withdrawal")
X <- raw[vars]
cat("\n\n========== M3: cleaned variable set ==========\n")
cat("variables:", paste(vars, collapse = ", "), "\n")
cat("n =", nrow(X), "\n")

## --- M2: is the matrix factorable? -----------------------------------
cat("\n\n========== M2: factorability ==========\n")
R <- cor(X)
cat("\n-- correlation matrix --\n"); print(round(R, 3))
cat("\n-- KMO (overall MSA + per variable) --\n"); print(try(psych::KMO(X)))
cat("\n-- Bartlett's test of sphericity --\n")
print(try(psych::cortest.bartlett(R, n = nrow(X))))

## --- retention: Horn's parallel analysis (answers Kaiser critique) ---
cat("\n\n========== retention: parallel analysis ==========\n")
pa <- try(psych::fa.parallel(X, fa = "pc", n.iter = 100, plot = FALSE), silent = TRUE)
if (!inherits(pa, "try-error")) {
  cat("suggested # components (parallel analysis):", pa$ncomp, "\n")
  cat("component eigenvalues:\n"); print(round(pa$pc.values, 3))
}

## --- M2: varimax-rotated loadings + communalities --------------------
cat("\n\n========== M2: 2-factor PCA, varimax rotated ==========\n")
pc2 <- try(psych::principal(X, nfactors = 2, rotate = "varimax"), silent = TRUE)
if (!inherits(pc2, "try-error")) {
  print(pc2$loadings)
  cat("\ncommunalities:\n"); print(round(pc2$communality, 3))
  cat("\nuniquenesses:\n"); print(round(pc2$uniquenesses, 3))
  cat("\nvariance accounted for:\n"); print(round(pc2$Vaccounted, 3))
}

## qrfactor on the cleaned set (consistency with the package's own basis)
m <- qrfactor(X)
cat("\n-- qrfactor eigenvalues (cleaned set) --\n"); print(round(m$eigen.value, 4))
cat("-- variance % --\n");     print(round(m$variance, 2))
cat("-- cumulative % --\n");   print(round(m$cumvariance, 2))

## --- M7: how many clusters? ------------------------------------------
cat("\n\n========== M7: cluster-number diagnostics ==========\n")
set.seed(1)
sc <- scale(X)
d  <- dist(sc)
for (k in 2:6) {
  km  <- kmeans(sc, k, nstart = 25)
  sil <- cluster::silhouette(km$cluster, d)
  cat("k =", k, "  avg silhouette width =", round(mean(sil[, 3]), 3),
      "  within-SS =", round(km$tot.withinss, 1), "\n")
}

## --- M4: spatial autocorrelation of the factor indices ---------------
cat("\n\n========== M4: Moran's I / LISA ==========\n")
sm <- qrfactor(src, layer = "Africanfreshwater", var = vars)
gd <- sf::st_as_sf(sm$gisdata)                 # Spatial* -> sf for spdep
nb <- spdep::poly2nb(gd, queen = TRUE)         # queen contiguity
lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
cat("-- neighbour structure --\n"); print(nb)

for (col in c("Factor1","Factor2","index","index2")) {
  if (col %in% names(gd)) {
    v  <- as.numeric(gd[[col]])
    mt <- try(spdep::moran.test(v, lw, zero.policy = TRUE, na.action = na.omit),
              silent = TRUE)
    if (!inherits(mt, "try-error")) {
      cat(sprintf("\nGlobal Moran's I  [%s]:  I = %.3f   E[I] = %.3f   p = %s\n",
                  col, mt$estimate[[1]], mt$estimate[[2]],
                  format.pval(mt$p.value, digits = 3)))
    }
  }
}

## Local Moran (LISA) for the Factor 1 index — the strongest-clustered countries
if ("index" %in% names(gd)) {
  lisa <- try(spdep::localmoran(as.numeric(gd$index), lw, zero.policy = TRUE),
              silent = TRUE)
  if (!inherits(lisa, "try-error")) {
    tab <- data.frame(COUNTRY = gd$COUNTRY,
                      Ii = round(lisa[, 1], 3),
                      p  = round(lisa[, ncol(lisa)], 4))
    cat("\n-- Local Moran (LISA) for Factor 1 index, most significant 15 --\n")
    print(head(tab[order(tab$p), ], 15), row.names = FALSE)
  }
}

cat("\n\nDONE\n")
sink()
cat("\nWrote", file.path(OUT, "paper3_diagnostics_output.txt"), "\n")
