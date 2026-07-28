## =====================================================================
## paper3_robustness.R  — referee-proofing checks for the longitudinal paper
##  (1) log-transform sensitivity of the ordination
##  (2) Egypt-excluded sensitivity (the dominant Factor-1 outlier)
##  (3) k-nearest-neighbour spatial weights vs queen contiguity (Moran's I)
## Output -> paper/paper3_robustness_output.txt
## =====================================================================
## Set the working directory to the folder holding longitudinal_scores.csv.
OUT <- "."
setwd(OUT)
suppressWarnings(suppressMessages({library(qrfactor); library(psych); library(sf); library(spdep)}))

dat  <- read.csv("longitudinal_scores.csv", stringsAsFactors=FALSE, check.names=FALSE)
years<- c(2000,2010,2021)
avar <- c("Domestic","Industry","Agricultur","Resources")
congr<- function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))

varimax_load <- function(X) unclass(principal(X, nfactors=2, rotate="varimax")$loadings)
var_pct      <- function(X){ m<-qrfactor(X); as.numeric(m$variance) }
maxdecouple  <- function(X) max(abs(c(cor(X$Resources,X$Domestic),
                                      cor(X$Resources,X$Industry),
                                      cor(X$Resources,X$Agricultur))))

sink("paper3_robustness_output.txt", split=TRUE)
cat("================ ROBUSTNESS CHECKS (qrfactor Paper 3) ================\n")

## baseline 2021 loadings for congruence comparisons
X21 <- dat[paste0(avar,"_2021")]; names(X21) <- avar
base_load <- varimax_load(X21)

## ---- (1) LOG TRANSFORM ------------------------------------------------
cat("\n### (1) log1p-transform sensitivity ###\n")
for (ty in years) {
  X <- dat[paste0(avar,"_",ty)]; names(X) <- avar
  Xl <- as.data.frame(lapply(X, log1p))
  vp <- var_pct(Xl); L <- varimax_load(Xl)
  cat(sprintf("%d  log: F1=%.1f%% F2=%.1f%%  KMO=%.2f  max|r(Res,use)|=%.3f  congruence(F1 vs raw-2021)=%.3f\n",
      ty, vp[1], vp[2], KMO(Xl)$MSA, maxdecouple(Xl),
      congr(L[,1], base_load[,1])))
}
cat("raw 2021 for reference: F1=68.8% F2=25.0% KMO=0.74 max|r|=0.045\n")
cat("log 2021 correlation matrix (use vars still collinear? resources still decoupled?):\n")
Xl21 <- as.data.frame(lapply(X21, log1p)); print(round(cor(Xl21),3))

## ---- (2) EGYPT EXCLUDED ----------------------------------------------
cat("\n### (2) Egypt-excluded sensitivity (2021) ###\n")
noeg <- dat$COUNTRY != "Egypt"
Xn <- dat[noeg, paste0(avar,"_2021")]; names(Xn) <- avar
vp <- var_pct(Xn); L <- varimax_load(Xn)
cat(sprintf("n=%d (Egypt removed)  F1=%.1f%% F2=%.1f%%  KMO=%.2f  max|r(Res,use)|=%.3f  congruence(F1 vs full)=%.3f\n",
    nrow(Xn), vp[1], vp[2], KMO(Xn)$MSA, maxdecouple(Xn), congr(L[,1], base_load[,1])))
cat("Egypt-excluded correlation matrix:\n"); print(round(cor(Xn),3))
cat("Egypt-excluded varimax loadings:\n"); print(round(L,3))

## ---- (3) k-NEAREST-NEIGHBOUR WEIGHTS vs QUEEN ------------------------
cat("\n### (3) spatial weights robustness: queen vs k-NN (k=4) ###\n")
gd <- st_read(system.file("external", package="qrfactor"), layer="Africanfreshwater", quiet=TRUE)
gd$CODE <- trimws(as.character(gd$CODE)); gd <- gd[gd$CODE %in% dat$CODE,]
ord <- match(gd$CODE, dat$CODE)
ctr <- suppressWarnings(st_coordinates(st_centroid(st_geometry(gd))))
nbQ <- poly2nb(gd, queen=TRUE);          lwQ <- nb2listw(nbQ, style="W", zero.policy=TRUE)
knn <- knearneigh(ctr, k=4);             lwK <- nb2listw(knn2nb(knn), style="W")
n_isl <- sum(card(nbQ)==0)
cat("island states with 0 queen neighbours (get neighbours under k-NN):", n_isl, "\n")
for (ty in years) {
  for (fac in c("F1","F2")) {
    v  <- dat[[paste0(fac,"_",ty)]][ord]
    mQ <- moran.test(v, lwQ, zero.policy=TRUE, na.action=na.omit)
    mK <- moran.test(v, lwK, na.action=na.omit)
    cat(sprintf("%d %s: queen I=%.3f (p=%s) | kNN4 I=%.3f (p=%s)\n", ty, fac,
        mQ$estimate[[1]], format.pval(mQ$p.value,digits=2),
        mK$estimate[[1]], format.pval(mK$p.value,digits=2)))
  }
}
cat("\nDONE\n"); sink()
