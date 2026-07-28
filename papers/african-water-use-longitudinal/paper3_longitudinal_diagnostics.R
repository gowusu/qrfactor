## =====================================================================
## paper3_longitudinal_diagnostics.R
## Full diagnostics for the LONGITUDINAL Paper 3: Q/R-mode factor analysis
## of African water use at three decadal snapshots (2000, 2010, 2021),
## plus the change/stability metrics that make it a comparison study.
##
## Analysis variable set (per snapshot): Domestic, Industry, Agricultur,
## Resources.  Total withdrawal is reported for context but excluded from
## the ordination (it is the exact sum of the three sectors).
##
## Output -> paper/paper3_longitudinal_output.txt   and
##           paper/longitudinal_scores.csv  (per-country F1/F2 each year).
## =====================================================================

## Set the working directory to the folder holding Africanfreshwater_longitudinal.csv.
OUT <- "."
setwd(OUT)
for (p in c("psych","spdep","sf","cluster"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressWarnings(suppressMessages({library(qrfactor); library(sf); library(spdep)}))

dat <- read.csv("Africanfreshwater_longitudinal.csv", stringsAsFactors = FALSE,
                check.names = FALSE)
years <- c(2000, 2010, 2021)
avar  <- c("Domestic","Industry","Agricultur","Resources")

## bundled geometry (queen-contiguity weights), built once
gd <- sf::st_read(system.file("external", package = "qrfactor"),
                  layer = "Africanfreshwater", quiet = TRUE)
gd$CODE <- trimws(as.character(gd$CODE))
gd <- gd[gd$CODE %in% dat$CODE, ]                 # keep the 49 panel countries
nb <- spdep::poly2nb(gd, queen = TRUE)
lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

sink("paper3_longitudinal_output.txt", split = TRUE)
cat("=====================================================================\n")
cat("Paper 3 LONGITUDINAL diagnostics — AQUASTAT 2000 / 2010 / 2021\n")
cat("n =", nrow(dat), "African countries; analysis vars:",
    paste(avar, collapse=", "), "\n")
cat("=====================================================================\n")

congr <- function(a, b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))   # Tucker congruence
store <- list()

for (ty in years) {
  cat("\n\n########################### YEAR", ty, "###########################\n")
  V  <- paste0(avar, "_", ty)
  X  <- dat[V]; names(X) <- avar
  Tw <- dat[[paste0("withdrawal_", ty)]]

  cat("\n-- continental totals (10^9 m3/yr) --\n")
  cat(sprintf("Total withdrawal SUM = %.1f | Domestic %.1f | Industry %.1f | Agricultur %.1f\n",
              sum(Tw), sum(X$Domestic), sum(X$Industry), sum(X$Agricultur)))
  cat(sprintf("sector shares: Dom %.1f%% Ind %.1f%% Agr %.1f%%\n",
              100*sum(X$Domestic)/sum(Tw), 100*sum(X$Industry)/sum(Tw),
              100*sum(X$Agricultur)/sum(Tw)))

  cat("\n-- correlation matrix --\n"); print(round(cor(X), 3))
  cat("decoupling  cor(Resources, use): Dom",
      round(cor(X$Resources,X$Domestic),3), " Ind", round(cor(X$Resources,X$Industry),3),
      " Agr", round(cor(X$Resources,X$Agricultur),3), "\n")

  cat("\n-- KMO --\n");      print(psych::KMO(X)$MSA); print(round(psych::KMO(X)$MSAi,3))
  bt <- psych::cortest.bartlett(cor(X), n = nrow(X))
  cat("\n-- Bartlett -- chisq=", round(bt$chisq,1), " df=", bt$df,
      " p=", format.pval(bt$p.value, digits=3), "\n")
  pa <- try(psych::fa.parallel(X, fa="pc", n.iter=100, plot=FALSE), silent=TRUE)
  if(!inherits(pa,"try-error")) cat("parallel analysis suggests", pa$ncomp, "component(s)\n")

  pc2 <- psych::principal(X, nfactors = 2, rotate = "varimax")
  cat("\n-- varimax 2-factor loadings --\n"); print(round(unclass(pc2$loadings),3))
  cat("communalities:\n"); print(round(pc2$communality,3))

  m <- qrfactor(X)
  cat("\n-- eigenvalues / variance% / cum% --\n")
  ev <- as.numeric(m$eigen.value); vv <- as.numeric(m$variance); cv <- as.numeric(m$cumvariance)
  print(round(rbind(eigenvalue=ev, variance=vv, cumulative=cv), 3))

  ## per-country factor scores, sign-oriented (use+ with Agricultur, resource+ with Resources)
  F1 <- as.numeric(m$q.loading[,1]); F2 <- as.numeric(m$q.loading[,2])
  if (cor(F1, X$Agricultur) < 0) F1 <- -F1
  if (cor(F2, X$Resources)  < 0) F2 <- -F2
  dat[[paste0("F1_",ty)]] <- F1; dat[[paste0("F2_",ty)]] <- F2

  ## spatial tests
  gd$F1 <- F1[match(gd$CODE, dat$CODE)]; gd$F2 <- F2[match(gd$CODE, dat$CODE)]
  for (fac in c("F1","F2")) {
    mt <- spdep::moran.test(gd[[fac]], lw, zero.policy=TRUE, na.action=na.omit)
    cat(sprintf("Moran's I (%s, %s): I=%.3f E[I]=%.3f p=%s\n", fac, ty,
                mt$estimate[[1]], mt$estimate[[2]], format.pval(mt$p.value,digits=3)))
  }
  li <- spdep::localmoran(gd$F2, lw, zero.policy=TRUE, na.action=na.omit)
  lt <- data.frame(COUNTRY=gd$COUNTRY, Ii=round(li[,1],3), p=round(li[,ncol(li)],4))
  cat("LISA on F2 (resource) — most significant 6:\n")
  print(head(lt[order(lt$p),],6), row.names=FALSE)

  ## clustering diagnostics (is it a gradient or discrete groups?)
  set.seed(1); sc <- scale(X); d <- dist(sc)
  sil <- sapply(2:6, function(k){ km<-kmeans(sc,k,nstart=25)
                mean(cluster::silhouette(km$cluster,d)[,3])})
  cat("silhouette k=2..6:", paste(round(sil,3),collapse=" "), "\n")

  ## top abstractors this year
  rk <- data.frame(COUNTRY=dat$COUNTRY, Total=Tw, F1=F1)
  cat("top 6 abstractors (by total withdrawal):\n")
  print(head(rk[order(-rk$Total),c("COUNTRY","Total")],6), row.names=FALSE)

  store[[as.character(ty)]] <- list(load=unclass(pc2$loadings), F1=F1, F2=F2,
                                    ev=ev, vv=vv, Tw=Tw, X=X,
                                    moranF1=spdep::moran.test(gd$F1,lw,zero.policy=TRUE)$estimate[[1]],
                                    moranF1p=spdep::moran.test(gd$F1,lw,zero.policy=TRUE)$p.value,
                                    moranF2=spdep::moran.test(gd$F2,lw,zero.policy=TRUE)$estimate[[1]],
                                    moranF2p=spdep::moran.test(gd$F2,lw,zero.policy=TRUE)$p.value)
}

## ============================ CHANGE 2000 -> 2021 ==========================
cat("\n\n===================== LONGITUDINAL CHANGE =====================\n")
s0 <- store[["2000"]]; s1 <- store[["2010"]]; s2 <- store[["2021"]]

cat("\n-- continental total withdrawal growth --\n")
tot <- sapply(store, function(s) sum(s$Tw))
cat(sprintf("2000=%.1f  2010=%.1f  2021=%.1f  (x%.2f over the period, +%.0f%%)\n",
            tot[1],tot[2],tot[3], tot[3]/tot[1], 100*(tot[3]/tot[1]-1)))
for (v in c("Domestic","Industry","Agricultur")) {
  g <- sapply(store, function(s) sum(s$X[[v]]))
  cat(sprintf("  %-11s 2000=%.1f 2010=%.1f 2021=%.1f  (x%.2f)\n", v,g[1],g[2],g[3],g[3]/g[1]))
}

cat("\n-- structural stability: Tucker congruence of varimax loadings --\n")
for (fac in 1:2) {
  cat(sprintf("  Factor %d congruence: 2000~2010=%.3f  2010~2021=%.3f  2000~2021=%.3f\n", fac,
      congr(s0$load[,fac],s1$load[,fac]), congr(s1$load[,fac],s2$load[,fac]),
      congr(s0$load[,fac],s2$load[,fac])))
}

cat("\n-- variance explained by Factor 1 (use) over time --\n")
cat(sprintf("  2000=%.1f%%  2010=%.1f%%  2021=%.1f%%\n", s0$vv[1],s1$vv[1],s2$vv[1]))

cat("\n-- rank persistence: correlation of F1 (use) scores across years --\n")
cat(sprintf("  r(2000,2010)=%.3f  r(2010,2021)=%.3f  r(2000,2021)=%.3f (Spearman)\n",
    cor(s0$F1,s1$F1,method="spearman"), cor(s1$F1,s2$F1,method="spearman"),
    cor(s0$F1,s2$F1,method="spearman")))

cat("\n-- spatial signature over time (Moran's I) --\n")
cat(sprintf("  Factor 1 (use):     2000 I=%.3f (p=%s)  2010 I=%.3f (p=%s)  2021 I=%.3f (p=%s)\n",
    s0$moranF1,format.pval(s0$moranF1p,digits=2), s1$moranF1,format.pval(s1$moranF1p,digits=2),
    s2$moranF1,format.pval(s2$moranF1p,digits=2)))
cat(sprintf("  Factor 2 (resource):2000 I=%.3f (p=%s)  2010 I=%.3f (p=%s)  2021 I=%.3f (p=%s)\n",
    s0$moranF2,format.pval(s0$moranF2p,digits=2), s1$moranF2,format.pval(s1$moranF2p,digits=2),
    s2$moranF2,format.pval(s2$moranF2p,digits=2)))

cat("\n-- fastest-growing abstractors 2000->2021 (absolute rise in total withdrawal) --\n")
gr <- data.frame(COUNTRY=dat$COUNTRY,
                 t2000=s0$Tw, t2021=s2$Tw,
                 abs_rise=s2$Tw - s0$Tw,
                 pct = ifelse(s0$Tw>0, 100*(s2$Tw/s0$Tw-1), NA))
cat("by absolute rise:\n");  print(head(gr[order(-gr$abs_rise),],8), row.names=FALSE)
cat("\nby % growth (base >= 0.3 km3/yr):\n")
grp <- gr[s0$Tw>=0.3,]; print(head(grp[order(-grp$pct),c("COUNTRY","t2000","t2021","pct")],8), row.names=FALSE)

cat("\n-- decoupling persistence: max |cor(Resources, any use var)| each year --\n")
for (ty in years) {
  X <- store[[as.character(ty)]]$X
  cc <- c(cor(X$Resources,X$Domestic),cor(X$Resources,X$Industry),cor(X$Resources,X$Agricultur))
  cat(sprintf("  %d: cor range [%.3f, %.3f]\n", ty, min(cc), max(cc)))
}

## largest abstractors table (for the paper) with endowment
cat("\n-- profile table: top abstractors 2021, with 2000 value and resources --\n")
prof <- data.frame(COUNTRY=dat$COUNTRY,
                   Total_2000=round(s0$Tw,2), Total_2021=round(s2$Tw,2),
                   Resources=dat$Resources_2021)
print(head(prof[order(-prof$Total_2021),],8), row.names=FALSE)
cat("\nmost water-rich (by resources), with their 2021 use:\n")
print(head(prof[order(-prof$Resources),c("COUNTRY","Total_2021","Resources")],6), row.names=FALSE)

## save per-country scores + key values for figures
write.csv(dat, "longitudinal_scores.csv", row.names = FALSE)
cat("\n\nwrote longitudinal_scores.csv\nDONE\n")
sink()
