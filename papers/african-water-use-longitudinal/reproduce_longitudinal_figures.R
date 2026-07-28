###############################################################################
#  reproduce_longitudinal_figures.R
#  Base-R figures for the LONGITUDINAL Paper 3 (the non-qrfactor-native ones):
#    fig1_growth            — sectoral growth 2000/2010/2021
#    fig4_moran_over_time   — Moran's I of the two factors across the three dates
#    fig5_cluster           — k-means (k=2) cluster plot: an outlier-led split
#    fig6_change_dumbbell   — total withdrawal 2000 vs 2021, who is rising
#  (Figures 2 and 3 are qrfactor's OWN plots — see reproduce_qrfactor_native_plots.R)
#  Reads longitudinal_scores.csv (from paper3_longitudinal_diagnostics.R).
###############################################################################

## Set the working directory to the folder holding longitudinal_scores.csv.
OUT <- "."
FIG <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE)
setwd(OUT)
png_ <- function(name, w=1100, h=850, res=130) png(file.path(FIG,name), width=w, height=h, res=res)

dat   <- read.csv("longitudinal_scores.csv", stringsAsFactors=FALSE, check.names=FALSE)
years <- c(2000, 2010, 2021)
accent<- "#c1462f"
short <- c("Democratic Republic of the Congo"="DR Congo","Central African Republic"="CAR",
           "United Republic of Tanzania"="Tanzania")
disp  <- ifelse(dat$COUNTRY %in% names(short), short[dat$COUNTRY], dat$COUNTRY)

## ============================================================================
## Figure 1 — continental withdrawal grows and shifts composition, 2000->2021
## ============================================================================
sect <- sapply(years, function(ty) c(
  Agricultural = sum(dat[[paste0("Agricultur_",ty)]]),
  Municipal    = sum(dat[[paste0("Domestic_",ty)]]),
  Industrial   = sum(dat[[paste0("Industry_",ty)]])))
colnames(sect) <- years
png_("fig1_growth.png", w=1150, h=780)
par(mar=c(4.5,4.8,4,1))
bp <- barplot(sect, col=c("#4a7c59","#e9c46a","#8d99ae"), border="white",
              names.arg=years, ylim=c(0,230),
              ylab=expression("Water withdrawal ("*10^9*" m"^3*" yr"^-1*")"),
              main="African water withdrawal grew 23% over 2000–2021,\nled by municipal and industrial use")
tot <- colSums(sect)
text(bp, tot+8, sprintf("%.0f", tot), font=2, cex=0.95)
legend("topleft", rev(rownames(sect)), fill=rev(c("#4a7c59","#e9c46a","#8d99ae")),
       bty="n", cex=0.9, border="white")
text(bp[1], 120, "Agriculture\n81% → 77%", cex=0.72, col="white", font=2)
dev.off()

## ============================================================================
## Figure 4 — the asymmetry is stable: Moran's I of both factors over time
## ============================================================================
moranF1 <- c(-0.001,-0.015,-0.010); moranF2 <- c(0.202,0.208,0.199)
png_("fig4_moran_over_time.png", w=1050, h=760)
par(mar=c(4.4,4.6,4,1))
plot(years, moranF2, type="o", pch=19, col="#023047", lwd=2, ylim=c(-0.08,0.28),
     xaxt="n", xlab="Year", ylab="Global Moran's I of factor scores",
     main="A stable spatial asymmetry:\nsupply always clusters, use never does")
axis(1, at=years)
lines(years, moranF1, type="o", pch=17, col=accent, lwd=2)
abline(h=-0.023, lty=3, col="grey60")
text(2010, -0.045, "E[I] under no autocorrelation", cex=0.62, col="grey45")
text(2021, 0.199, "  renewable availability\n  (p < 0.01 throughout)", pos=2, cex=0.66, col="#023047", font=3)
text(2021, -0.010, "  scale of use\n  (n.s. throughout)", pos=2, cex=0.66, col=accent, font=3)
legend("topleft", c("Factor 2 — availability","Factor 1 — use"),
       pch=c(19,17), col=c("#023047",accent), lwd=2, bty="n", cex=0.82)
dev.off()

## ============================================================================
## Figure 5 — k-means (k=2) cluster plot on the interpretable factor axes:
##            the split isolates Egypt (an outlier-led gradient, not two groups)
## ============================================================================
avar <- c("Domestic","Industry","Agricultur","Resources")
X21  <- dat[paste0(avar,"_2021")]; names(X21) <- avar
set.seed(1); km <- kmeans(scale(X21), 2, nstart=25)
## label clusters so cluster "A" is the big group, "B" the Egypt-led extreme
big  <- which.max(table(km$cluster)); cl <- ifelse(km$cluster==big, 1L, 2L)
F1 <- dat$F1_2021; F2 <- dat$F2_2021
pal <- c("#1f7a8c","#c1462f")
png_("fig5_cluster.png", w=1080, h=900)
par(mar=c(4.6,4.6,4,1))
xlim <- range(F1) + c(-0.05, 0.18)*diff(range(F1))
plot(F1, F2, type="n", xlim=xlim,
     xlab="Factor 1  —  scale of water use", ylab="Factor 2  —  renewable availability",
     main=paste0("k-means (k = 2) isolates the extreme abstractor:\n",
                 "an outlier-led gradient, not two balanced groups"))
abline(h=0, v=0, col="grey85", lty=3)
## convex hulls per cluster (shaded), singleton clusters skipped
for (g in 1:2) { idx <- which(cl==g)
  if (length(idx) >= 3) { h <- idx[chull(F1[idx], F2[idx])]
    polygon(F1[h], F2[h], col=adjustcolor(pal[g], 0.12), border=adjustcolor(pal[g],0.5)) } }
points(F1, F2, pch=19, col=adjustcolor(pal[cl],0.9), cex=1.05)
lab <- dat$COUNTRY %in% c("Egypt","South Africa","Nigeria","Madagascar","Algeria",
                          "Morocco","Ethiopia","Democratic Republic of the Congo","Congo","Cameroon")
lpos <- ifelse(dat$COUNTRY=="Egypt", 2, 4)
text(F1[lab], F2[lab], disp[lab], pos=lpos[lab], cex=0.6, col="#0a3d52")
szB <- sum(cl==2); szA <- sum(cl==1)
legend("topright", bty="n", cex=0.82, pch=19, col=pal,
       legend=c(sprintf("main group (n = %d)", szA),
                sprintf("high-use extreme (n = %d)", szB)))
text(mean(range(F1)), max(F2)*0.9,
     "silhouette peaks at k = 2 (0.86) —\nthe split peels off the outlier, it does not\ndivide the mass into balanced classes",
     cex=0.62, col="grey35", font=3)
dev.off()
cat("cluster sizes (1=main, 2=extreme):", szA, szB, " | extreme members:",
    paste(dat$COUNTRY[cl==2], collapse=", "), "\n")

## ============================================================================
## Figure 6 — who is rising: dumbbell of total withdrawal 2000 -> 2021
## ============================================================================
tw0 <- dat$withdrawal_2000; tw2 <- dat$withdrawal_2021
o   <- order(tw2); top <- tail(o, 18)
png_("fig6_change_dumbbell.png", w=1050, h=1050)
par(mar=c(4.5,7.5,4,1))
plot(NA, xlim=c(0, max(tw2[top])*1.02), ylim=c(1,length(top)),
     yaxt="n", xlab=expression("Total water withdrawal ("*10^9*" m"^3*" yr"^-1*")"),
     ylab="", main="Who is driving the growth, 2000 → 2021")
grew <- tw2[top] >= tw0[top]
segments(tw0[top], seq_along(top), tw2[top], seq_along(top),
         col=ifelse(grew,"#adb5bd","#e07a5f"), lwd=2.4)
points(tw0[top], seq_along(top), pch=19, col="#8ecae6", cex=1.1)
points(tw2[top], seq_along(top), pch=19, col="#023047", cex=1.1)
axis(2, at=seq_along(top), labels=disp[top], las=1, cex.axis=0.62)
fast <- c("Kenya","Ethiopia","South Africa","Algeria","Niger","Senegal")
for (ct in fast) { i <- which(disp[top]==ct | dat$COUNTRY[top]==ct)
  if(length(i)) text(tw2[top][i], i, sprintf("  +%.0f%%",100*(tw2[top][i]/tw0[top][i]-1)),
                     pos=4, cex=0.6, col="#023047", font=3, xpd=NA) }
legend("bottomright", c("2000","2021","decline"), pch=c(19,19,NA), lwd=c(NA,NA,2.4),
       col=c("#8ecae6","#023047","#e07a5f"), bty="n", cex=0.85)
dev.off()

## ============================================================================
## Figure 7 — trajectory in factor space 2000 -> 2010 -> 2021: stable geometry,
##   changing magnitude. All snapshots scored on ONE pooled varimax solution
##   (a common set of axes), so movement reflects real change, not per-year
##   rescaling. log used for legibility (raw volumes are heavily right-skewed).
## ============================================================================
suppressWarnings(suppressMessages(library(psych)))
avarT <- c("Domestic","Industry","Agricultur","Resources")
long <- do.call(rbind, lapply(years, function(ty){
  d <- dat[paste0(avarT,"_",ty)]; names(d) <- avarT
  data.frame(CODE=dat$CODE, COUNTRY=dat$COUNTRY, year=ty, d, stringsAsFactors=FALSE)
}))
pc <- principal(log1p(long[avarT]), nfactors=2, rotate="varimax", scores=TRUE)
L  <- unclass(pc$loadings); sc <- as.data.frame(pc$scores); names(sc) <- c("A","B")
useax <- if (abs(L["Resources","RC1"]) > abs(L["Resources","RC2"])) "B" else "A"
long$F1 <- sc[[useax]]; long$F2 <- sc[[setdiff(c("A","B"), useax)]]
if (cor(long$F1, long$Agricultur) < 0) long$F1 <- -long$F1
if (cor(long$F2, long$Resources)  < 0) long$F2 <- -long$F2

w1 <- reshape(long[c("CODE","COUNTRY","year","F1")], idvar=c("CODE","COUNTRY"),
              timevar="year", direction="wide")
w2 <- reshape(long[c("CODE","year","F2")], idvar="CODE", timevar="year", direction="wide")
w  <- merge(w1, w2, by="CODE")
dispw <- ifelse(w$COUNTRY %in% names(short), short[w$COUNTRY], w$COUNTRY)

png_("fig7_trajectory.png", w=1180, h=1000)
par(mar=c(4.6,4.8,4.4,1))
plot(range(long$F1)+c(-0.1,0.5), range(long$F2)+c(-0.15,0.15), type="n",
     xlab="Factor 1  —  scale of water use (common pooled axis)",
     ylab="Factor 2  —  renewable availability (time-invariant)",
     main="Stable geometry, changing magnitude (2000 → 2021):\ncountries slide along use while availability holds them in place")
abline(h=0, v=0, col="grey88", lty=3)
grew <- w$F1.2021 >= w$F1.2000
arrows(w$F1.2000, w$F2.2000, w$F1.2021, w$F2.2021, length=0.06,
       col=ifelse(grew,"#1f7a8c","#e07a5f"), lwd=1.8)
points(w$F1.2010, w$F2.2010, pch=19, col="grey65", cex=0.5)
points(w$F1.2000, w$F2.2000, pch=1,  col="grey45", cex=0.7)
lab <- w$COUNTRY %in% c("Egypt","South Africa","Nigeria","Madagascar","Algeria","Ethiopia",
                        "Kenya","Morocco","Democratic Republic of the Congo","Congo","Cameroon")
lpz <- ifelse(dispw=="Egypt", 2, 4)
text(w$F1.2021[lab], w$F2.2021[lab], dispw[lab], pos=lpz[lab], cex=0.58, col="#0a3d52")
legend("topright", bty="n", cex=0.78,
       legend=c("2000 (start)","2010","2021 (arrowhead)","grew","declined"),
       pch=c(1,19,NA,NA,NA), lwd=c(NA,NA,NA,1.8,1.8),
       col=c("grey45","grey65",NA,"#1f7a8c","#e07a5f"))
dev.off()

message("Wrote fig1_growth, fig4_moran_over_time, fig5_cluster, fig6_change_dumbbell, fig7_trajectory")
