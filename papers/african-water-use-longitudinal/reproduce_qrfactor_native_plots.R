###############################################################################
#  reproduce_qrfactor_native_plots.R
#  Regenerates qrfactor's OWN signature graphics — the simultaneous Q/R biplot
#  and the "Factor 1 / Factor 2 Index" choropleth maps — on the longitudinal
#  panel (2000 & 2021), so the paper showcases the package's native output.
#  Requires the spatial model (shapefile joined to the snapshot values).
###############################################################################

## Set the working directory to the folder holding longitudinal_scores.csv.
OUT <- "."
FIG <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE)
suppressWarnings(suppressMessages({library(qrfactor); library(sf); library(sp); library(lattice)}))

dat  <- read.csv(file.path(OUT, "longitudinal_scores.csv"), stringsAsFactors=FALSE, check.names=FALSE)
avar <- c("Domestic","Industry","Agricultur","Resources")

## local copies of qrfactor's internal label helpers (they live inside plot.qrfactor)
sp.label     <- function(x, label, cex, pos) list("sp.text", coordinates(x), label, cex=cex, pos=pos)
ISO.sp.label <- function(x, var, cex, pos) if (is.null(var)) sp.label(x, row.names(x), cex, pos) else sp.label(x, x[[var]], cex, pos)
label        <- function(x, var=NULL, cex=0.6, pos=1) do.call("list", ISO.sp.label(x, var, cex, pos))

## build a spatial qrfactor model for one snapshot (49 panel countries) ---------
build_year <- function(ty) {
  gd <- st_read(system.file("external", package="qrfactor"), layer="Africanfreshwater", quiet=TRUE)
  gd$CODE <- trimws(as.character(gd$CODE))
  gd <- gd[gd$CODE %in% dat$CODE, ]
  mm <- match(gd$CODE, dat$CODE)
  for (v in avar) gd[[v]] <- dat[[paste0(v,"_",ty)]][mm]
  m <- qrfactor(source = as(gd, "Spatial"), layer = "gisobject", var = avar)
  ## label biplot points by 3-letter CODE (cleaner than bare row numbers)
  rownames(m$q.loading) <- gd$CODE
  m
}
m2000 <- build_year(2000)
m2021 <- build_year(2021)

## ---- Figure 2 (native): simultaneous Q/R biplots, 2000 and 2021 --------------
png(file.path(FIG,"fig2_qr_biplot.png"), width=1600, height=820, res=140)
par(mfrow=c(1,2), mar=c(4.4,4.4,3.2,1))
plot(m2000, main="(a) 2000")
plot(m2021, main="(b) 2021")
dev.off()

## ---- Figure 3 (native): Factor 1 & Factor 2 Index maps, 2021 -----------------
sp21 <- m2021$gisdata
if (inherits(sp21, "sf")) sp21 <- as(sp21, "Spatial")
row.names(sp21) <- as.character(sp21$CODE)

p1 <- spplot(sp21, "index",
             sp.layout = list(label(sp21, var="index", pos=2), label(sp21, pos=3)),
             col.regions = bpy.colors(100), scales = list(draw=TRUE),
             main = "Factor 1 Index — scale of water use (2021)")
p2 <- spplot(sp21, "index2",
             sp.layout = list(label(sp21, var="index2", pos=2), label(sp21, pos=3)),
             col.regions = bpy.colors(100), scales = list(draw=TRUE),
             main = "Factor 2 Index — renewable availability (2021)")
png(file.path(FIG,"fig3_qr_index_maps.png"), width=1650, height=880, res=140)
print(p1, split=c(1,1,2,1), more=TRUE)
print(p2, split=c(2,1,2,1))
dev.off()

## report the index orientation so captions match what is drawn
sp_df <- as.data.frame(sp21)
cat("Factor 1 Index — highest (top 4):\n")
print(head(sp_df[order(-sp_df$index), c("CODE","COUNTRY","index")], 4), row.names=FALSE)
cat("Factor 2 Index — highest (top 4):\n")
print(head(sp_df[order(-sp_df$index2), c("CODE","COUNTRY","index2")], 4), row.names=FALSE)
message("Wrote fig2_qr_biplot.png and fig3_qr_index_maps.png")
