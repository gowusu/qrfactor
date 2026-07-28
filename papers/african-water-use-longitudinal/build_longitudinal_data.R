## =====================================================================
## build_longitudinal_data.R
## Builds the LONGITUDINAL AQUASTAT water-use panel for qrfactor Paper 3.
##
## Design: three decadal snapshots — 2000, 2010, 2021 — for the African
## countries, taking the observation nearest each benchmark within a +/-3
## year window. Four analysis variables per snapshot:
##   Domestic   = Municipal water withdrawal        [10^9 m3/yr]
##   Industry   = Industrial water withdrawal       [10^9 m3/yr]
##   Agricultur = Agricultural water withdrawal     [10^9 m3/yr]
##   Resources  = Total renewable water resources   [10^9 m3/yr]  (time-invariant)
## plus Total water withdrawal (reported, = sum of the 3 sectors) for context.
##
## Countries are matched to AQUASTAT by their official numeric area code
## (M49-style REF_AREA), verified against the bulk file's own (code,name)
## table — NOT by hand-typed guesses — so no country can silently receive
## another's data.
##
## SOURCE: FAO AQUASTAT bulk export (bulk_eng(in).csv), 1960-2022,
##   downloaded from https://data.apps.fao.org/aquastat/ .
## Everything (incl. renewable resources) comes from this one file, so the
## panel is internally consistent and reproducible from a single source.
## =====================================================================

## Run this script with the working directory set to the folder that holds the
## data file bulk_eng(in).csv (in RStudio: Session > Set Working Directory >
## To Source File Location, or edit the path below).
pdir <- "."
setwd(pdir)

b <- read.csv("bulk_eng(in).csv", stringsAsFactors = FALSE, check.names = FALSE)
names(b)[1:9] <- c("elemCode","element","refArea","area","year","year2","value","f1","f2")
b$year    <- suppressWarnings(as.integer(b$year))
b$value   <- suppressWarnings(as.numeric(b$value))
b$refArea <- suppressWarnings(as.integer(b$refArea))

## ---- verified crosswalk: COUNTRY, CODE (shapefile key), AQUASTAT refArea ----
## refArea codes read directly from bulk_eng's own (REF_AREA, AREA) table.
xw <- read.table(header = TRUE, sep = "|", strip.white = TRUE, quote = "",
                 stringsAsFactors = FALSE, text =
"COUNTRY|CODE|refArea
Algeria|ALG|12
Angola|ANG|24
Benin|BEN|204
Botswana|BOT|72
Burkina Faso|BUF|854
Burundi|BUR|108
Cameroon|CAM|120
Cabo Verde|CAP|132
Central African Republic|CAR|140
Chad|CHA|148
Comoros|COM|174
Congo|CNG|178
Cote d'Ivoire|CDI|384
Democratic Republic of the Congo|ZAI|180
Djibouti|DJI|262
Egypt|EGY|818
Equatorial Guinea|EQG|226
Eritrea|ERI|232
Ethiopia|ETH|231
Gambia|GAM|270
Ghana|GHA|288
Guinea|GIN|324
Guinea-Bissau|GUB|624
Kenya|KEN|404
Lesotho|LES|426
Liberia|LIB|430
Libya|LAJ|434
Madagascar|MAD|450
Malawi|MAA|454
Mali|MAL|466
Mauritania|MAU|478
Morocco|MOR|504
Mozambique|MOZ|508
Namibia|NAM|516
Niger|NIG|562
Nigeria|NIR|566
Rwanda|RWA|646
Senegal|SEN|686
Sierra Leone|SIL|694
Somalia|SOM|706
South Africa|SOU|710
Sudan|SUD|729
Eswatini|SWA|748
United Republic of Tanzania|TAN|834
Togo|TOG|768
Tunisia|TUN|788
Uganda|UGA|800
Zambia|ZAM|894
Zimbabwe|ZIM|716
Gabon|GAB|266")

b <- b[b$refArea %in% xw$refArea, ]
b$COUNTRY <- xw$COUNTRY[match(b$refArea, xw$refArea)]
b$CODE    <- xw$CODE   [match(b$refArea, xw$refArea)]

elems <- c(Domestic   = "Municipal water withdrawal [10^9 m3/year]",
           Industry   = "Industrial water withdrawal [10^9 m3/year]",
           Agricultur = "Agricultural water withdrawal [10^9 m3/year]",
           withdrawal = "Total water withdrawal [10^9 m3/year]",
           Resources  = "Total renewable water resources [10^9 m3/year]")

nearest <- function(varElem, target, win = 3) {
  s <- b[b$element == varElem & !is.na(b$value) & !is.na(b$year) &
           b$year >= target - win & b$year <= target + win, ]
  if (!nrow(s)) return(data.frame(CODE=character(), val=numeric(), yr=integer()))
  do.call(rbind, lapply(split(s, s$CODE), function(x){
    i <- which.min(abs(x$year - target))
    data.frame(CODE = x$CODE[1], val = x$value[i], yr = x$year[i],
               stringsAsFactors = FALSE)
  }))
}

benchmarks <- c(2000, 2010, 2021)
panel <- xw[, c("COUNTRY","CODE")]
for (ty in benchmarks) for (v in names(elems)) {
  n <- nearest(elems[[v]], ty, win = 3)
  panel[[paste0(v, "_", ty)]]    <- n$val[match(panel$CODE, n$CODE)]
  panel[[paste0(v, "_yr_", ty)]] <- n$yr [match(panel$CODE, n$CODE)]
}

## ---- matched panel: complete on the 4 analysis vars at all 3 snapshots ------
avar <- c("Domestic","Industry","Agricultur","Resources")
need <- as.vector(outer(avar, benchmarks, function(a,b) paste0(a,"_",b)))
ok   <- complete.cases(panel[need])
cat("countries complete on all 4 vars x 3 snapshots:", sum(ok), "of", nrow(panel), "\n")
if (any(!ok)) cat("dropped (incomplete panel):", paste(panel$COUNTRY[!ok], collapse=", "), "\n")
panel <- panel[ok, ]

for (ty in benchmarks) {
  yrs <- panel[[paste0("withdrawal_yr_", ty)]]
  cat(sprintf("benchmark %d: actual years used  min=%d median=%g max=%d\n",
              ty, min(yrs), median(yrs), max(yrs)))
}

cat("\nRenewable resources spot-check (AQUASTAT: Congo~832, DR Congo~1283):\n")
print(panel[panel$COUNTRY %in% c("Congo","Democratic Republic of the Congo"),
            c("COUNTRY","Resources_2021")], row.names = FALSE)

write.csv(panel, "Africanfreshwater_longitudinal.csv", row.names = FALSE)
cat("\nwrote Africanfreshwater_longitudinal.csv (", nrow(panel), " countries )\n", sep="")

if (requireNamespace("sf", quietly = TRUE)) {
  gd <- sf::st_read(system.file("external", package = "qrfactor"),
                    layer = "Africanfreshwater", quiet = TRUE)
  gd$CODE <- trimws(as.character(gd$CODE))
  miss <- setdiff(panel$CODE, gd$CODE)
  cat("panel CODEs missing from shapefile:", if(length(miss)) paste(miss,collapse=", ") else "none", "\n")
}
cat("DONE\n")
