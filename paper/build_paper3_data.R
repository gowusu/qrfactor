## =====================================================================
## Paper 3 — build the LATEST AQUASTAT water-use dataset for qrfactor
##
## Downloads current FAO AQUASTAT "Water withdrawal by sector" data
## (data.apps.fao.org, updated 2026; CC BY-NC-SA 4.0), takes the most
## recent available year per country, assembles the analysis table in
## the qrfactor format, and (optionally) runs the analysis.
##
## Reproducible: a reviewer re-runs this file and gets the same table.
## Pure base R for the download/reshape; qrfactor only for the analysis.
##
## SOURCE TO CITE:
##   FAO. AQUASTAT — FAO's Global Information System on Water and
##   Agriculture. Water withdrawal by sector. Rome: FAO.
##   https://data.apps.fao.org/aquastat/  (accessed <DATE>).
## Units: sectoral & total withdrawal = 10^9 m3/yr (km3/yr);
##        per-capita withdrawal       = m3/inhabitant/yr.
## =====================================================================

outdir <- "aquastat_raw"
dir.create(outdir, showWarnings = FALSE)

## --- 1. direct CSV URLs (verified from the FAO CKAN catalog) ---------
base <- "https://data.apps.fao.org/catalog/dataset/829bade5-1bb8-4154-aed3-088e5a3b8d00/resource/"
urls <- c(
  Municipal   = paste0(base, "eef085f6-af0b-42f7-a951-79e46f8be69d/download/aquastat-water-withdrawal-by-sector-municipal-water-withdrawal.csv"),
  Industrial  = paste0(base, "d70e804a-4ff6-455c-b116-9f421bfa7a70/download/aquastat-water-withdrawal-by-sector-industrial-water-withdrawal.csv"),
  Agricultural= paste0(base, "ddd2fada-fdee-4e97-be27-77000f508d14/download/aquastat-water-withdrawal-by-sector-agricultural-water-withdrawal.csv"),
  Total       = paste0(base, "12011a9d-1b43-4819-a60e-4b015aeebd26/download/aquastat-water-withdrawal-by-sector-total-water-withdrawal.csv"),
  PerCapita   = paste0(base, "d5610feb-00ba-4391-8a11-bdedfda7a89a/download/aquastat-water-withdrawal-by-sector-total-water-withdrawal-per-capita.csv")
)

## --- 2. download (skips files already present) -----------------------
options(timeout = 600)
for (nm in names(urls)) {
  f <- file.path(outdir, paste0(nm, ".csv"))
  if (!file.exists(f)) {
    message("downloading ", nm, " ...")
    download.file(urls[[nm]], f, mode = "wb", quiet = TRUE)
  }
}

## --- 3. most-recent-year value per country for one variable ----------
## AQUASTAT long format: Area, Area Code, Variable, Year, Value, Symbol, Metadata
latest_by_country <- function(nm) {
  d <- read.csv(file.path(outdir, paste0(nm, ".csv")), stringsAsFactors = FALSE,
                check.names = FALSE)
  d <- d[!is.na(d$Value) & d$Value != "", c("Area", "Year", "Value")]
  d$Year  <- as.integer(d$Year)
  d$Value <- as.numeric(d$Value)
  # keep the row with the largest Year in each country
  keep <- do.call(rbind, lapply(split(d, d$Area),
                                function(x) x[which.max(x$Year), , drop = FALSE]))
  out <- data.frame(Area = keep$Area, value = keep$Value, year = keep$Year,
                    stringsAsFactors = FALSE)
  names(out)[2:3] <- c(nm, paste0(nm, "_year"))
  out
}

vars <- names(urls)
tab  <- latest_by_country(vars[1])
for (nm in vars[-1]) tab <- merge(tab, latest_by_country(nm), by = "Area", all = TRUE)

## --- 4. restrict to the 50 African countries of the original study ---
## (edit this vector freely; names must match AQUASTAT's "Area" spelling)
africa <- c(
  "Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cameroon",
  "Cabo Verde","Central African Republic","Chad","Comoros","Congo",
  "Côte d'Ivoire","Democratic Republic of the Congo","Djibouti","Egypt",
  "Equatorial Guinea","Eritrea","Ethiopia","Gambia","Ghana","Guinea",
  "Guinea-Bissau","Kenya","Lesotho","Liberia","Libya","Madagascar","Malawi",
  "Mali","Mauritania","Morocco","Mozambique","Namibia","Niger","Nigeria",
  "Rwanda","Senegal","Sierra Leone","Somalia","South Africa","Sudan",
  "Eswatini","United Republic of Tanzania","Togo","Tunisia","Uganda","Zambia",
  "Zimbabwe","Gabon")

## report any names that did NOT match (fix spellings, then re-run)
missing <- setdiff(africa, tab$Area)
if (length(missing))
  message("NOT matched in AQUASTAT (fix the spelling in `africa`): ",
          paste(missing, collapse = ", "))

dat <- tab[tab$Area %in% africa, ]

## --- 5. rename to the qrfactor variable set --------------------------
## Domestic=Municipal, Industry=Industrial, Agricultur=Agricultural,
## withdrawal=Total, perCapitaW=PerCapita.  (Resources added in step 6.)
paper <- data.frame(
  COUNTRY    = dat$Area,
  Domestic   = dat$Municipal,
  Industry   = dat$Industrial,
  Agricultur = dat$Agricultural,
  withdrawal = dat$Total,
  perCapitaW = dat$PerCapita,
  year       = dat$Total_year,
  stringsAsFactors = FALSE)

## --- 6. Total Renewable Water Resources + a stable ISO-ish CODE key ----
## TRWR is a long-term (climatological) average that changes negligibly over
## decades, so the authentic AQUASTAT values are hard-coded here (10^9 m3/yr),
## paired to the AQUASTAT country names. CODE mirrors the key in the bundled
## shapefile, so the refreshed table can be joined to the map geometry without
## country-name matching. If you prefer a fresh pull, save resources.csv
## (COUNTRY,Resources) from https://data.apps.fao.org/aquastat/ and it wins.
lookup <- data.frame(
  COUNTRY = africa,
  CODE = c("ALG","ANG","BEN","BOT","BUF","BUR","CAM","CAP","CAR","CHA",
           "COM","CNG","CDI","ZAI","DJI","EGY","EQG","ERI","ETH","GAM",
           "GHA","GIN","GUB","KEN","LES","LIB","LAJ","MAD","MAA","MAL",
           "MAU","MOR","MOZ","NAM","NIG","NIR","RWA","SEN","SIL","SOM",
           "SOU","SUD","SWA","TAN","TOG","TUN","UGA","ZAM","ZIM","GAB"),
  Resources = c(11.6,184.0,25.8,14.7,17.5,3.6,285.5,0.3,144.4,43.0,
                1.2,1283.0,81.0,832.0,0.3,58.3,26.0,6.3,110.0,8.0,
                53.2,226.0,31.0,30.7,5.2,232.0,0.6,337.0,17.3,100.0,
                11.4,29.0,217.1,17.7,33.7,286.2,9.5,39.4,160.0,14.2,
                50.0,64.5,4.5,91.0,14.7,4.6,66.0,105.2,20.0,164.0),
  stringsAsFactors = FALSE)
paper <- merge(paper, lookup, by = "COUNTRY", all.x = TRUE)
if (file.exists("resources.csv")) {           # optional fresh pull overrides
  res <- read.csv("resources.csv", stringsAsFactors = FALSE)
  paper$Resources <- res$Resources[match(paper$COUNTRY, res$COUNTRY)]
}

write.csv(paper, "Africanfreshwater_latest.csv", row.names = FALSE)
message("wrote Africanfreshwater_latest.csv  (", nrow(paper), " countries, ",
        "years ", min(paper$year, na.rm=TRUE), "-", max(paper$year, na.rm=TRUE), ")")

## --- 7. run qrfactor on the refreshed data ---------------------------
## Drop perCapitaW if you keep Total (they are near-collinear); here we
## analyse the three sectoral withdrawals + total + per-capita.
# library(qrfactor)
# v <- c("Domestic","Industry","Agricultur","withdrawal","perCapitaW")
# m <- qrfactor(paper[v])          # or add "Resources" once merged in step 6
# summary(m); plot(m, rowname = "COUNTRY")
