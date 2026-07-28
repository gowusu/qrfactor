# Stable structure, growing use: a longitudinal Q/R-mode factor analysis of African water use, 2000–2021

Analysis code and derived data for the paper (George Owusu, *Water Resources Management*).
Everything reproduces from public FAO AQUASTAT data using the [`qrfactor`](../../) R package.

## Run order

Set the working directory to this folder, then run:

1. `build_longitudinal_data.R` — reads the AQUASTAT bulk export (`bulk_eng(in).csv`, downloaded
   from <https://data.apps.fao.org/aquastat/>), extracts the 2000/2010/2021 snapshots for 49
   countries, and writes `Africanfreshwater_longitudinal.csv`.
2. `paper3_longitudinal_diagnostics.R` — all Results quantities (correlations, KMO, Bartlett,
   parallel analysis, varimax loadings, Tucker congruence, silhouette, Moran's I / LISA); writes
   `longitudinal_scores.csv`.
3. `reproduce_qrfactor_native_plots.R` — Figures 2–3 (qrfactor native biplot and Factor-Index maps).
4. `reproduce_longitudinal_figures.R` — Figures 1, 4, 5, 6, 7.
5. `paper3_robustness.R` — the log-transform, Egypt-excluded and k-nearest-neighbour checks (Table 6).

## Files here

| File | Contents |
|---|---|
| `build_longitudinal_data.R` | builds the matched 3-snapshot panel |
| `paper3_longitudinal_diagnostics.R` | every table/figure statistic |
| `reproduce_qrfactor_native_plots.R` | Figures 2–3 |
| `reproduce_longitudinal_figures.R` | Figures 1, 4–7 |
| `paper3_robustness.R` | robustness checks (Table 6) |
| `Africanfreshwater_longitudinal.csv` | the assembled 49-country panel |
| `longitudinal_scores.csv` | per-country factor scores at each date |

The 84 MB AQUASTAT bulk export `bulk_eng(in).csv` is not stored here; download it from AQUASTAT
(link above) into this folder before running step 1.
