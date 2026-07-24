---
title: "Simultaneous Q- and R-Mode Factor Analysis of Water Use Across African Countries: A Reproducible, Map-First Workflow"
author: "George Owusu — Department of Geography and Resource Development, University of Ghana"
date: "2026"
---

**Keywords:** Q-mode factor analysis; R-mode factor analysis; water withdrawal;
spatial ordination; Africa; AQUASTAT; reproducible workflow

---

## Abstract

Africa is described, in the same breath, as *water-rich* and *water-stressed*,
yet the two claims are rarely tested together because the analysis has two
sides — the structure among water-use *variables* and the grouping of *countries*
— usually handled by separate tools. This study applies simultaneous Q- and
R-mode factor analysis, in a single reproducible call, to six water-use
indicators for 50 African countries, and writes the results back to the map. Two
factors explain 74% of the variation: Factor 1 (56%) is an agriculture-driven
scale of overall use; Factor 2 (18%) is a resource-versus-use contrast. Renewable
availability is essentially uncorrelated with use (r between −0.04 and −0.19).
Countries partition into two profiles whose heavy users draw roughly eight times
more agricultural water yet hold less than half the renewable resource (mean 59
versus 128). The contribution is a *workflow*: one function delivers the
ordination, the classification, cluster tests and choropleth maps, so applied
users can move from table to map without stitching five packages together.

---

## 1. Introduction

Water security is among the defining development questions for sub-Saharan and
North Africa, where irrigation dominates abstraction, renewable supply is
unevenly distributed, and per-capita availability is falling with population
growth. A recurring paradox frames the debate: the continent is portrayed as
water-abundant — the Congo Basin alone holds a substantial share of Africa's
renewable freshwater — and simultaneously as one of the most water-stressed
regions on Earth. Both statements can be true, but only if *use* and
*availability* are governed by different processes and different geographies.
Testing that proposition requires a method that classifies **variables** (which
kinds of water use move together?) and **samples** (which countries resemble one
another?) at the same time.

Factor analysis provides both lenses. R-mode analysis groups variables by their
shared variation; Q-mode analysis groups observations by their similarity. In the
natural sciences the two are classically run in tandem — Reyment and Jöreskog
(1993) and Davis (2002) treat Q/R-mode analysis as a standard pairing in geology
and geochemistry — but in practice they are computed separately, in different
software, and reconciled by hand. For an applied water scientist this is a
barrier: the analysis fragments across packages, the results are tables rather
than maps, and reproducibility suffers.

This study has two aims. First, substantively, it characterises the structure of
water use across African countries and tests whether use is decoupled from
availability — the empirical core of the water paradox. Second, methodologically,
it demonstrates a **map-first, single-call workflow** in which simultaneous Q- and
R-mode factor analysis, factor scoring, clustering, cluster testing and choropleth
mapping are produced from one command and written back to the geographic
attribute table. The workflow is implemented in the open-source R package
`qrfactor` (Owusu, 2026) and is fully reproducible from public data.

The remainder of the paper reviews Q/R-mode analysis and spatial ordination
(Section 2), describes the data and workflow (Section 3), reports the ordination,
maps and clusters (Section 4), and discusses the paradox, its policy reading and
the honest limits of a single cross-section (Section 5), before concluding
(Section 6).

## 2. Literature review

**Q- and R-mode factor analysis.** R-mode factor analysis summarises a set of
variables by a smaller set of latent factors, each a weighted combination of the
originals; the weights are *loadings* and the factors are ordered by the variance
(eigenvalue) they capture (Rencher, 2002). Q-mode analysis applies the same
algebra to the transpose, classifying observations. Reyment and Jöreskog (1993)
formalised the joint use of both modes in the natural sciences, and Davis (2002)
made Q/R-mode ordination a staple of quantitative geology. In hydrochemistry the
pairing is routine: R-mode identifies the processes controlling solute variation
while Q-mode groups samples of common origin, typically with varimax rotation and
the Kaiser (1960) eigenvalue-greater-than-one rule for retention, supported by the
scree test (Cattell, 1966).

**Spatial ordination.** Linking ordination to geography has a long lineage.
Wartenberg (1985) introduced multivariate spatial correlation, maximising the
product of variance and spatial autocorrelation; Dray, Saïd and Débias (2008)
generalised it into `multispati`, the basis of modern spatial-ordination tools in
R. This literature is mature, and it is important not to overstate novelty: linking
multivariate structure to spatial pattern is not new. What remains awkward for the
*applied* user is integration — reading a shapefile, running the ordination,
classifying the units, testing the classes and drawing the maps still spans several
packages (`sf` and `sp` for geometry; `stats`/`psych` for factors; `cluster` and
`pvclust` for grouping), each with its own idioms (Bivand, Pebesma & Gómez-Rubio,
2013; Pebesma, 2018).

**The gap and the framework.** Two gaps motivate this study. First, most applied
Q/R-mode work reports variable structure *or* sample grouping, seldom both on a
shared set of axes; and it rarely writes the result back to a map. Second, the
workflow is not reproducible end-to-end from a single object. We therefore adopt a
**pipeline framework** (Table 1): one call ingests a shapefile or table; produces
the correlation structure, eigen-decomposition, R- and Q-mode loadings and scores
on shared axes; derives per-unit indices and clusters; and appends all of these to
the attribute table so that mapping is the final, not a separate, step. The
framework's contribution is integration and reproducibility for applied users, not
a new estimator.

*Table 1. The map-first Q/R workflow: one call, six stages.*

| Stage | Operation | Output |
|---|---|---|
| 1 | Ingest table or shapefile | numeric matrix + geometry |
| 2 | Scale + correlate | correlation/covariance matrix |
| 3 | Eigen-decompose | eigenvalues, variance, PCA loadings |
| 4 | R- and Q-mode loadings/scores | shared-axis biplot |
| 5 | Index, rank, cluster | per-country indices + `kmeans` classes |
| 6 | Write back + map | attribute table + choropleths |

## 3. Data and methods

**Study units and variables.** The analysis covers 50 African countries described
by six water-use indicators derived from FAO's AQUASTAT global water information
system (FAO, 2026): municipal (`Domestic`), industrial (`Industry`) and
agricultural (`Agricultur`) water withdrawal; total renewable water resources
(`Resources`); total water withdrawal (`withdrawal`); and per-capita withdrawal
(`perCapitaW`). Withdrawal quantities follow AQUASTAT conventions
(10⁹ m³ yr⁻¹ for volumes; m³ inhabitant⁻¹ yr⁻¹ for per-capita measures). The
present analysis uses the cross-section bundled with the `qrfactor` package (a
circa-2000 AQUASTAT compilation); a refreshed extract is obtained reproducibly from
the current AQUASTAT dissemination service (Appendix A), and the substantive
findings should be re-confirmed on it (Section 5).

**Analytical workflow.** All six stages of Table 1 are executed by a single call,
`qrfactor()`. Variables are standardised (correlation-based scaling) before
decomposition, so that all indicators contribute equally rather than the
largest-magnitude variable dominating. Factor retention follows the Kaiser (1960)
rule with the scree test (Cattell, 1966) as a check. R- and Q-mode loadings and
scores are placed on a common set of axes to permit joint interpretation. Country
groupings use `kmeans`; group separation on each variable is tested with one-way
ANOVA and the non-parametric Kruskal–Wallis test; cluster stability is assessed by
multiscale bootstrap resampling (`pvclust`; Suzuki & Shimodaira, 2006), and
multivariate outliers by the adjusted-quantile method (Filzmoser, Garrett &
Reimann, 2005). Shapefiles are read with `sf` (Pebesma, 2018) and mapped with `sp`
(Pebesma & Bivand, 2005). Because the factor scores, indices and cluster labels are
written back to the attribute table, the final products are choropleth maps of the
retained factors. Analyses were run in R 4.6 (R Core Team, 2026) with `qrfactor`
1.5; the complete script is given in Appendix A.

## 4. Results

**Correlation structure.** The six indicators show the pattern that underlies the
paradox (Table 2). Agricultural and per-capita withdrawal are nearly
interchangeable (r = 0.99): in Africa, per-capita water use *is* agricultural water
use. Municipal and industrial use rise together (r = 0.75), the signature of a
developing modern economy. Critically, renewable `Resources` correlates only
−0.04 to −0.19 with every use variable: **availability and use are decoupled.**

*Table 2. Pearson correlations among the six water-use indicators (n = 50).*

| | Domestic | Industry | Agricultur | Resources | withdrawal | perCapitaW |
|---|---:|---:|---:|---:|---:|---:|
| Domestic | 1.00 | 0.75 | 0.33 | −0.19 | 0.29 | 0.44 |
| Industry | 0.75 | 1.00 | 0.47 | −0.15 | 0.65 | 0.55 |
| Agricultur | 0.33 | 0.47 | 1.00 | −0.12 | 0.63 | 0.99 |
| Resources | −0.19 | −0.15 | −0.12 | 1.00 | −0.04 | −0.14 |
| withdrawal | 0.29 | 0.65 | 0.63 | −0.04 | 1.00 | 0.65 |
| perCapitaW | 0.44 | 0.55 | 0.99 | −0.14 | 0.65 | 1.00 |

**Factor retention.** The first two eigenvalues exceed unity (3.36, 1.08) and
together explain 74.0% of the variation; a third (0.91) sits just below the Kaiser
threshold, and the sixth collapses to zero — the arithmetic signature of the
`Agricultur`–`perCapitaW` redundancy (Table 3). Two factors are retained.

*Table 3. Eigenvalues and variance explained (correlation-based scaling).*

| Factor | Eigenvalue | Variance (%) | Cumulative (%) |
|---|---:|---:|---:|
| 1 | 3.360 | 56.01 | 56.01 |
| 2 | 1.079 | 17.98 | 73.98 |
| 3 | 0.906 | 15.10 | 89.08 |
| 4 | 0.521 | 8.69 | 97.77 |
| 5 | 0.134 | 2.23 | 100.00 |
| 6 | ~0 | ~0 | 100.00 |

**What the factors mean.** R-mode loadings (Table 4) show that five of six
variables load strongly and negatively on Factor 1 (−0.67 to −0.91): it is a
general *scale of water use*, dominated by agriculture. `Resources` alone loads on
Factor 2 (−0.68), which is therefore a *resource-versus-use contrast*. The
simultaneous Q/R biplot (Figure 1) places variables and countries on these axes
together. It resolves into interpretable regions: the heavy agricultural users —
Egypt, Sudan, Madagascar, Libya, Mali, Mauritania and Swaziland — lie at the
high-use pole; the water-rich Congo and the Democratic Republic of the Congo lie
alone by the `Resources` vector, with vast supply and negligible withdrawal;
Equatorial Guinea, Algeria and Botswana define a domestic/industrial tilt; and a
dense cluster of small, low-use economies fills the remainder. As Figure 1 shows,
Congo and Egypt land almost as far apart as the data allow — one holds the water,
the other uses it.

*Table 4. R-mode loadings on the first two retained factors.*

| Variable | Factor 1 | Factor 2 |
|---|---:|---:|
| Domestic | −0.667 | 0.529 |
| Industry | −0.824 | 0.292 |
| Agricultur | −0.859 | −0.353 |
| Resources | 0.217 | −0.676 |
| withdrawal | −0.791 | −0.249 |
| perCapitaW | −0.909 | −0.264 |

![Figure 1. Simultaneous Q/R-mode biplot. Variables (labelled diamonds) and countries (numbered points) share the first two factor axes (Factor 1 = 56%, Factor 2 = 18%). Country 16 = Egypt, 42 = Sudan, 12 = Congo, 14 = DR Congo, 17 = Equatorial Guinea.](../manual/figures/fig_simultaneous.png)

**The map-first result.** Because the factor indices are written to the attribute
table, the ordination becomes a map. Figure 2 maps the Factor 1 index (overall
water use): the high-use band across Egypt, Sudan, Libya and the Sahel contrasts
with the deep lows of the Congo Basin. Figure 3 maps the Factor 2 index (the
resource contrast), a distinct geography. Together they render the paradox
spatially: the places using the most water are not the places that have it.

![Figure 2. Factor 1 index (overall, agriculture-driven water use) mapped across Africa. Warmer colours denote higher use; the Congo Basin (dark) is lowest.](../manual/figures/fig_map_20.png)

![Figure 3. Factor 2 index (the resource-versus-use contrast) mapped across Africa, showing a geography distinct from Figure 2.](../manual/figures/fig_map_22.png)

**Country profiles.** `kmeans` on the scores yields two robust groups whose mean
profiles name them (Table 5). Cluster B — the heavy users (Egypt, Sudan, Libya,
Mali, South Africa, Morocco and similar) — withdraws roughly eight times more
agricultural water per country than Cluster A, yet sits on **less than half** the
renewable resource (59 versus 128). One-way ANOVA confirms that every *use*
variable separates the clusters at p < 0.001, whereas `Resources` does not
(F = 0.90, p = 0.35): the split is defined by consumption, not endowment, exactly
as Factor 1 implied. Kruskal–Wallis tests agree. The clustering is visualised in
Figure 4.

*Table 5. Cluster mean profiles and one-way ANOVA of group separation.*

| Variable | Cluster A (light) | Cluster B (heavy) | F (1,48) | p |
|---|---:|---:|---:|---:|
| Domestic | 10.7 | 52.9 | 48.79 | 7.8×10⁻⁹ |
| Industry | 3.3 | 17.6 | 44.38 | 2.4×10⁻⁸ |
| Agricultur | 56 | 448 | 63.14 | 2.7×10⁻¹⁰ |
| withdrawal | 1.2 | 14.1 | 16.60 | 1.7×10⁻⁴ |
| perCapitaW | 70 | 519 | 91.52 | 1.1×10⁻¹² |
| **Resources** | **128** | **59** | 0.90 | 0.35 (n.s.) |

![Figure 4. Cluster ordination (clusplot) of the 50 countries on the retained factors, with group hulls. Cluster stability was confirmed by multiscale bootstrap resampling.](../manual/figures/fig_cluster_01.png)

## 5. Discussion

**The paradox, measured.** The analysis converts a rhetorical claim into a
measurement. Water use in Africa is a single dominant axis (56% of the variance),
driven overwhelmingly by agriculture; water *availability* is a nearly independent
second axis (18%); and the two are statistically decoupled (Table 2). The clustering
sharpens the point: the countries drawing the most water are, on average, the least
resource-endowed (Table 5, Figures 2–3). Egypt, Sudan, Libya and Mali abstract
heavily against limited renewable supply, while Congo and the DRC hold abundant
water they scarcely use. For the water paradox, both halves are true — of different
countries.

**Policy reading.** The decoupling implies that continental averages mislead:
"Africa has enough water" is a statement about endowment (Cluster A, the Congo
Basin) that says nothing about the stressed heavy users (Cluster B). Interventions —
irrigation efficiency, transfer infrastructure, demand management — should be
targeted to the Factor 1 high-use, low-resource band that Figures 2–3 delineate,
not applied at a continental scale. The map-first output is directly usable for such
targeting.

**Methodological position.** We deliberately do not claim a new estimator. Linking
multivariate structure to spatial pattern dates to Wartenberg (1985) and is well
served by `multispati` and its descendants (Dray et al., 2008). The contribution
here is *workflow*: simultaneous Q- and R-mode analysis, indices, clustering,
cluster testing and choropleth mapping from one reproducible call, with results
written back to the geography (Table 1). For the applied water scientist that
integration — not a new algorithm — is the practical advance.

**Limitations.** Four are material. (i) *Data vintage.* The cross-section analysed
here is a circa-2000 AQUASTAT compilation; the substantive numbers must be
re-confirmed on the current AQUASTAT extract, obtained reproducibly by the script in
Appendix A, before any policy weight is placed on specific countries. (ii)
*Collinearity.* `Agricultur` and `perCapitaW` are near-duplicates (r = 0.99),
inflating agriculture's role and collapsing the final factor; one of the pair should
be dropped in confirmatory work. (iii) *Cross-sectional inference.* These are
associations in a single snapshot of 50 countries, not causal claims; the drivers of
the decoupling (climate, economy, irrigation policy, data year) are for domain
analysis, not the factor model. (iv) *Spatial autocorrelation not yet tested.* The
maps display clear regional pattern, but we have not formally tested whether the
factor indices are spatially clustered.

**Future work.** The natural next step is a global Moran's I and local (LISA) test of
the factor indices, built from the country adjacency graph (`spdep`), to establish
spatial-clustering significance and identify hot- and cold-spots — a capability
planned for the workflow. A refreshed, multi-year AQUASTAT extract would also permit
tracking whether the use–availability decoupling is widening.

## 6. Conclusion

This study set out to test the African water paradox and to demonstrate a
reproducible, map-first workflow for doing so. Simultaneous Q- and R-mode factor
analysis of six water-use indicators across 50 countries reduces to two
interpretable axes — an agriculture-driven scale of use (56%) and a
resource-versus-use contrast (18%) — that are statistically decoupled, and to two
country profiles whose heavy users hold less than half the renewable resource of
the light users. The paradox is therefore real, and geographic: abundance and
stress describe different countries, mapped in Figures 2 and 3. Methodologically,
the value is integration: one call moved the analysis from table to map, with every
quantity written back to the geography and the whole pipeline reproducible from
public data. As African water planning grows more data-driven, workflows that reach
a map — not just a table — will make the difference between analysis and action.

## Acknowledgments

The author thanks the FAO AQUASTAT programme for open access to the underlying data,
and the R spatial community for the `sf`, `sp`, `cluster` and `pvclust` packages on
which the workflow builds. The `qrfactor` package and the scripts reproducing every
table and figure in this paper are openly available (Appendix A). The author declares
no conflict of interest. *[Add funding/grant details if applicable.]*

## References

*[The list below is a real, curated starting bibliography. Verify each against the
target journal's style, and expand to the journal's expected 30–80 entries — in
particular, add the applied Q/R-mode hydrochemistry studies (e.g. northern Volta
region, Ghana; Nethravathi catchment, India; Cariri Valley, Brazil), whose exact
bibliographic details should be confirmed.]*

1. Bivand, R. S., Pebesma, E. J., & Gómez-Rubio, V. (2013). *Applied Spatial Data Analysis with R* (2nd ed.). Springer.
2. Cattell, R. B. (1966). The scree test for the number of factors. *Multivariate Behavioral Research*, 1(2), 245–276.
3. Davis, J. C. (2002). *Statistics and Data Analysis in Geology* (3rd ed.). Wiley.
4. Dray, S., Saïd, S., & Débias, F. (2008). Spatial ordination of vegetation data using a generalization of Wartenberg's multivariate spatial correlation. *Journal of Vegetation Science*, 19(1), 45–56.
5. FAO (2026). *AQUASTAT — FAO's Global Information System on Water and Agriculture: Water withdrawal by sector.* Rome: Food and Agriculture Organization of the United Nations. https://data.apps.fao.org/aquastat/
6. Filzmoser, P., Garrett, R. G., & Reimann, C. (2005). Multivariate outlier detection in exploration geochemistry. *Computers & Geosciences*, 31(5), 579–587.
7. Kaiser, H. F. (1960). The application of electronic computers to factor analysis. *Educational and Psychological Measurement*, 20(1), 141–151.
8. Owusu, G. (2026). *qrfactor: Simultaneous Q-mode and R-mode factor analysis for spatial data.* R package version 1.5. https://github.com/gowusu/qrfactor
9. Pebesma, E. (2018). Simple Features for R: Standardized support for spatial vector data. *The R Journal*, 10(1), 439–446.
10. Pebesma, E. J., & Bivand, R. S. (2005). Classes and methods for spatial data in R. *R News*, 5(2), 9–13.
11. R Core Team (2026). *R: A Language and Environment for Statistical Computing.* Vienna: R Foundation for Statistical Computing.
12. Rencher, A. C. (2002). *Methods of Multivariate Analysis* (2nd ed.). Wiley.
13. Reyment, R. A., & Jöreskog, K. G. (1993). *Applied Factor Analysis in the Natural Sciences.* Cambridge University Press.
14. Suzuki, R., & Shimodaira, H. (2006). Pvclust: an R package for assessing the uncertainty in hierarchical clustering. *Bioinformatics*, 22(12), 1540–1542.
15. Wartenberg, D. (1985). Multivariate spatial correlation: A method for exploratory geographical analysis. *Geographical Analysis*, 17(4), 263–283.

## Appendix A — Reproducibility

The complete analysis reproduces from public data. The refreshed dataset is built by
`paper/build_paper3_data.R`, which downloads the current AQUASTAT withdrawal series
and assembles the analysis table; the analysis itself is a single call:

```r
library(qrfactor)
# `dat` = the AQUASTAT table (columns: Domestic, Industry, Agricultur,
#          Resources, withdrawal, perCapitaW) for the 50 countries
m <- qrfactor(dat[c("Domestic","Industry","Agricultur",
                    "Resources","withdrawal","perCapitaW")])
summary(m)                     # Tables 3 and 4
m$correlation                  # Table 2
plot(m, rowname = "COUNTRY")   # Figure 1 (biplot)

# spatial version: scores/indices/clusters written back to the map
sm <- qrfactor(shapefile_dir, layer = "Africanfreshwater",
               var = c("Domestic","Industry","Agricultur",
                       "Resources","withdrawal","perCapitaW"))
plot(sm, plot = "map")         # Figures 2-3 (factor-index choropleths)
plot(sm, plot = "cluster")     # Figure 4 (clusplot) + Table 5 (ANOVA)
```
