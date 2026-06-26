# Changelog

## dissmapr 0.2.0

- Aligned the package with the [B-Cubed software development
  guide](https://docs.b-cubed.eu/guides/software-development/);
  `R CMD check` now passes with 0 errors, 0 warnings and 0 notes.
- Rewrote the README with a self-contained, runnable minimal example;
  moved the nine tutorials to pkgdown *Articles* and added a “Get
  started” vignette so the documentation compiles cleanly with `b3doc`.
- Replaced whole-namespace `@import` declarations with explicit
  `@importFrom` imports and removed all
  [`library()`](https://rdrr.io/r/base/library.html)/[`require()`](https://rdrr.io/r/base/library.html)
  calls from package code.
- Removed the unused `optimized_compute_orderwise()` development
  duplicate and the leftover `hello()` template function.
- Documented all exported-function arguments and fixed the example data
  so every example runs; added a real `testthat` test for
  [`get_occurrence_data()`](https://b-cubed-eu.github.io/dissmapr/reference/get_occurrence_data.md).
- Trimmed ~11 MB of unused bundled data (development caches and
  workspace files).
- Now requires R (\>= 4.1.0) (uses the native pipe `|>`).

## dissmapr 0.1.0

- First public release, archived on Zenodo with a citable DOI
  ([10.5281/zenodo.20842434](https://doi.org/10.5281/zenodo.20842434)).
- Core workflow functions:
  [`get_occurrence_data()`](https://b-cubed-eu.github.io/dissmapr/reference/get_occurrence_data.md),
  [`generate_grid()`](https://b-cubed-eu.github.io/dissmapr/reference/generate_grid.md),
  [`assign_mapsheet()`](https://b-cubed-eu.github.io/dissmapr/reference/assign_mapsheet.md),
  [`get_enviro_data()`](https://b-cubed-eu.github.io/dissmapr/reference/get_enviro_data.md),
  [`format_df()`](https://b-cubed-eu.github.io/dissmapr/reference/format_df.md),
  [`compute_orderwise()`](https://b-cubed-eu.github.io/dissmapr/reference/compute_orderwise.md),
  [`rm_correlated()`](https://b-cubed-eu.github.io/dissmapr/reference/rm_correlated.md),
  [`predict_dissim()`](https://b-cubed-eu.github.io/dissmapr/reference/predict_dissim.md),
  [`map_bioreg()`](https://b-cubed-eu.github.io/dissmapr/reference/map_bioreg.md)
  and
  [`map_bioregDiff()`](https://b-cubed-eu.github.io/dissmapr/reference/map_bioregDiff.md),
  plus the ζ-MSGDM helpers
  [`run_ispline_models()`](https://b-cubed-eu.github.io/dissmapr/reference/run_ispline_models.md),
  [`plot_ispline_lines()`](https://b-cubed-eu.github.io/dissmapr/reference/plot_ispline_lines.md)
  and
  [`plot_ispline_boxplots()`](https://b-cubed-eu.github.io/dissmapr/reference/plot_ispline_boxplots.md).
- Added machine-readable metadata: `codemeta.json` and a `CITATION.cff`
  (Citation File Format 1.2.0) with ORCID and DOI.
- Added a data dictionary (`inst/extdata/data_dictionary.csv`) aligning
  key input/output fields with Darwin Core terms where applicable.
- Documentation website built with `pkgdown`; tutorials available under
  “Articles”.
