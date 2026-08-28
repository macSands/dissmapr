# dissmapr 0.2.2

* Corrected the package version metadata following the 0.2.1 release.
* Synchronized release metadata across `DESCRIPTION`, `CITATION.cff` and
  `.zenodo.json`.
* Improved vignette rendering compatibility for the B-Cubed documentation
  workflow.

# dissmapr 0.2.1

* Updated package documentation and vignettes.
* Rebuilt the pkgdown documentation site.

# dissmapr 0.2.0

* Aligned the package with the
  [B-Cubed software development guide](https://docs.b-cubed.eu/guides/software-development/);
  `R CMD check` now passes with 0 errors, 0 warnings and 0 notes.
* Rewrote the README with a self-contained, runnable minimal example; moved the
  nine tutorials to pkgdown *Articles* and added a "Get started" vignette so the
  documentation compiles cleanly with `b3doc`.
* Replaced whole-namespace `@import` declarations with explicit `@importFrom`
  imports and removed all `library()`/`require()` calls from package code.
* Removed the unused `optimized_compute_orderwise()` development duplicate and
  the leftover `hello()` template function.
* Documented all exported-function arguments and fixed the example data so every
  example runs; added a real `testthat` test for `get_occurrence_data()`.
* Trimmed ~11 MB of unused bundled data (development caches and workspace files).
* Now requires R (>= 4.1.0) (uses the native pipe `|>`).

# dissmapr 0.1.0

* First public release, archived on Zenodo with a citable DOI
  ([10.5281/zenodo.20842434](https://doi.org/10.5281/zenodo.20842434)).
* Core workflow functions: `get_occurrence_data()`, `generate_grid()`,
  `assign_mapsheet()`, `get_enviro_data()`, `format_df()`,
  `compute_orderwise()`, `rm_correlated()`, `predict_dissim()`,
  `map_bioreg()` and `map_bioregDiff()`, plus the ζ-MSGDM helpers
  `run_ispline_models()`, `plot_ispline_lines()` and
  `plot_ispline_boxplots()`.
* Added machine-readable metadata: `codemeta.json` and a `CITATION.cff`
  (Citation File Format 1.2.0) with ORCID and DOI.
* Added a data dictionary (`inst/extdata/data_dictionary.csv`) aligning
  key input/output fields with Darwin Core terms where applicable.
* Documentation website built with `pkgdown`; tutorials available under
  "Articles".
