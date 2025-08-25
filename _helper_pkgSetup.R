
install.packages(setdiff(
  c("usethis","devtools","desc","pkgdown","rcmdcheck","testthat","covr",
    "lintr","styler","goodpractice","codemetar","cffr","withr"),
  rownames(installed.packages())
))

# Find files with non-ASCII characters in your package directory:
# tools::showNonASCIIfile("D:/Methods/R/myR_Packages/b-cubed-versions/dissmapr/R")
files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
# Or, to search all R files
lapply(files, tools::showNonASCIIfile)

# All your functions explained in single “reference” page
# Good roxygen2‐style documentation for every function (so that .Rd files get created in man/),
# A site‐builder (we’ll use pkgdown) to turn those .Rd files into an HTML index.

# if you’re in the project root:
devtools::document()

rmarkdown::render("README.Rmd",
                  output_format = rmarkdown::github_document(html_preview = FALSE),
                  output_file   = "README.md")


# in R
roxygen2::roxygenise()
devtools::document()
devtools::check(document = FALSE)   # as you did

devtools::document()
devtools::test()


# if you’re in the project root:
devtools::document()

rmarkdown::render("README.Rmd",
                  output_format = rmarkdown::github_document(html_preview = FALSE),
                  output_file   = "README.md")

# rebuild the site
pkgdown::build_site()

# # If you haven’t already done so - Set up pkgdown
# usethis::use_pkgdown()      # creates a _pkgdown.yml, adds pkgdown to suggests

# Build the site
pkgdown::build_site()
# -Render each .Rd into its own .html under docs/reference/.
# - Produce docs/reference/index.html, with links to all functions (and any groupings you specified in _pkgdown.yml).
# >> Now push your docs/ directory up to GitHub Pages (e.g. via gh-pages), and voilà—you get exactly the style of page you linked.

# Save that as _pkgdown.yml, then rebuild:
devtools::document()
pkgdown::build_site()


# Pkgdown will skip anything that your package’s .Rbuildignore file tells R to ignore.
# So the easiest way is to add _oldREADME.md to .Rbuildignore.
# You can do it by hand—or let usethis do it for you:
usethis::use_build_ignore("_oldREADME.md")

# That writes a line like
# ^_oldREADME\.md$
# into .Rbuildignore.
# Once that’s in place, pkgdown::build_site() will no longer copy or process _oldREADME.md.

# ------------------------------------------------------

# Quick local checks (run in the package root)
# Does testthat scaffolding exist?
file.exists("tests/testthat.R")
length(list.files("tests/testthat", pattern = "^test-.*\\.R$", recursive = TRUE))

# Run tests
devtools::test()          # or: testthat::test_dir("tests/testthat")

# (Optional) Coverage
# install.packages("covr")
covr::report()            # opens HTML report

# If any of these are missing, add them with:
usethis::use_testthat()       # scaffolds tests/
usethis::use_test("my_func")  # creates tests/testthat/test-my_func.R

renv::status()
renv::snapshot()


options(pkgdown.internet = FALSE)         # skip CRAN lookups
pkgdown::build_site(preview = TRUE, new_process = FALSE)


devtools::document()
devtools::install(upgrade = "never", dependencies = TRUE)

pkgdown::build_site(preview = TRUE)  # default new_process = TRUE is OK now
