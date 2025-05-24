
# All your functions explained in single “reference” page
# Good roxygen2‐style documentation for every function (so that .Rd files get created in man/),
# A site‐builder (we’ll use pkgdown) to turn those .Rd files into an HTML index.

# if you’re in the project root:
devtools::document()

# If you haven’t already done so - Set up pkgdown
usethis::use_pkgdown()      # creates a _pkgdown.yml, adds pkgdown to suggests

# Build the site
pkgdown::build_site()
# -Render each .Rd into its own .html under docs/reference/.
# - Produce docs/reference/index.html, with links to all functions (and any groupings you specified in _pkgdown.yml).
# >> Now push your docs/ directory up to GitHub Pages (e.g. via gh-pages), and voilà—you get exactly the style of page you linked.

