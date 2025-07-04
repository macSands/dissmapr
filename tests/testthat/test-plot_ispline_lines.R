# tests/testthat/test-plot_ispline_lines.R

library(testthat)
library(zetadiv)
library(dplyr)
library(ggplot2)

test_that("plot_ispline_lines returns a ggplot and handles errors", {
  skip_on_cran()

  # 1) Prepare a small run_ispline_models example
  data(bird.spec.fine)
  data(bird.env.fine)

  xy   <- bird.spec.fine[, 1:2]
  spp  <- as.matrix(bird.spec.fine[, 3:102])
  env  <- bird.env.fine[, 3:9]

  out <- run_ispline_models(
    spp_df       = spp,
    env_df       = env,
    xy_df        = xy,
    orders       = 2:3,
    sam          = 10,
    normalize    = "Jaccard",
    reg_type     = "ispline"
  )
  df <- out$ispline_table

  # 2) Successful plot for 'dist'
  p <- plot_ispline_lines(
    ispline_data = df,
    x_var        = "dist",
    orders       = c("Order2", "Order3"),
    cols         = c("red", "blue"),
    shapes       = c(1, 2)
  )
  expect_s3_class(p, "ggplot")
  # Should have three layers: lines + quantile pts + start pts
  expect_length(p$layers, 3)

  # 3) Ambiguous x_var (matches multiple spline cols) triggers error
  expect_error(
    plot_ispline_lines(df, x_var = "temp"),
    "`x_var = 'temp'` must match exactly one"
  )

  # 4) Nonexistent x_var triggers raw-covariate-not-found error
  expect_error(
    plot_ispline_lines(df, x_var = "foobar"),
    "Raw covariate 'foobar' not found"
  )
})
