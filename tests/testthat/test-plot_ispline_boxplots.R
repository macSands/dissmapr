# tests/testthat/test-plot_ispline_boxplots.R

library(testthat)
library(zetadiv)
library(dplyr)
library(tidyr)
library(ggplot2)

test_that("plot_ispline_boxplots returns a ggplot and handles edge cases", {
  skip_on_cran()

  # Prepare minimal ispline_data via run_ispline_models
  data(bird.spec.fine)
  data(bird.env.fine)

  xy   <- bird.spec.fine[,1:2]
  spp  <- as.matrix(bird.spec.fine[,3:102])
  env  <- bird.env.fine[,3:9]

  out <- run_ispline_models(
    spp_df    = spp,
    env_df    = env,
    xy_df     = xy,
    orders    = 2:3,
    sam       = 20,
    normalize = "Jaccard",
    reg_type  = "ispline"
  )
  df <- out$ispline_table

  # 1. Successful boxplot
  p <- plot_ispline_boxplots(
    ispline_data   = df,
    ispline_suffix = "_is",
    order_col      = "zOrder",
    palette        = "cividis",
    direction      = 1,
    ncol           = 2,
    outlier_size   = 0.2
  )
  expect_s3_class(p, "ggplot")
  # Should have one geom_boxplot layer
  geoms <- vapply(p$layers, function(x) class(x$geom)[1], character(1))
  expect_true("GeomBoxplot" %in% geoms)

  # 2. All detected spline columns appear in facets
  spline_cols <- grep("_is$", names(df), value = TRUE)
  facet_vars  <- p$facet$params$facets[[1]]
  expect_true(all(spline_cols %in% facet_vars))

  # 3. Missing order_col triggers error
  df_bad1 <- df %>% rename(other = zOrder)
  expect_error(
    plot_ispline_boxplots(ispline_data = df_bad1),
    "`order_col` must be a column in `ispline_data`"
  )

  # 4. No spline columns triggers error
  df_bad2 <- df %>% select(-ends_with("_is"))
  expect_error(
    plot_ispline_boxplots(ispline_data = df_bad2),
    "no spline columns detected"
  )
})
