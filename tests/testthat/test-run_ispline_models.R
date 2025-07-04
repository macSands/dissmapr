# tests/testthat/test-run_ispline_models.R
library(testthat)

test_that("run_ispline_models returns both models and table", {
  skip_on_cran()
  library(zetadiv)

  # load example data
  data(bird.spec.fine)
  data(bird.env.fine)

  # split into coords, species, env
  xy   <- bird.spec.fine[, 1:2]
  spp  <- as.matrix(bird.spec.fine[, 3:102])
  env  <- bird.env.fine[, 3:9]

  # run on a small set of orders and low sam to keep test speed reasonable
  out <- run_ispline_models(
    spp_df       = spp,
    env_df       = env,
    xy_df        = xy,
    orders       = 2:3,
    sam          = 10,
    distance.type= "Euclidean",
    normalize    = "Jaccard",
    reg_type     = "ispline"
  )

  # 1) top-level structure
  expect_type(out, "list")
  expect_named(out, c("zeta_gdm_list", "ispline_table"))

  # 2) zeta_gdm_list
  expect_type(out$zeta_gdm_list, "list")
  expect_length(out$zeta_gdm_list, 2)
  expect_named(out$zeta_gdm_list, c("Order2", "Order3"))
  # each element should be a zetadiv msgdm object
  expect_s3_class(out$zeta_gdm_list[[1]], c("glm", "lm"))

  # 3) ispline_table
  tbl <- out$ispline_table
  expect_s3_class(tbl, "tbl_df")
  # should have at least: all env columns, plus a *_is column for each, plus zOrder
  env_cols    <- names(env)
  ispline_cols<- paste0(env_cols, "_is")
  expect_true(all(env_cols    %in% names(tbl)))
  expect_true(all(ispline_cols%in% names(tbl)))
  expect_true("zOrder"       %in% names(tbl))
  # should tag each row with the correct order
  expect_true(all(tbl$zOrder %in% c("Order2", "Order3")))

  # 4) number of rows: each order returns same number of samples
  n_samples <- nrow(Return.ispline(out$zeta_gdm_list$Order2, env, TRUE)$env)
  expect_equal(nrow(tbl) / 2, n_samples)
})
