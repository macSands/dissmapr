test_that("get_occurrence_data() harmonises a simple in-memory data.frame", {
  df <- data.frame(
    site_id = 1:4,
    x       = c(18.4, 18.5, 18.6, 18.7),
    y       = c(-33.9, -33.8, -33.7, -33.6),
    sp_name = c("sp_a", "sp_b", "sp_a", "sp_c"),
    pa      = c(1, 1, 1, 1)
  )

  out <- get_occurrence_data(data = df, source_type = "data_frame")

  expect_s3_class(out, "data.frame")
  expect_true(all(c("site_id", "x", "y", "sp_name") %in% names(out)))
  expect_equal(nrow(out), nrow(df))
})

test_that("get_occurrence_data() errors on an invalid source_type", {
  df <- data.frame(x = 1, y = 1, sp_name = "sp_a", pa = 1)
  expect_error(get_occurrence_data(data = df, source_type = "not_a_type"))
})

test_that("get_occurrence_data() rejects out-of-range coordinates", {
  df <- data.frame(x = 999, y = 0, sp_name = "sp_a", pa = 1)
  expect_error(
    get_occurrence_data(data = df, source_type = "data_frame"),
    "Longitude"
  )
})
