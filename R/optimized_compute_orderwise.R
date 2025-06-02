#' Compute Order-wise Pairwise Metrics with Optional Parallel Sampling
#'
#' @description
#' A flexible, high-performance engine for calculating dissimilarity,
#' distance, or any custom metric between a *focal* site and every
#' combination of \eqn{(order-1)} other sites.  The routine supports:
#' \itemize{
#'   \item special fast paths for \code{order = 1} (trivial) and
#'         \code{order = 2} when \code{func = orderwise_diss_gower};
#'   \item great-circle distances via \code{geodist_helper}
#'         (requires \code{coord_cols});
#'   \item arbitrary metrics acting on numeric “site vectors” defined by
#'         \code{sp_cols};
#'   \item parallel execution with **future.apply** and progress bars
#'         from **pbapply**;
#'   \item random or complete sampling of higher-order combinations.
#' }
#'
#' @details
#' For each site the relevant data vector is pre-extracted once and
#' re-used, avoiding repeated sub-setting.  When \code{order = 2} and
#' \code{func = orderwise_diss_gower} the full pairwise Gower distance
#' matrix is built in one call to \code{cluster::daisy()}, melted, and
#' returned, yielding a major speed-up over looping.
#'
#' @param df Data frame or data.table containing the site data.
#' @param func Function of two numeric vectors returning a single numeric
#'   value (e.g. \code{orderwise_diss_gower}, \code{geodist_helper}).
#' @param site_col Character; column with the unique site identifiers.
#' @param sp_cols Character or integer vector giving the columns used to
#'   build each site’s numeric vector **when \code{func} is not
#'   \code{geodist_helper}**.
#' @param order Integer scalar or vector of orders (≥ 1).  For values
#'   \eqn{> 2} the function iterates through each order.
#' @param sample_no Integer. Maximum *number* of combinations to sample
#'   **per focal site** for orders ≥ 3. If \code{NULL} (default) all
#'   combinations are used (subject to \code{sample_portion}).
#' @param sample_portion Numeric in (0, 1]. Proportion of combinations to
#'   retain when \code{sample_no = NULL}.  Ignored when
#'   \code{sample_no} is supplied.
#' @param parallel Logical. Run computations in parallel with
#'   **future.apply** (default \code{TRUE}).
#' @param n_workers Integer. Number of background workers; default
#'   \code{parallel::detectCores() - 1}.
#' @param coord_cols Character vector of coordinate columns
#'   (e.g. \code{c("x","y")}) required when
#'   \code{func = geodist_helper}.
#'
#' @return A `data.table` with four columns
#' \describe{
#'   \item{site_from}{Focal site ID.}
#'   \item{site_to}{Comma-separated list of comparison sites
#'         (\code{NA} when \code{order = 1}).}
#'   \item{value}{Computed metric.}
#'   \item{order}{The order used for the calculation.}
#' }
#'
#' @section Dependencies:
#' Requires **pbapply**, **data.table**, **future.apply**, **reshape2**,
#' and **dplyr** (checked at run-time).  For Gower distances the example
#' assumes \code{orderwise_diss_gower()} is exported by *dissmapr*.
#'
#' @importFrom data.table as.data.table rbindlist :=
#' @importFrom pbapply pblapply
#' @importFrom future plan multisession sequential
#' @importFrom stats setNames
#'
#' @keywords internal
#' @examples
#' ## ---------------------------------------------------------------
#' ## Example 1 – Gower dissimilarity between species vectors (order 2)
#' ## ---------------------------------------------------------------
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(1)
#'   n_sites  <- 8
#'   n_spp    <- 20
#'   test_df  <- data.frame(
#'     site_id = paste0("s", seq_len(n_sites)),
#'     x       = runif(n_sites, 23, 24),
#'     y       = runif(n_sites, -34, -33),
#'     matrix(rpois(n_sites * n_spp, 1),
#'            nrow = n_sites,
#'            dimnames = list(NULL, paste0("sp", seq_len(n_spp))))
#'   )
#'
#'   res2 <- optimized_compute_orderwise(
#'     df        = test_df,
#'     func      = orderwise_diss_gower,   # dissmapr helper
#'     site_col  = "site_id",
#'     sp_cols   = paste0("sp", 1:n_spp),
#'     order     = 2,
#'     parallel  = FALSE
#'   )
#'   head(res2)
#' }
#'
#' ## ---------------------------------------------------------------
#' ## Example 2 – Great-circle distance (km) between site centroids
#' ## ---------------------------------------------------------------
#' if (requireNamespace("geosphere", quietly = TRUE)) {
#'   # geodist_helper is expected to be exported by dissmapr
#'   res_geo <- optimized_compute_orderwise(
#'     df         = test_df,
#'     func       = geodist_helper,
#'     site_col   = "site_id",
#'     order      = 2,
#'     coord_cols = c("x", "y"),
#'     parallel   = FALSE
#'   )
#'   head(res_geo)
#' }
#'
optimized_compute_orderwise <- function(df,
                                        func,
                                        site_col,
                                        sp_cols,      # used for non-geodist functions
                                        order = 2,
                                        sample_no = NULL,
                                        sample_portion = 1,  # proportion of combinations to sample (default 1 = 100%)
                                        parallel = TRUE,
                                        n_workers = parallel::detectCores() - 1,
                                        coord_cols = NULL) {  # new parameter for coordinate columns (e.g. c("x", "y"))
  # Required packages check
  required_packages <- c("pbapply", "data.table", "future.apply", "reshape2", "dplyr")
  invisible(lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required but not installed.")
  }))

  suppressPackageStartupMessages({
    library(pbapply)       # progress bars in loops
    library(data.table)    # efficient manipulation
    library(future.apply)  # parallel processing
    library(dplyr)         # data manipulation
    library(reshape2)      # melt() function
  })

  # Convert input data frame to data.table and ensure site IDs are characters.
  dt <- data.table::as.data.table(df)
  dt[[site_col]] <- as.character(dt[[site_col]])
  site_ids <- unique(dt[[site_col]])

  # Precompute each site's "vector":
  # If using geodist_helper, ignore sp_cols and use coord_cols.
  if (identical(func, geodist_helper)) {
    if (is.null(coord_cols)) {
      stop("For geodist_helper, please supply coord_cols.")
    }
    site_vectors <- lapply(site_ids, function(site) {
      vec <- as.numeric(unlist(dt[dt[[site_col]] == site, ..coord_cols, with = FALSE]))
      names(vec) <- coord_cols  # so the helper can extract "x" and "y"
      return(vec)
    })
  } else if (identical(func, orderwise_diss_gower)) {
    # When using orderwise_diss_gower, use sp_cols.
    site_vectors <- lapply(site_ids, function(site) {
      as.numeric(unlist(dt[dt[[site_col]] == site, ..sp_cols, with = FALSE]))
    })
  } else {
    # Otherwise default to sp_cols.
    site_vectors <- lapply(site_ids, function(site) {
      as.numeric(unlist(dt[dt[[site_col]] == site, ..sp_cols, with = FALSE]))
    })
  }
  names(site_vectors) <- site_ids

  start_time <- Sys.time()
  results_list <- list()

  ### --- Order 1: Single-Site Computations --- ###
  if (order == 1) {
    result_dt <- data.table::rbindlist(lapply(site_ids, function(site) {
      # For geodist_helper, we define distance from a site to itself as 0.
      val <- if (identical(func, geodist_helper)) 0 else func(site_vectors[[site]])
      data.table::data.table(site_from = site, site_to = NA, value = val, order = 1L)
    }), use.names = TRUE)
    message("Order 1 (single-site) computations done.")
    return(result_dt)
  }

  ### --- Order 2: Fast Path or Loop over Pairs --- ###
  if (order == 2) {
    if (identical(func, orderwise_diss_gower)) {
      # Fast, vectorized computation for orderwise_diss_gower using sp_cols.
      species_mat <- as.matrix(dt[, ..sp_cols])
      rownames(species_mat) <- dt[[site_col]]  # Set site IDs as row names.
      diss_matrix <- as.matrix(cluster::daisy(species_mat, metric = "gower", stand = FALSE))
      diag(diss_matrix) <- NA  # Exclude self-comparisons.
      melted <- reshape2::melt(diss_matrix,
                               varnames = c("site_from", "site_to"),
                               value.name = "value",
                               na.rm = TRUE)
      melted <- data.table::as.data.table(melted)
      melted[, order := 2L]
      message("Fast vectorized computation used for order 2 with orderwise_diss_gower.")
      return(melted)
    } else if (identical(func, geodist_helper)) {
      # For geodist_helper order 2, loop over pairs.
      comb_data <- expand.grid(site_from = site_ids, site_to = site_ids, stringsAsFactors = FALSE)
      comb_data <- comb_data[comb_data$site_from != comb_data$site_to, ]
      comb_indices <- data.table::as.data.table(comb_data)
      comb_indices[, order := 2L]

      compute_order2 <- function(i) {
        site_from <- comb_indices$site_from[i]
        site_to <- comb_indices$site_to[i]
        value <- func(site_vectors[[site_from]], site_vectors[[site_to]], coord_cols = coord_cols)
        return(value)
      }

      if (parallel && n_workers > 0) {
        future::plan(future::multisession, workers = n_workers)
        values <- pbapply::pblapply(1:nrow(comb_indices), compute_order2)
        future::plan(future::sequential)
      } else {
        values <- pbapply::pblapply(1:nrow(comb_indices), compute_order2)
      }
      comb_indices[, value := unlist(values)]
      results_list[["2"]] <- comb_indices
    } else {
      # Generic looped approach for order 2.
      comb_data <- expand.grid(site_from = site_ids, site_to = site_ids, stringsAsFactors = FALSE)
      comb_data <- comb_data[comb_data$site_from != comb_data$site_to, ]
      comb_indices <- data.table::as.data.table(comb_data)
      comb_indices[, order := 2L]

      compute_order2 <- function(i) {
        site_from <- comb_indices$site_from[i]
        site_to <- comb_indices$site_to[i]
        value <- func(site_vectors[[site_from]], site_vectors[[site_to]])
        return(value)
      }
      if (parallel && n_workers > 0) {
        future::plan(future::multisession, workers = n_workers)
        values <- pbapply::pblapply(1:nrow(comb_indices), compute_order2)
        future::plan(future::sequential)
      } else {
        values <- pbapply::pblapply(1:nrow(comb_indices), compute_order2)
      }
      comb_indices[, value := unlist(values)]
      results_list[["2"]] <- comb_indices
    }
  }

  ### --- Higher Orders (>= 3) --- ###
  for (ord in order) {
    if (ord <= 2) next

    comb_list <- list()
    for (site in site_ids) {
      remaining_sites <- setdiff(site_ids, site)
      total_possible <- choose(length(remaining_sites), ord - 1)
      max_samples <- if (is.null(sample_no)) total_possible else min(sample_no, total_possible)

      if (total_possible <= max_samples) {
        combinations <- combn(remaining_sites, ord - 1, simplify = FALSE)
      } else {
        set.seed(123)
        combinations <- replicate(max_samples, sort(sample(remaining_sites, ord - 1)), simplify = FALSE)
        combinations <- unique(combinations)
        while (length(combinations) < max_samples) {
          additional <- replicate(max_samples - length(combinations),
                                  sort(sample(remaining_sites, ord - 1)),
                                  simplify = FALSE)
          combinations <- unique(c(combinations, additional))
        }
      }

      if (length(combinations) > 0) {
        comb_to <- sapply(combinations, function(x) paste(x, collapse = ","))
        comb_list[[site]] <- data.table::data.table(
          site_from = site,
          site_to = comb_to,
          order = ord
        )
      }
    }
    if (length(comb_list) == 0) next
    comb_indices <- data.table::rbindlist(comb_list, use.names = TRUE)

    compute_higher_order <- function(i) {
      site_from <- comb_indices$site_from[i]
      site_to <- comb_indices$site_to[i]
      sites <- unlist(strsplit(as.character(site_to), ","))
      if (identical(func, geodist_helper)) {
        # For geodist_helper, compute centroid (mean coordinates) of the 'to' sites.
        vec_list <- lapply(sites, function(s) site_vectors[[s]])
        M <- do.call(rbind, vec_list)
        vec_to <- colMeans(M)
        value <- func(site_vectors[[site_from]], vec_to, coord_cols = coord_cols)
      } else {
        vec_list <- lapply(sites, function(s) site_vectors[[s]])
        vec_to <- Reduce("+", vec_list)
        value <- func(site_vectors[[site_from]], vec_to)
      }
      return(value)
    }

    if (parallel && n_workers > 0) {
      future::plan(future::multisession, workers = n_workers)
      values <- pbapply::pblapply(1:nrow(comb_indices), compute_higher_order)
      future::plan(future::sequential)
    } else {
      values <- pbapply::pblapply(1:nrow(comb_indices), compute_higher_order)
    }
    comb_indices[, value := unlist(values)]
    results_list[[as.character(ord)]] <- comb_indices

    elapsed_ord <- difftime(Sys.time(), start_time, units = "secs")
    message(sprintf("Time elapsed for order %d: %.2f seconds", ord, as.numeric(elapsed_ord)))
  }

  if (length(results_list) > 0) {
    final_results <- data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE)
  } else {
    final_results <- NULL
  }

  total_elapsed <- difftime(Sys.time(), start_time, units = "secs")
  message(sprintf("Total computation time: %.2f seconds", as.numeric(total_elapsed)))

  return(final_results)
}
