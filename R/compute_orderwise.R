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
