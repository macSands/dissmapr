#' Compute Order-wise Metrics
#'
#' This function computes metrics for ecological data across specified order levels.
#' It supports single-site, pairwise, and higher-order calculations, and allows for
#' parallel processing for efficiency.
#'
#' @param df A data frame containing the ecological data.
#' @param func A function to compute metrics. It must accept inputs in the form
#'   of species vectors or site information depending on the order.
#' @param site_col A character string specifying the column name in `df` representing site IDs.
#' @param sp_cols A vector of column names in `df` representing species data (default: NULL).
#' @param order An integer or vector of integers specifying the order(s) of computation.
#'   - `1`: Single-site computations.
#'   - `2`: Pairwise computations.
#'   - `>= 3`: Higher-order computations.
#' @param sample_no An integer specifying the maximum number of combinations to sample for
#'   higher-order computations (default: NULL for all combinations).
#' @param sample_portion A numeric value between 0 and 1 indicating the proportion of
#'   combinations to sample for higher-order computations (default: 1, meaning 100%).
#' @param parallel A logical value indicating whether to enable parallel computation
#'   (default: TRUE).
#' @param n_workers An integer specifying the number of parallel workers to use
#'   (default: one less than the number of available cores).
#'
#' @return A `data.table` containing the results of computations. Columns include:
#'   - `site_from`: The source site.
#'   - `site_to`: The target site(s) (NA for order = 1).
#'   - `order`: The computation order.
#'   - `value`: The computed metric value.
#'
#' @examples
#' # Example usage with a custom metric function
#' metric_func = function(vec1, vec2 = NULL) {
#'   if (is.null(vec2)) {
#'     return(sum(vec1))  # Example: species richness
#'   } else {
#'     return(sum((vec1 - vec2)^2))  # Example: pairwise dissimilarity
#'   }
#' }
#'
#' data = data.frame(
#'   site = rep(letters[1:3], each = 3),
#'   sp1 = c(1, 0, 2, 1, 2, 0, 0, 1, 1),
#'   sp2 = c(0, 1, 1, 2, 0, 0, 1, 0, 2)
#' )
#'
#'# SPECIES RICHNESS
#'# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'richness = function(vec_from, vec_to = NULL) {
#'  if (is.null(vec_to)) {
#'    # Handle single-site calculations (order = 1)
#'    return(sum(vec_from != 0, na.rm = TRUE))
#'  } else if (length(vec_from) > 1 && length(vec_to) > 1) {
#'    # Handle pairwise or higher-order comparisons
#'    return(abs(sum(vec_from != 0, na.rm = TRUE) - sum(vec_to != 0, na.rm = TRUE)))
#'  } else {
#'    # Invalid input case
#'    return(NA)
#'  }
#'}
#'
#'rich_o12 = compute_orderwise(
#'  df = block_sp,
#'  func = richness,
#'  site_col = 'grid_id',
#'  sp_cols = sp_cols,
#'  # sample_no = 1000,
#'  # sample_portion = 0.5,  # Default is 1 (100%)
#'  order = 1:2,  # Compute for pairwise and higher-order comparisons
#'  parallel = TRUE,
#'  n_workers = 4
#')
#'head(rich_o12)
#'
#' @export
compute_orderwise = function(df,# Optimized Compute_Orderwise Function
                              func,
                              site_col,
                              sp_cols = NULL,
                              order = 2,
                              sample_no = NULL,
                              sample_portion = 1,  # Default is 1 (100%)
                              parallel = TRUE,
                              n_workers = parallel::detectCores() - 1) {

  # Load required packages
  required_packages = c("pbapply", "data.table", "future.apply")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required but not installed.")
  })

  # Load Required Libraries
  suppressPackageStartupMessages({
    library(pbapply)       # For progress bar
    library(data.table)    # For data manipulation
    library(cluster)       # For dissimilarity calculations
    library(future.apply)  # For parallel processing
    library(dplyr)         # For data manipulation
  })

  # Convert to data.table for efficiency
  dt = data.table::as.data.table(df)

  # Ensure `site_col` is character
  dt[[site_col]] = as.character(dt[[site_col]])

  # Get all unique site IDs
  site_ids = unique(dt[[site_col]])

  # Precompute site vectors if sp_cols are provided
  if (!is.null(sp_cols)) {
    site_vectors = lapply(site_ids, function(site) {
      as.numeric(unlist(dt[dt[[site_col]] == site, ..sp_cols, with = FALSE]))
    })
    names(site_vectors) = site_ids
  } else {
    site_vectors = NULL
  }

  # Initialize a list to store results
  results_list = list()

  # Start timing
  start_time = Sys.time()

  # Define Compute Value Function
  compute_value = function(site_from, site_to, ord) {
    if (is.null(sp_cols)) {
      # If sp_cols is NULL, func should handle df, site_col, site_from, site_to directly
      return(func(df, site_col, site_from, site_to))
    }

    vec_from = site_vectors[[site_from]]

    if (ord == 1) {
      # Single-site calculation
      return(func(vec_from))
    } else if (ord == 2) {
      # Pairwise comparison
      vec_to = site_vectors[[site_to]]
      return(func(vec_from, vec_to))
    } else {
      # Higher-order comparison
      sites = unlist(strsplit(as.character(site_to), ","))
      vec_list = lapply(sites, function(s) site_vectors[[s]])
      vec_to = Reduce("+", vec_list)
      return(func(vec_from, vec_to))
    }
  }

  # Iterate over each order value
  if (length(order) > 1) {
    ord_sequence = order
  } else {
    ord_sequence = order
  }

  for (ord in ord_sequence) {
    if (ord == 1) {
      # Single-site calculations
      comb_indices = data.table::data.table(site_from = site_ids, site_to = NA, order = ord)

      # Define a wrapper function for progress bar
      compute_order1 = function(i) {
        site_from = comb_indices$site_from[i]
        value = compute_value(site_from, NA, ord)
        return(value)
      }

      # Parallel Processing Setup
      if (parallel && n_workers > 0) {
        plan(future::multisession, workers = n_workers)
        values = pbapply::pblapply(1:nrow(comb_indices), compute_order1)
        plan(future::sequential)  # Reset to sequential
      } else {
        values = pbapply::pblapply(1:nrow(comb_indices), compute_order1)
      }

      # Assign values
      comb_indices[, value := unlist(values)]

      # Append to results
      results_list[[as.character(ord)]] = comb_indices

    } else if (ord == 2) {
      # Pairwise comparisons
      comb_data = expand.grid(site_from = site_ids, site_to = site_ids, stringsAsFactors = FALSE)
      # Exclude self-pairs
      comb_data = comb_data[comb_data$site_from != comb_data$site_to, ]
      comb_indices = data.table::as.data.table(comb_data)
      comb_indices[, order := ord]

      # Define a wrapper function for progress bar
      compute_order2 = function(i) {
        site_from = comb_indices$site_from[i]
        site_to = comb_indices$site_to[i]
        value = compute_value(site_from, site_to, ord)
        return(value)
      }

      # Parallel Processing Setup
      if (parallel && n_workers > 0) {
        plan(future::multisession, workers = n_workers)
        values = pbapply::pblapply(1:nrow(comb_indices), compute_order2)
        plan(future::sequential)  # Reset to sequential
      } else {
        values = pbapply::pblapply(1:nrow(comb_indices), compute_order2)
      }

      # Assign values
      comb_indices[, value := unlist(values)]

      # Append to results
      results_list[[as.character(ord)]] = comb_indices

    } else if (ord >= 3) {
      # Higher-order comparisons
      comb_list = list()

      for (site in site_ids) {
        remaining_sites = setdiff(site_ids, site)
        total_possible_combinations = choose(length(remaining_sites), ord - 1)
        max_samples = ifelse(is.null(sample_no), total_possible_combinations, sample_no)
        max_samples = min(max_samples, total_possible_combinations)

        if (total_possible_combinations <= max_samples) {
          # Generate all combinations
          combinations = combn(remaining_sites, ord - 1, simplify = FALSE)
        } else {
          # Sample combinations directly
          set.seed(123)  # For reproducibility (optional)
          combinations = replicate(max_samples, sort(sample(remaining_sites, ord - 1)), simplify = FALSE)
          # Remove duplicates
          combinations = unique(combinations)
          # If duplicates were removed and we have fewer than max_samples, resample
          while (length(combinations) < max_samples) {
            additional_combos = replicate(max_samples - length(combinations), sort(sample(remaining_sites, ord - 1)), simplify = FALSE)
            combinations = unique(c(combinations, additional_combos))
          }
        }

        if (length(combinations) > 0) {
          comb_to = sapply(combinations, function(x) paste(x, collapse = ","))
          comb_list[[site]] = data.table::data.table(
            site_from = site,
            site_to = comb_to,
            order = ord
          )
        }
      }

      comb_indices = data.table::rbindlist(comb_list, use.names = TRUE)

      # Define a wrapper function for progress bar
      compute_higher_order = function(i) {
        site_from = comb_indices$site_from[i]
        site_to = comb_indices$site_to[i]
        value = compute_value(site_from, site_to, ord)
        return(value)
      }

      # Parallel Processing Setup
      if (parallel && n_workers > 0) {
        plan(future::multisession, workers = n_workers)
        values = pbapply::pblapply(1:nrow(comb_indices), compute_higher_order)
        plan(future::sequential)  # Reset to sequential
      } else {
        values = pbapply::pblapply(1:nrow(comb_indices), compute_higher_order)
      }

      # Assign values
      comb_indices[, value := unlist(values)]

      # Append to results
      results_list[[as.character(ord)]] = comb_indices
    } else {
      warning(paste("Order", ord, "is not supported. Skipping."))
      next
    }

    # Print elapsed time for this order
    end_time_ord = Sys.time()
    elapsed_time_ord = difftime(end_time_ord, start_time, units = "secs")
    elapsed_minutes_ord = floor(as.numeric(elapsed_time_ord) / 60)
    elapsed_seconds_ord = as.numeric(elapsed_time_ord) %% 60
    cat(sprintf("Time elapsed for order %d: %d minutes and %.2f seconds\n", ord, elapsed_minutes_ord, elapsed_seconds_ord))
  }

  # Combine all results into a single data.table
  final_results = data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE)

  # End timing
  end_time = Sys.time()
  total_elapsed = difftime(end_time, start_time, units = "secs")
  total_minutes = floor(as.numeric(total_elapsed) / 60)
  total_seconds = as.numeric(total_elapsed) %% 60
  cat(sprintf("Total computation time: %d minutes and %.2f seconds\n", total_minutes, total_seconds))

  return(final_results)
}
