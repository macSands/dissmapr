#' Helper functions for compute_orderwise
#'
#' A suite of internal helper functions to compute ecological indices used by
#' `compute_orderwise()`, including geographic distance, dissimilarity metrics,
#' correlations, and mutual information.
#'
#' @name required_packages
#' @keywords internal
#' @importFrom data.table as.data.table rbindlist
#' @importFrom geosphere distHaversine
#' @importFrom vegan vegdist
#' @importFrom cluster daisy
#' @importFrom reshape2 melt
#' @importFrom entropy mi.plugin
#' @keywords internal
NULL

# Verify required packages for helpers
required_packages <- c("geosphere", "vegan", "cluster", "reshape2", "entropy")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
}))

# DISTANCE BETWEEN SITES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate geographic distance via Haversine formula (helper)
#'
#' Computes the geographic (Haversine) distance in meters between two coordinate
#' pairs, or returns 0 for single-site (order=1) calculations when `vec_to` is `NULL`.
#'
#' @param vec_from Numeric vector of length 2, or a named vector containing
#'   coordinates. If named, should include names given by `coord_cols`.
#' @param vec_to Optional numeric vector of length 2 (destination coordinates);
#'   if `NULL` (default), the function returns 0.
#' @param coord_cols Character vector of length 2 giving the names of the longitude
#'   and latitude elements in `vec_from`/`vec_to` when those are named vectors.
#'   Defaults to `c("x", "y")`.
#'
#' @return Numeric distance in meters between the two points, or 0 if `vec_to` is `NULL`.
#'
#' @importFrom geosphere distHaversine
#' @export
#' @keywords internal
#'
geodist_helper <- function(vec_from, vec_to = NULL, coord_cols = c("x", "y")) {
  # For single-site (order 1) calculations, return 0.
  if (is.null(vec_to)) {
    return(0)
  }

  # Extract coordinates from named vectors or assume order if unnamed
  if (!is.null(names(vec_from)) && all(coord_cols %in% names(vec_from))) {
    vec_from_coords <- as.numeric(vec_from[coord_cols])
  } else {
    vec_from_coords <- as.numeric(vec_from)
  }

  if (!is.null(names(vec_to)) && all(coord_cols %in% names(vec_to))) {
    vec_to_coords <- as.numeric(vec_to[coord_cols])
  } else {
    vec_to_coords <- as.numeric(vec_to)
  }

  # Validate coordinate lengths
  if (length(vec_from_coords) != 2 || length(vec_to_coords) != 2) {
    stop("Both input vectors must be coordinate pairs of length 2 according to coord_cols: ",
         paste(coord_cols, collapse = ", "))
  }

  # Compute and return Haversine distance
  geosphere::distHaversine(vec_from_coords, vec_to_coords)
}




# !! CURRENTLY DOESN'T WORK WITH `compute_orderwise` !!

#' Calculate Distance Between Sites
#'
#' @param df A data frame containing site coordinates.
#' @param site_col The column name representing site IDs.
#' @param vec_from The site ID or vector for the starting site.
#' @param vec_to The site ID or vector for the destination site(s).
#'
#' @return Distance(s) in meters between the specified sites.
#' @export
#' @keywords internal
#'
distance <- function(df, site_col, vec_from, vec_to = NULL) {
  library(geosphere)

  # Subset 'from' data
  site_from_data <- df[df[[site_col]] == vec_from, c("x", "y"), drop = FALSE]

  # Handle missing or invalid inputs
  if (nrow(site_from_data) == 0) stop("Invalid 'from' site ID: ", vec_from)

  if (is.null(vec_to)) {
    stop("Order = 1 calculations are not supported for the distance function.")
  } else if (length(vec_to) == 1) {
    # Pairwise comparison
    site_to_data <- df[df[[site_col]] == vec_to, c("x", "y"), drop = FALSE]
    if (nrow(site_to_data) == 0) stop("Invalid 'to' site ID: ", vec_to)
    return(distHaversine(site_from_data, site_to_data))
  } else {
    # Higher-order comparisons
    site_to_data <- df[df[[site_col]] %in% vec_to, c("x", "y"), drop = FALSE]
    if (nrow(site_to_data) == 0) stop("Invalid 'to' site IDs: ", paste(vec_to, collapse = ", "))
    distances <- apply(site_to_data, 1, function(to_coords) {
      distHaversine(site_from_data, to_coords)
    })
    # Ensure output is scalar
    return(sum(distances, na.rm = TRUE))
  }
}

# !! WORKS AS STANDALONE FUNCTION I.E. NOT WITH `compute_orderwise` !!
# Function to calculate pairwise distances using distm
#' Calculate Pairwise Distance Matrix
#'
#' @param data A data frame containing site coordinates and IDs.
#' @param distance_fun The distance function to use (default: `distGeo`).
#'
#' @return A data frame containing pairwise distances between sites.
#' @export
#' @keywords internal
#'
calculate_pairwise_distances_matrix <- function(data, distance_fun = distGeo) {
  if (!all(c("grid_id", "x", "y") %in% colnames(data))) {
    stop("Data frame must contain 'grid_id', 'x', and 'y' columns.")
  }

  distance_matrix <- distm(data[, c("x", "y")], fun = distance_fun) / 1000  # Convert to km
  distances <- as.data.frame(as.table(distance_matrix))
  colnames(distances) <- c("site_from_index", "site_to_index", "value")

  distances$site_from <- data$grid_id[distances$site_from_index]
  distances$site_to <- data$grid_id[distances$site_to_index]

  distances <- distances[distances$site_from != distances$site_to, ]
  return(distances[, c("site_from", "site_to", "value")])
}

# SPECIES RICHNESS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate Species Richness
#'
#' @param vec_from A numeric vector representing species counts at the first site.
#' @param vec_to (Optional) A numeric vector for pairwise or higher-order comparisons.
#'
#' @return The species richness value.
#' @export
#' @keywords internal
#'
richness <- function(vec_from, vec_to = NULL) {
  if (is.null(vec_to)) {
    return(sum(vec_from != 0, na.rm = TRUE))
  } else if (length(vec_from) > 1 && length(vec_to) > 1) {
    return(abs(sum(vec_from != 0, na.rm = TRUE) - sum(vec_to != 0, na.rm = TRUE)))
  } else {
    return(NA)
  }
}


# SPECIES TURNOVER (BETA DIVERSITY)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate Species Turnover or Beta Diversity
#'
#' @param vec_from A numeric vector representing species counts at the first site.
#' @param vec_to (Optional) A numeric vector for pairwise or higher-order comparisons.
#'
#' @return The species turnover value.
#' @export
#' @keywords internal
#'
turnover <- function(vec_from, vec_to = NULL) {
  if (is.null(vec_to)) {
    stop("Turnover calculation requires both vec_from and vec_to.")
  } else if (length(vec_from) > 1 && length(vec_to) > 1) {
    # Calculate species turnover
    total_species <- sum((vec_from != 0 | vec_to != 0), na.rm = TRUE) # Identifies all species present in either vec_from or vec_to
    shared_species <- sum((vec_from != 0 & vec_to != 0), na.rm = TRUE) # Identifies species shared between vec_from and vec_to
    turnover <- (total_species - shared_species) / total_species # calculates the proportion of species not shared relative to the total number of species present
    return(turnover)
  } else {
    return(NA)
  }
}

# Min: 0.0000: Indicates complete similarity (no turnover).
# >> All species present at site_from are also present at site_to, and vice versa.
# Max: 1.0000: Indicates complete turnover (no shared species).
# >> All species present at site_from are absent at site_to and vice versa.

# 1st Qu.: 0.9778, Median: 1.0000, 3rd Qu.: 1.0000:
# >> The species turnover is very high in most site pairs, as the majority of values
# are close to or equal to 1.
# This suggests that most sites have few or no shared species, which could indicate
# high species heterogeneity across the landscape.
# Mean: 0.9807:
# >> The average turnover across all site pairs is approximately 98%.
# >> This reinforces the observation that most sites have a high degree of dissimilarity
# >> in their species composition.

# SPECIES ABUNDANCE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate Species Abundance
#'
#' @param vec_from A numeric vector representing species counts at the first site.
#' @param vec_to (Optional) A numeric vector for pairwise or higher-order comparisons.
#'
#' @return The species abundance value.
#' @export
#' @keywords internal
#'
abund <- function(vec_from, vec_to = NULL) {
  if (is.null(vec_to)) {
    return(sum(vec_from, na.rm = TRUE))
  } else if (length(vec_from) > 1 && length(vec_to) > 1) {
    return(abs(sum(vec_from, na.rm = TRUE) - sum(vec_to, na.rm = TRUE)))
  } else {
    return(NA)
  }
}

# PHI COEFFICIENT (PRESENCE-ABSENCE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate Phi Coefficient
#'
#' @param vec_from A binary presence-absence vector for the first site.
#' @param vec_to A binary presence-absence vector for the second site.
#'
#' @return The Phi Coefficient value.
#' @export
#' @keywords internal
#'
phi_coef <- function(vec_from, vec_to) {
  data_i <- as.numeric(vec_from > 0)
  data_j <- as.numeric(vec_to > 0)

  A <- sum(data_i == 1 & data_j == 1)
  B <- sum(data_i == 1 & data_j == 0)
  C <- sum(data_i == 0 & data_j == 1)
  D <- sum(data_i == 0 & data_j == 0)

  denominator <- sqrt((A + B) * (A + C) * (B + D) * (C + D))
  if (is.na(denominator) || denominator == 0) {
    return(NA)
  }

  return((A * D - B * C) / denominator)
}

# Phi Coefficient (Presence-Absence Data):
# Measures the strength of association between species pairs,
# ranging from -1 (perfect negative association) to +1 (perfect positive association)

# --> Intrepretting plots:
# -1 = perfect negative association (species never co-occur)
# 0 = no association (species co-occur randomly).
# +1 = perfect positive association (species always co-occur)
# In plot, if mean Phi values are all very close to 0 = on average,
# there is little to no strong co-occurrence signal across sites.

# SPEARMAN'S RANK CORRELATION (ABUNDANCES)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate Spearman's Rank Correlation
#'
#' @param vec_from A numeric vector representing species abundances at the first site.
#' @param vec_to A numeric vector representing species abundances at the second site.
#'
#' @return Spearman's rank correlation coefficient.
#' @export
#' @keywords internal
#'
cor_spear <- function(vec_from, vec_to) {
  if (length(vec_from) > 1 && length(vec_to) > 1) {
    return(cor(vec_from, vec_to, method = "spearman", use = "pairwise.complete.obs"))
  } else {
    return(NA)
  }
}

# Spearman’s Rank Correlation (Abundance Data):
# Measures the rank-based association between species pairs, also ranging from -1 to +1.
# >> Description: Measures the strength and direction of a monotonic association between two species' abundances.
# >> When to Use: When data is non-parametric/doesn't meet assumptions of normality.
# It is based on ranks i.e. is robust to outliers.
# >> Interpretation: Ranges from −1 (perfect negative correlation) to
# 1 (perfect positive correlation), with 0 indicating no association.

# --> Intrepretting plots:
# High Mean Spearman Values (Red) = strong positive associations between
# species abundances at the site (e.g., species tend to have similar abundance patterns).
# Low Mean Spearman Values (Green) = weak or negative associations
# (e.g., species have dissimilar abundance patterns).
# Near-Zero Mean Spearman Values = species abundances not correlated at site
# (random or neutral associations).

# PEARSON'S CORRELATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate Pearson's Correlation
#'
#' @param vec_from A numeric vector representing species abundances at the first site.
#' @param vec_to A numeric vector representing species abundances at the second site.
#'
#' @return Pearson's correlation coefficient.
#' @export
#' @keywords internal
#'
cor_pears <- function(vec_from, vec_to) {
  if (length(vec_from) > 1 && length(vec_to) > 1) {
    return(cor(vec_from, vec_to, method = "pearson", use = "pairwise.complete.obs"))
  } else {
    return(NA)
  }
}

# >> Description: Measures the linear association between two species' abundances.
# >> When to Use: When the data is normally distributed and you are interested in
# linear relationships.
# >> Interpretation: Similar to Spearman’s, it ranges from −1 to 1.


# BRAY-CURTIS DISSIMILARITY
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate Bray-Curtis Dissimilarity
#'
#' @param vec_from A numeric vector representing species abundances at the first site.
#' @param vec_to A numeric vector for species abundances at the second site.
#'
#' @return The Bray-Curtis dissimilarity value.
#' @export
#' @keywords internal
#'
diss_bcurt <- function(vec_from, vec_to) {
  if (length(vec_from) > 1 && length(vec_to) > 1) {
    return(vegan::vegdist(rbind(vec_from, vec_to), method = "bray")[1])
  } else {
    return(NA)
  }
}

# >> Description: Measures the dissimilarity between two samples based on species abundances.
# Ranges from 0 (identical) to 1 (completely dissimilar).
# >> When to Use: Quantify dissimilarity based on abundance data, taking into account
# the differences in species counts.
# >> Interpretation: Close to 0 indicates high similarity, value close to
# 1 indicates high dissimilarity.

# GOWER'S SIMILARITY
# !! CURRENTLY DOESN'T WORK WITH `compute_orderwise` !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate Gower's Dissimilarity
#'
#' @param vec_from A numeric or categorical vector for the first site.
#' @param vec_to A numeric or categorical vector for the second site.
#'
#' @return The Gower dissimilarity value.
#' @export
#' @keywords internal
#'
gower_dissimilarity <- function(vec_from, vec_to) {
  if (length(vec_from) > 1 && length(vec_to) > 1) {
    comb_df <- data.frame(t(rbind(vec_from, vec_to)))
    diss_mat <- as.matrix(cluster::daisy(comb_df, metric = "gower"))
    return(diss_mat[1, 2])
  } else {
    return(NA)
  }
}


orderwise_diss_gower <- function(vec_from, vec_to = NULL) {
  # If the second vector is missing, we cannot compute a dissimilarity between two groups.
  # (For order=1 this function returns NA because comparing a single site to nothing isn’t meaningful.)
  if (is.null(vec_to)) {
    return(NA_real_)
  }

  # Combine the two input vectors (vec_from and vec_to) into a matrix.
  # Each row represents one group (or one site vector).
  M <- rbind(vec_from, vec_to)

  # Use daisy() from the cluster package to compute the Gower dissimilarity.
  # Gower's metric works well with mixed data types. Here we assume that vec_from
  # and vec_to are numeric vectors (for example, species counts or aggregated values)
  # and we set stand = FALSE so that no standardization occurs in daisy().
  diss <- cluster::daisy(M, metric = "gower", stand = FALSE)

  # The daisy() function returns a "dissimilarity" object.
  # For a two-row matrix, it contains one dissimilarity value.
  # We convert that to a plain numeric value and return it.
  return(as.numeric(diss))
}

# >> Description: Versatile measure that can handle both continuous and categorical data.
# Calculates similarity between two samples based on attributes.
# >> When to Use: When your data includes different types of variables
# (e.g., abundance, presence/absence, and categorical data)
# >> Interpretation: Ranges from 0 (no similarity) to 1 (complete similarity).

# Low Dissimilarity (Close to 0): The two vectors (vec_from and vec_to) are very similar
# For numeric attributes, their values are close, and categorical attributes match frequently.
# High Dissimilarity (Close to 1): The two vectors are very different
# Intermediate Values (Between 0 and 1)


# Function to calculate pairwise Gower dissimilarities
#' Calculate Pairwise Gower Dissimilarity Matrix
#'
#' @param df A data frame containing site information.
#' @param sp_cols A vector of column names for species data.
#'
#' @return A data frame containing pairwise Gower dissimilarities between sites.
#' @export
#' @keywords internal
#'
calculate_pairwise_gower_dist_matrix <- function(df, sp_cols) {
  sbs_gower_df <- as.data.frame(as.matrix(cluster::daisy(df[, sp_cols], metric = "gower", stand = FALSE)))
  sbs_gower_df$site_from <- row.names(sbs_gower_df)
  sbs_gower_melt <- reshape2::melt(sbs_gower_df, id.vars = "site_from", variable.name = "site_to", value.name = "value")

  # Exclude self-pairs
  sbs_gower_melt <- sbs_gower_melt[sbs_gower_melt$site_from != sbs_gower_melt$site_to, ]
  sbs_gower_melt$x <- df$x[match(sbs_gower_melt$site_from, df$grid_id)]
  sbs_gower_melt$y <- df$y[match(sbs_gower_melt$site_from, df$grid_id)]

  return(sbs_gower_melt)
}


# MUTUAL INFORMATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate Mutual Information
#'
#' @param vec_from A numeric or categorical vector for the first variable.
#' @param vec_to A numeric or categorical vector for the second variable.
#'
#' @return The mutual information value between the two variables.
#' @export
#' @keywords internal
#'
mutual_info <- function(vec_from, vec_to) {
  library(entropy)
  if (length(vec_from) > 1 && length(vec_to) > 1) {
    joint_dist <- table(vec_from, vec_to)
    mi <- mi.plugin(joint_dist)
    return(mi)
  } else {
    return(NA)
  }
}

# >> Description: A non-parametric measure of mutual dependence between two variables
# i.e. captures both linear and non-linear associations
# >> When to Use: When you suspect non-linear relationships or want a flexible measure of association.
# >> Interpretation: Higher values indicate stronger associations.
# Does not have fixed range but is always non-negative.

# Non-Negative Values (≥ 0)
# A value of 0 means the two variables are independent (no shared information).
# Higher values indicate greater dependency or shared information between the two variables.
