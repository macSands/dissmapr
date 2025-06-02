#' Map Bioregional Change Metrics Between Categorical Raster Layers
#'
#' @description
#' Calculates five complementary indices that quantify how the categorical
#' (e.g. bioregion / cluster) label of each raster cell changes across a
#' temporal or scenario stack of rasters:
#'
#' * **Difference count** – total number of times a cell’s label differs from
#'   the first layer.
#' * **Shannon entropy** – information‐theoretic diversity of labels within the
#'   cell’s time-series.
#' * **Stability** – proportion of layers in which the label is identical to
#'   the first layer *(1 = always unchanged, 0 = always different)*.
#' * **Transition frequency** – sum of binary change maps between successive
#'   layers (*how often a change occurs between any pair of neighbours*).
#' * **Weighted change index** – cumulative dissimilarity-weighted change where
#'   the weight is derived from the empirical frequency of transitions between
#'   all pairs of labels.
#'
#' @details
#' The dissimilarity weights for the **weighted change index** are built from
#' the observed transition table of successive layers, normalised to lie
#' between 0 and 1 (larger = rarer transition).  The function accepts either a
#' multi-layer `SpatRaster` or a plain `list` of single-layer `SpatRaster`s,
#' which is internally concatenated with **terra**.
#'
#' @param raster_input A multi-layer `SpatRaster` **or** a `list` of single-layer
#'   `SpatRaster` objects representing the same spatial extent/resolution.
#' @param approach Character string specifying the metric to return:
#'   `"difference_count"`, `"shannon_entropy"`, `"stability"`,
#'   `"transition_frequency"`, `"weighted_change_index"`, or `"all"` (default)
#'   for a five-layer stack containing every metric.
#'
#' @return A `SpatRaster`:
#' * **single-layer** if `approach` is one of the named metrics;
#' * **five-layer** (names: `Difference_Count`, `Shannon_Entropy`, `Stability`,
#'   `Transition_Frequency`, `Weighted_Change_Index`) if `approach = "all"`.
#'
#' @seealso \code{\link[terra]{app}}, \code{\link[terra]{rast}}
#'
#' @import terra
#' @import purrr
#' @importFrom stats dist
#'
#' @examples
#' ## -------------------------------------------------------------
#' ## Minimal reproducible example with three random categorical
#' ## rasters (four classes, identical geometry)
#' ## -------------------------------------------------------------
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   set.seed(42)
#'
#'   r1 <- terra::rast(nrows = 40, ncols = 40,
#'                     vals  = sample(1:4, 40 * 40, TRUE))
#'   r2 <- terra::rast(r1, vals = sample(1:4, ncell(r1), TRUE))
#'   r3 <- terra::rast(r1, vals = sample(1:4, ncell(r1), TRUE))
#'
#'   r_stack <- terra::rast(list(r1, r2, r3))
#'   names(r_stack) <- paste0("t", 1:3)
#'
#'   ## 1. All five metrics
#'   diff_all <- map_bioregDiff(r_stack, approach = "all")
#'   print(diff_all)
#'
#'   ## 2. Just the Shannon-entropy layer
#'   ent <- map_bioregDiff(r_stack, approach = "shannon_entropy")
#'   terra::plot(ent, main = "Shannon entropy of label sequence")
#' }
#'
#' @aliases map_bioregDiff
#' @export
map_bioregDiff <- function(raster_input, approach = "all") {
  library(terra)
  library(viridis)
  library(purrr)
  ## if it's a plain list, assume it's a list of SpatRasters and catenate them
  if (is.list(raster_input)) {
    raster_input <- rast(raster_input)
  }

  ## now require a SpatRaster
  if (!inherits(raster_input, "SpatRaster")) {
    stop("Input must be a SpatRaster object (or a list of SpatRaster layers).")
  }

  # Get the number of layers
  nlyr <- nlyr(raster_input)

  # Ensure there are at least two layers for comparison
  if (nlyr < 2) {
    stop("The SpatRaster must have at least two layers for comparison.")
  }

  # Calculate the dissimilarity matrix dynamically
  unique_clusters <- sort(unique(values(raster_input)))
  n_clusters <- length(unique_clusters)

  # Initialize the dissimilarity matrix
  dissimilarity_matrix <- matrix(0, nrow = n_clusters, ncol = n_clusters,
                                 dimnames = list(as.character(unique_clusters), as.character(unique_clusters)))

  # Calculate pairwise transitions across layers
  for (i in 1:(nlyr - 1)) {
    layer1 <- values(raster_input[[i]])
    layer2 <- values(raster_input[[i + 1]])

    transitions <- table(layer1, layer2)

    for (row in 1:nrow(transitions)) {
      for (col in 1:ncol(transitions)) {
        dissimilarity_matrix[rownames(transitions)[row], colnames(transitions)[col]] <-
          dissimilarity_matrix[rownames(transitions)[row], colnames(transitions)[col]] +
          transitions[row, col]
      }
    }
  }

  # Normalize the dissimilarity matrix
  max_transition <- max(dissimilarity_matrix)
  dissimilarity_matrix <- max_transition - dissimilarity_matrix
  dissimilarity_matrix <- (dissimilarity_matrix + t(dissimilarity_matrix)) / 2
  dissimilarity_matrix <- dissimilarity_matrix / max(dissimilarity_matrix)

  # Function: Difference Count
  diff_count <- app(raster_input, fun = function(x) sum(x != x[1]))

  # Function: Shannon Entropy
  shannon_entropy <- function(values) {
    p <- table(values) / length(values)
    p <- p[p > 0]
    -sum(p * log(p))
  }
  entropy_map <- app(raster_input, fun = shannon_entropy)

  # Function: Stability Map
  stable_map <- app(raster_input, fun = function(x) all(x == x[1]))
  stability_map <- 1 - stable_map

  # Function: Transition Frequency
  transition_maps <- list()
  for (i in 1:(nlyr - 1)) {
    transition_maps[[i]] <- raster_input[[i]] != raster_input[[i + 1]]
  }
  total_transitions <- reduce(transition_maps, `+`)

  # Function: Weighted Change Index
  weighted_change_index <- function(values) {
    wci <- 0
    for (i in 1:(length(values) - 1)) {
      cluster1 <- as.character(values[i])
      cluster2 <- as.character(values[i + 1])
      if (cluster1 %in% rownames(dissimilarity_matrix) && cluster2 %in% colnames(dissimilarity_matrix)) {
        weight <- dissimilarity_matrix[cluster1, cluster2]
      } else {
        weight <- 0
      }
      wci <- wci + weight
    }
    return(wci)
  }
  wci_map <- app(raster_input, fun = weighted_change_index)

  # Return selected approach or all
  if (approach == "difference_count") return(diff_count)
  if (approach == "shannon_entropy") return(entropy_map)
  if (approach == "stability") return(stability_map)
  if (approach == "transition_frequency") return(total_transitions)
  if (approach == "weighted_change_index") return(wci_map)

  # Return all as a SpatRaster with multiple layers
  result <- c(diff_count, entropy_map, stability_map, total_transitions, wci_map)
  names(result) <- c("Difference_Count", "Shannon_Entropy", "Stability", "Transition_Frequency", "Weighted_Change_Index")

  return(result)
}
