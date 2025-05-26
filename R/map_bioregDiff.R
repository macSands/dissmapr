#' Map bioregional difference metrics across raster layers
#'
#' Takes two or more categorical raster layers (or a list of them) and computes how
#' each cell's region label changes over time.  You can return:
#' - the **difference count** (number of times a cell's label changes),
#' - **Shannon entropy** of its label sequence,
#' - **stability** (proportion of time unchanged),
#' - **transition frequency** (sum of per-step changes),
#' - **weighted change index** (sum of pairwise-change weights),
#' or all five as a multi‐layer `SpatRaster`.
#'
#' @aliases map_bioregDiff
#' @param raster_input A `SpatRaster` or a **list** of single‐layer `SpatRaster` objects.
#' @param approach Character; one of
#'   \describe{
#'     \item{`"difference_count"`}{Count of label changes per cell.}
#'     \item{`"shannon_entropy"`}{Entropy of the cell's label distribution.}
#'     \item{`"stability"`}{Proportion of layers where the label stays the same.}
#'     \item{`"transition_frequency"`}{Sum of binary change maps across layers.}
#'     \item{`"weighted_change_index"`}{Sum of dissimilarity‐weighted transitions.}
#'     \item{`"all"`}{All five metrics combined into a multi‐layer `SpatRaster`.}
#'   }
#'   Default: `"all"`.
#' @return A `SpatRaster`:
#' - A single layer if `approach` is one metric,
#' - A 5‐layer raster (`Difference_Count`, `Shannon_Entropy`, `Stability`,
#'   `Transition_Frequency`, `Weighted_Change_Index`) if `approach = "all"`.
#'
#' @examples
#' \dontrun{
#' # Combine four single‐layer rasters into one multi‐layer object:
#' rlist <- list(km=km_rast, pam=pam_rast, hc=hc_rast, gmm=gmm_rast)
#' multi <- rast(rlist)
#'
#' # Compute all metrics:
#' diff_metrics <- map_bioregDiff(multi, approach="all")
#'
#' # Just the entropy map:
#' ent_map <- map_bioregDiff(multi, approach="shannon_entropy")
#' }
#'
#' @import terra
#' @importFrom purrr reduce
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
