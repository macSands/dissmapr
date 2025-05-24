#' Map bioregional difference metrics across raster layers
#'
#' Computes a suite of change indices—including difference count, Shannon entropy,
#' stability, transition frequency, and weighted change index—across sequential layers
#' of a `SpatRaster`. Returns either a single selected metric or all metrics as a
#' multi-layer `SpatRaster`.
#'
#' @param raster_input A `SpatRaster` with at least two layers representing categorical
#'   cluster assignments (e.g., bioregions) over time or scenarios.
#' @param approach Character specifying which metric to return. Options:
#'   \describe{
#'     \item{"difference_count"}{Count of pixel value changes across layers.}
#'     \item{"shannon_entropy"}{Shannon entropy of pixel value distribution.}
#'     \item{"stability"}{Proportion of layers where pixel value remains constant.}
#'     \item{"transition_frequency"}{Number of transitions (binary change) per pixel.}
#'     \item{"weighted_change_index"}{Sum of transition weights from the normalized dissimilarity matrix.}
#'     \item{"all"}{All metrics combined into a multi-layer `SpatRaster`.}
#'   }
#'   Default is "all".
#'
#' @return A `SpatRaster`: a single-layer raster if `approach` is specific;
#'   or a multi-layer raster with layers:
#'   `Difference_Count`, `Shannon_Entropy`, `Stability`,
#'   `Transition_Frequency`, `Weighted_Change_Index` when `approach = "all"`.
#'
#' @details
#' 1. Builds a transition matrix of cluster changes between consecutive layers and
#'    normalizes it to create a dissimilarity matrix.
#' 2. `difference_count`: counts how many times each pixel changes value.
#' 3. `shannon_entropy`: computes entropy of pixel value distribution across layers.
#' 4. `stability`: calculates 1 minus the indicator of constant pixel values.
#' 5. `transition_frequency`: sums per-layer binary change maps.
#' 6. `weighted_change_index`: accumulates weighted transitions using the dissimilarity matrix.
#'
#' @importFrom terra nlyr values app
#' @importFrom stats dist
#' @importFrom purrr reduce
#' @export
map_bioregDiff <- function(raster_input,
                           approach = "all") {
  if (!inherits(raster_input, "SpatRaster")) {
    stop("Input must be a SpatRaster object.")
  }

  nlyr_val <- terra::nlyr(raster_input)
  if (nlyr_val < 2) {
    stop("The SpatRaster must have at least two layers for comparison.")
  }

  # Build transition dissimilarity matrix
  vals <- terra::values(raster_input)
  unique_clusters <- sort(unique(vals))
  dissim <- matrix(0,
                   nrow = length(unique_clusters),
                   ncol = length(unique_clusters),
                   dimnames = list(unique_clusters, unique_clusters))

  for (i in seq_len(nlyr_val - 1)) {
    layer1 <- vals[, i]
    layer2 <- vals[, i + 1]
    trans <- stats::table(layer1, layer2)
    for (r in seq_len(nrow(trans))) {
      for (c in seq_len(ncol(trans))) {
        dissim[rownames(trans)[r], colnames(trans)[c]] <-
          dissim[rownames(trans)[r], colnames(trans)[c]] + trans[r, c]
      }
    }
  }
  max_t <- max(dissim)
  dissim <- (max_t - dissim + t(dissim)) / 2
  dissim <- dissim / max(dissim)

  # difference count
  diff_count <- terra::app(raster_input, fun = function(x) sum(x != x[1]))

  # Shannon entropy
  entropy_map <- terra::app(raster_input, fun = function(x) {
    p <- stats::table(x) / length(x)
    p <- p[p > 0]
    -sum(p * log(p))
  })

  # stability
  stable_map <- terra::app(raster_input, fun = function(x) all(x == x[1]))
  stability_map <- 1 - stable_map

  # transition frequency
  transition_layers <- lapply(
    seq_len(nlyr_val - 1),
    function(i) raster_input[[i]] != raster_input[[i + 1]]
  )
  total_transitions <- purrr::reduce(transition_layers, `+`)

  # weighted change index
  wci_map <- terra::app(raster_input, fun = function(x) {
    wci <- 0
    for (j in seq_len(length(x) - 1)) {
      c1 <- as.character(x[j]); c2 <- as.character(x[j + 1])
      wt <- if (c1 %in% rownames(dissim) && c2 %in% colnames(dissim))
        dissim[c1, c2] else 0
      wci <- wci + wt
    }
    wci
  })

  # return selected
  if (approach == "difference_count") return(diff_count)
  if (approach == "shannon_entropy")  return(entropy_map)
  if (approach == "stability")        return(stability_map)
  if (approach == "transition_frequency") return(total_transitions)
  if (approach == "weighted_change_index") return(wci_map)

  result <- terra::c(
    diff_count, entropy_map, stability_map,
    total_transitions, wci_map
  )
  names(result) <- c(
    "Difference_Count", "Shannon_Entropy",
    "Stability", "Transition_Frequency",
    "Weighted_Change_Index"
  )
  result
}
