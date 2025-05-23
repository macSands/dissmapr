# library(terra)
# library(viridis)
# library(purrr)

# Function to map bioregional differences
map_bioregDiff <- function(raster_input, approach = "all") {
  # Check if input is a SpatRaster
  if (!inherits(raster_input, "SpatRaster")) {
    stop("Input must be a SpatRaster object.")
  }

  library(terra)
  library(viridis)
  library(purrr)

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

# # Example Usage
# # Assuming 'bioreg_result$clusters' is a SpatRaster
# result_bioregDiff <- map_bioregDiff(bioreg_result$clusters, approach = "all")
# result_bioregDiff = mask(result_bioregDiff, masked_grid)
# # Plot all layers
# plot(result_bioregDiff, col = viridis(100, direction = -1))
