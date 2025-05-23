map_bioreg <- function(data,
                       scale_cols,
                       clus_method = "all",
                       plot = TRUE,
                       interp = "both") {
  library(NbClust)
  library(clValid)
  library(cluster)
  library(factoextra) # Required for hcut
  library(ggplot2)
  library(terra)
  library(fields)
  library(sf)
  library(mclust)  # Required for GMM
  library(patchwork)

  align_clusters <- function(data, ref_col, target_col) {
    # Calculate centroids for reference and target clustering
    ref_centroids <- aggregate(data[, c("x", "y")], by = list(data[[ref_col]]), FUN = mean)
    target_centroids <- aggregate(data[, c("x", "y")], by = list(data[[target_col]]), FUN = mean)

    # Compute pairwise distances between centroids
    dist_matrix <- as.matrix(dist(rbind(ref_centroids[, 2:3], target_centroids[, 2:3])))
    dist_matrix <- dist_matrix[1:nrow(ref_centroids), (nrow(ref_centroids) + 1):nrow(dist_matrix)]

    # Find the closest clusters
    mapping <- apply(dist_matrix, 2, which.min)

    # Reassign target cluster IDs
    new_ids <- sapply(data[[target_col]], function(x) mapping[x])
    return(new_ids)
  }

  # Step 1: Scale specified columns
  data_scaled <- scale(data[, scale_cols])

  # Step 2: Determine optimal clusters for k-means
  optimal_k <- 6 # Default value if not calculated
  if (clus_method %in% c("kmeans", "all")) {
    res <- NbClust(data_scaled,
                   distance = "euclidean",
                   min.nc = 2,
                   max.nc = 10,
                   method = "kmeans",
                   index = "silhouette")
    optimal_k_val <- res$Best.nc[1]
    print(paste0('Optimal number of clusters identified: ', optimal_k_val))
    optimal_k = 5
  }

  # Step 3: Perform clustering
  clusters <- list()

  if (clus_method %in% c("kmeans", "all")) {
    km.res <- kmeans(data_scaled, centers = optimal_k, nstart = 10)
    data$kmeans_opt <- km.res$cluster
    clusters$kmeans <- km.res
  }

  # if (clus_method %in% c("pam", "all")) {
  #   pam.res <- pam(data_scaled, k = optimal_k)
  #   data$pam_opt <- pam.res$clustering
  #   clusters$pam <- pam.res
  # }
  if (clus_method %in% c("pam", "all")) {
    pam.res <- pam(data_scaled, k = optimal_k)
    data$pam_opt <- pam.res$clustering
    clusters$pam <- pam.res
    data$pam_opt_aligned <- align_clusters(data, "kmeans_opt", "pam_opt")
  }

  # if (clus_method %in% c("hclust", "all")) {
  #   hc.cut <- hcut(data_scaled, k = optimal_k, hc_method = "complete")
  #   data$hc_opt <- hc.cut$cluster
  #   clusters$hclust <- hc.cut
  # }
  if (clus_method %in% c("hclust", "all")) {
    hc.cut <- hcut(data_scaled, k = optimal_k, hc_method = "complete")
    data$hc_opt <- hc.cut$cluster
    clusters$hclust <- hc.cut
    data$hc_opt_aligned <- align_clusters(data, "kmeans_opt", "hc_opt")
  }

  # if (clus_method %in% c("gmm", "all")) {
  #   # Perform GMM (Gaussian Mixture Models)
  #   gmm.res <- Mclust(data_scaled, G = optimal_k)
  #   data$gmm_opt <- gmm.res$classification
  #   clusters$gmm <- gmm.res
  # }
  if (clus_method %in% c("gmm", "all")) {
    gmm.res <- Mclust(data_scaled, G = optimal_k)
    data$gmm_opt <- gmm.res$classification
    clusters$gmm <- gmm.res
    data$gmm_opt_aligned <- align_clusters(data, "kmeans_opt", "gmm_opt")
  }

  # Step 4: Plot clusters
  if (plot) {
    plot_list <- list()

    if ("kmeans_opt" %in% names(data)) {
      plot_list$kmeans <- ggplot(data) +
        geom_tile(aes(x = x, y = y, fill = factor(kmeans_opt))) +
        scale_fill_brewer(palette = "Set1") +
        labs(title = "K-means clustering", x = "Longitude", y = "Latitude", fill = "Cluster") +
        theme_minimal()
    }

    if ("pam_opt_aligned" %in% names(data)) {
      plot_list$pam <- ggplot(data) +
        geom_tile(aes(x = x, y = y, fill = factor(pam_opt_aligned))) +
        scale_fill_brewer(palette = "Set1") +
        labs(title = "PAM clustering", x = "Longitude", y = "Latitude", fill = "Cluster") +
        theme_minimal()
    }

    if ("hc_opt_aligned" %in% names(data)) {
      plot_list$hclust <- ggplot(data) +
        geom_tile(aes(x = x, y = y, fill = factor(hc_opt_aligned))) +
        scale_fill_brewer(palette = "Set1") +
        labs(title = "Hierarchical clustering", x = "Longitude", y = "Latitude", fill = "Cluster") +
        theme_minimal()
    }

    # if ("gmm_opt" %in% names(data)) {
    #   plot_list$gmm <- ggplot(data) +
    #     geom_tile(aes(x = x, y = y, fill = factor(gmm_opt))) +
    #     scale_fill_brewer(palette = "Set1") +
    #     labs(title = "GMM clustering", x = "Longitude", y = "Latitude", fill = "Cluster") +
    #     theme_minimal()
    # }

    if ("gmm_opt_aligned" %in% names(data)) {
      plot_list$gmm <- ggplot(data) +
        geom_tile(aes(x = x, y = y, fill = factor(gmm_opt_aligned))) +
        scale_fill_brewer(palette = "Set1") +
        labs(title = "GMM Clustering", x = "Longitude", y = "Latitude", fill = "Cluster") +
        theme_minimal()
    }

    # Combine plots into a grid (2x2)
    combined_plot <- (plot_list$kmeans + plot_list$pam) /
      (plot_list$hclust + plot_list$gmm)
    print(combined_plot)
  }

  # Step 5: Interpolation (unchanged)
  # Interpolation code here (use the original function’s implementation)
  # Step 5: Interpolate missing values
  if (interp %in% c("NN", "both")) {
    # rsa <- st_read('D:/Data/South_Africa/southafrica_provinces_lesotho_swaziland.shp')
    pts <- st_as_sf(data, coords = c("x", "y"), crs = 4326)
    r <- rast(rsa, resolution = 0.05)
    coords <- st_coordinates(pts)

    # For each pixel, find nearest point (and thus cluster)
    # This can be done using, e.g., a KD-tree or similar approach:
    nn_cluster_kmeans = apply(as.matrix(crds(r)), 1, function(coord) {
      dists = sqrt((coords[,1]-coord[1])^2 + (coords[,2]-coord[2])^2)
      # pts$cluster_4[which.min(dists)]
      pts$kmeans_opt[which.min(dists)]
    })

    nn_cluster_pam = apply(as.matrix(crds(r)), 1, function(coord) {
      dists = sqrt((coords[,1]-coord[1])^2 + (coords[,2]-coord[2])^2)
      # pts$cluster_4[which.min(dists)]
      pts$pam_opt_aligned[which.min(dists)]
    })

    nn_cluster_hcut = apply(as.matrix(crds(r)), 1, function(coord) {
      dists = sqrt((coords[,1]-coord[1])^2 + (coords[,2]-coord[2])^2)
      # pts$cluster_4[which.min(dists)]
      pts$hc_opt_aligned[which.min(dists)]
    })

    nn_cluster_gmm = apply(as.matrix(crds(r)), 1, function(coord) {
      dists = sqrt((coords[,1]-coord[1])^2 + (coords[,2]-coord[2])^2)
      # pts$cluster_4[which.min(dists)]
      pts$gmm_opt_aligned[which.min(dists)]
    })

    # values(r) = nn_cluster_kmeans
    r$kmeans_nn = nn_cluster_kmeans
    r$pam_nn = nn_cluster_pam
    r$hcut_nn = nn_cluster_hcut
    r$gmm_nn = nn_cluster_gmm

    # Mask by SA boundary
    r = mask(r, vect(rsa))

    # Plot raster
    plot(r)
  }

  if (interp %in% c("Tps", "both")) {
    # rsa <- st_read('D:/Data/South_Africa/southafrica_provinces_lesotho_swaziland.shp')
    vect_pts <- vect(data, crs = "EPSG:4326", geom = c("x", "y"))
    grid <- rast(rsa, resolution = 0.25, crs = "EPSG:4326")
    values(grid) = 1  # Assign an initial value (e.g., 1) to all cells
    # Convert the South Africa boundary to a SpatVector
    vect_rsa = vect(rsa)
    # Mask the raster grid using the boundary of South Africa
    masked_grid = mask(grid, vect_rsa)
    tps_kmeans = Tps(data[,c('x','y')], data$kmeans_opt)
    tps_pam = Tps(data[,c('x','y')], data$pam_opt_aligned)
    tps_hc = Tps(data[,c('x','y')], data$hc_opt_aligned)
    tps_gmm = Tps(data[,c('x','y')], data$gmm_opt_aligned)
    p = rast(masked_grid)

    # use model to predict values at all locations
    p_kmeans = interpolate(p, tps_kmeans)
    p_pam = interpolate(p, tps_pam)
    p_hc = interpolate(p, tps_hc)
    p_gmm = interpolate(p, tps_gmm)
    p_clusters = c(p_kmeans,p_pam,p_hc,p_gmm)
    p_clusters = mask(p_clusters, masked_grid)
    names(p_clusters) = c('kmeans_tps','pam_tps','hcut_tps','gmm_tps')
    plot(p_clusters)

    # Define a reclassification matrix
    reclass_matrix = matrix(
      c(0, 1.5, 1,   # Values between 0 and 1.5 -> 1
        1.5, 2.5, 2, # Values between 1.5 and 2.5 -> 2
        2.5, 3.5, 3, # Values between 2.5 and 3.5 -> 3
        3.5, 4.5, 4,
        4.5, 5.5, 5,
        5.5, 6.5, 6),  # Values between 3.5 and 5 -> 4
      ncol = 3,
      byrow = TRUE
    )
    # Apply the classification
    p_kmeans_discrete = classify(p_clusters[[1]], rcl = reclass_matrix)
    p_pam_discrete = classify(p_clusters[[2]], rcl = reclass_matrix)
    p_hc_discrete = classify(p_clusters[[3]], rcl = reclass_matrix)
    p_gmm_discrete = classify(p_clusters[[4]], rcl = reclass_matrix)
    p_all_discrete = c(p_kmeans_discrete,p_pam_discrete,p_hc_discrete,p_gmm_discrete)
    names(p_all_discrete) = c('kmeans_rcl','pam_rcl','hcut_rcl','gmm_rcl')
    plot(p_all_discrete)
  }

  return(list(data=data, interpolated = p_clusters, clusters = p_all_discrete))
}

# # USAGE RESULTS
# # Example dataset
# bioreg_result <- map_bioreg(
#   data = predictors_df,
#   scale_cols = c("pred_zetaExp", "x", "y"),
#   clus_method = "all", # Includes DBSCAN
#   plot = TRUE,
#   interp = "both"
# )
#
# # Inspect results
# str(bioreg_result, max.level = 1)
# head(bioreg_result$data)
# bioreg_result$interpolated
# bioreg_result$clusters
