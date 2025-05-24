#' Map bioregion clusters and interpolate spatial assignments
#'
#' Performs clustering on scaled variables to delineate bioregions, aligns cluster labels
#' across methods, generates comparative cluster maps, and interpolates cluster assignments
#' over a raster grid using nearest-neighbor (NN) and/or thin-plate spline (Tps) methods.
#'
#' @param data A `data.frame` containing at least coordinate columns (`x`, `y`) and the
#'   variables specified in `scale_cols` to be used for clustering.
#' @param scale_cols Character vector of column names in `data` to scale before clustering.
#' @param clus_method Character; clustering method(s) to apply. One of:
#'   \describe{
#'     \item{"kmeans"}{k-means clustering}
#'     \item{"pam"}{Partitioning Around Medoids}
#'     \item{"hclust"}{Hierarchical clustering (complete linkage)}
#'     \item{"gmm"}{Gaussian Mixture Model clustering}
#'     \item{"all"}{All of the above methods}
#'   }
#'   Default is "all".
#' @param plot Logical; if `TRUE`, produces side-by-side cluster maps for each method. Default is `TRUE`.
#' @param interp Character; interpolation method for raster assignment. One of:
#'   \describe{
#'     \item{"NN"}{nearest-neighbor interpolation}
#'     \item{"Tps"}{thin-plate spline interpolation}
#'     \item{"both"}{both NN and Tps}
#'   }
#'   Default is "both".
#'
#' @return A list with elements:
#'   \item{data}{Original data with appended cluster assignment columns (`kmeans_opt`, `pam_opt_aligned`, etc.)}
#'   \item{interpolated}{A `SpatRaster` or stack of rasters of interpolated cluster assignments}
#'   \item{clusters}{A `SpatRaster` of discretized cluster classes after interpolation}
#'
#' @details
#' 1. Scales the columns in `scale_cols` using `scale()`.
#' 2. If `clus_method` includes "kmeans", automatically determines optimal k via
#'    `NbClust::NbClust()` (silhouette index, k = 2–10).
#' 3. Applies each requested clustering method, aligns non-kmeans labels to kmeans reference.
#' 4. If `plot = TRUE`, generates comparative cluster maps using `ggplot2` and `patchwork`.
#' 5. If `interp` includes "NN", assigns each raster cell the cluster of its nearest point.
#' 6. If `interp` includes "Tps", fits thin-plate spline surfaces (via `fields::Tps()`) for each
#'    cluster layer and predicts over the raster grid.
#'
#' @importFrom NbClust NbClust
#' @importFrom cluster pam
#' @importFrom factoextra hcut
#' @importFrom mclust Mclust
#' @importFrom stats kmeans dist aggregate
#' @importFrom terra rast mask classify crds vect interpolate
#' @importFrom sf st_as_sf st_coordinates
#' @importFrom ggplot2 ggplot geom_tile scale_fill_brewer labs theme_minimal
#' @importFrom patchwork wrap_plots
#' @importFrom fields Tps
#' @export
map_bioreg <- function(data,
                       scale_cols,
                       clus_method = "all",
                       plot        = TRUE,
                       interp      = "both") {
  # Ensure required columns
  if (!all(c("x", "y") %in% names(data))) {
    stop("`data` must contain columns `x` and `y`")
  }

  # 1. Scale data
  data_scaled <- scale(data[, scale_cols, drop = FALSE])

  # 2. Determine optimal k for kmeans
  optimal_k <- 6L
  if (clus_method %in% c("kmeans", "all")) {
    nc <- NbClust::NbClust(
      data    = data_scaled,
      distance= "euclidean",
      min.nc  = 2,
      max.nc  = 10,
      method  = "kmeans",
      index   = "silhouette"
    )
    optimal_k <- nc$Best.nc[1]
    message("Optimal number of clusters identified: ", optimal_k)
  }

  # Prepare storage
  clusters <- list()

  # 3a. kmeans
  if (clus_method %in% c("kmeans", "all")) {
    km_res <- stats::kmeans(data_scaled, centers = optimal_k, nstart = 10)
    data$kmeans_opt <- km_res$cluster
    clusters$kmeans <- km_res
  }

  # 3b. PAM
  if (clus_method %in% c("pam", "all")) {
    pam_res <- cluster::pam(data_scaled, k = optimal_k)
    data$pam_opt <- pam_res$clustering
    clusters$pam <- pam_res
    data$pam_opt_aligned <- align_clusters(data, "kmeans_opt", "pam_opt")
  }

  # 3c. Hierarchical clustering
  if (clus_method %in% c("hclust", "all")) {
    hc_cut <- factoextra::hcut(data_scaled, k = optimal_k, hc_method = "complete")
    data$hc_opt <- hc_cut$cluster
    clusters$hclust <- hc_cut
    data$hc_opt_aligned <- align_clusters(data, "kmeans_opt", "hc_opt")
  }

  # 3d. Gaussian Mixture Model
  if (clus_method %in% c("gmm", "all")) {
    gmm_res <- mclust::Mclust(data_scaled, G = optimal_k)
    data$gmm_opt <- gmm_res$classification
    clusters$gmm <- gmm_res
    data$gmm_opt_aligned <- align_clusters(data, "kmeans_opt", "gmm_opt")
  }

  # 4. Plot clusters
  if (plot) {
    plots <- list()
    if (!is.null(data$kmeans_opt)) {
      plots$kmeans <-
        ggplot2::ggplot(data) +
        ggplot2::geom_tile(
          aes(x = x, y = y, fill = factor(kmeans_opt))
        ) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::labs(
          title = "K-means clustering",
          x     = "Longitude",
          y     = "Latitude",
          fill  = "Cluster"
        ) +
        ggplot2::theme_minimal()
    }
    if (!is.null(data$pam_opt_aligned)) {
      plots$pam <-
        ggplot2::ggplot(data) +
        ggplot2::geom_tile(
          aes(x = x, y = y, fill = factor(pam_opt_aligned))
        ) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::labs(
          title = "PAM clustering",
          x     = "Longitude",
          y     = "Latitude",
          fill  = "Cluster"
        ) +
        ggplot2::theme_minimal()
    }
    if (!is.null(data$hc_opt_aligned)) {
      plots$hclust <-
        ggplot2::ggplot(data) +
        ggplot2::geom_tile(
          aes(x = x, y = y, fill = factor(hc_opt_aligned))
        ) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::labs(
          title = "Hierarchical clustering",
          x     = "Longitude",
          y     = "Latitude",
          fill  = "Cluster"
        ) +
        ggplot2::theme_minimal()
    }
    if (!is.null(data$gmm_opt_aligned)) {
      plots$gmm <-
        ggplot2::ggplot(data) +
        ggplot2::geom_tile(
          aes(x = x, y = y, fill = factor(gmm_opt_aligned))
        ) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::labs(
          title = "GMM clustering",
          x     = "Longitude",
          y     = "Latitude",
          fill  = "Cluster"
        ) +
        ggplot2::theme_minimal()
    }
    # combine in 2x2 grid
    print(patchwork::wrap_plots(plots, ncol = 2))
  }

  # 5. Interpolation
  # nearest-neighbor
  interp_rasters <- list()
  if (interp %in% c("NN", "both")) {
    pts <- sf::st_as_sf(data, coords = c("x", "y"), crs = 4326)
    r   <- terra::rast(terra::vect(pts), resolution = 0.05)
    coords_mat <- sf::st_coordinates(pts)
    nn_fun <- function(cluster_col) {
      apply(terra::crds(r), 1, function(coord) {
        d <- sqrt((coords_mat[,1] - coord[1])^2 + (coords_mat[,2] - coord[2])^2)
        data[[cluster_col]][which.min(d)]
      })
    }
    nn_rasters <- sapply(c("kmeans_opt", "pam_opt_aligned",
                           "hc_opt_aligned", "gmm_opt_aligned"), function(col) {
                             if (!is.null(data[[col]])) {
                               terra::classify(
                                 terra::setValues(r, nn_fun(col)),
                                 rcl = matrix(c(-Inf, Inf, NA), ncol = 3)
                               )
                             }
                           }, simplify = FALSE)
    interp_rasters$NN <- terra::rast(nn_rasters)
  }
  # thin-plate spline
  if (interp %in% c("Tps", "both")) {
    grid <- terra::rast(terra::vect(pts), resolution = 0.25)
    masked <- terra::mask(grid, terra::vect(pts))
    tps_fun <- function(cluster_col) {
      fields::Tps(sf::st_coordinates(pts), data[[cluster_col]])
    }
    tps_models <- lapply(c("kmeans_opt", "pam_opt_aligned",
                           "hc_opt_aligned", "gmm_opt_aligned"), function(col) {
                             if (!is.null(data[[col]])) tps_fun(col)
                           })
    tps_preds <- mapply(function(model, col) {
      terra::interpolate(masked, model)
    }, model = tps_models,
    SIMPLIFY = FALSE)
    interp_rasters$Tps <- terra::rast(tps_preds)
  }

  # Stack and classify
  all_interp <- do.call(c, interp_rasters)
  discrete <- terra::classify(all_interp, rcl = matrix(
    c(0, 1.5, 1,
      1.5, 2.5, 2,
      2.5, 3.5, 3,
      3.5, 4.5, 4,
      4.5, 5.5, 5,
      5.5, 6.5, 6), byrow = TRUE, ncol = 3
  ))

  return(list(
    data         = data,
    interpolated = all_interp,
    clusters     = discrete
  ))
}
