#' Map bioregion clusters and generate interpolated spatial assignments
#'
#' Performs environmental clustering (k-means, PAM, hierarchical, GMM), aligns cluster labels
#' to k-means reference, plots comparative cluster maps, and interpolates assignments over a raster
#' grid using nearest-neighbor and/or thin-plate spline methods.
#'
#' @param data      A data.frame containing coordinate columns and values for clustering.
#' @param scale_cols Character vector of column names in `data` to scale before clustering.
#' @param clus_method Character; one or more clustering methods: "kmeans", "pam", "hclust", "gmm", or "all". Default = "all".
#' @param plot      Logical; if TRUE, displays side-by-side cluster maps. Default = TRUE.
#' @param interp    Character; interpolation method(s): "NN" (nearest neighbor), "Tps" (thin-plate spline), or "both". Default = "both".
#' @param x_col     Name of the longitude column in `data`. Default = "x".
#' @param y_col     Name of the latitude  column in `data`. Default = "y".
#'
#' @return A list with elements:
#'   \item{data}{Original data with appended cluster assignment columns (`kmeans_opt`, `pam_opt_aligned`, etc.)}
#'   \item{interpolated}{SpatRaster of interpolated cluster values}
#'   \item{clusters}{SpatRaster of discrete clusters after interpolation}
#'
#' @import NbClust
#' @import clValid
#' @import cluster
#' @import purrr
#' @importFrom factoextra hcut
#' @importFrom stats kmeans dist
#' @import terra
#' @import fields
#' @import sf
#' @import mclust
#' @import patchwork
#' @import ggplot2
#' @export
map_bioreg <- function(data,
                       scale_cols,
                       clus_method = "all",
                       plot        = TRUE,
                       interp      = "both",
                       x_col       = "x",
                       y_col       = "y",
                       resolution = 0.5) {

  #--- Libraries
  library(NbClust); library(clValid); library(cluster)
  library(factoextra); library(ggplot2); library(terra)
  library(fields); library(sf);      library(mclust)
  library(patchwork)

  #--- sanity checks
  if (!all(c(x_col, y_col, scale_cols) %in% names(data))) {
    stop("`data` must contain columns: ",
         paste(c(x_col, y_col, scale_cols), collapse = ", "))
  }

  #--- helper to align cluster IDs by centroid proximity
  align_clusters <- function(df, ref_col, target_col) {
    ref_cntrs    <- aggregate(df[, c(x_col, y_col)],
                              by = list(cluster = df[[ref_col]]),
                              FUN = mean)
    target_cntrs <- aggregate(df[, c(x_col, y_col)],
                              by = list(cluster = df[[target_col]]),
                              FUN = mean)

    dm <- as.matrix(dist(rbind(ref_cntrs[,2:3], target_cntrs[,2:3])))
    dm <- dm[1:nrow(ref_cntrs),
             (nrow(ref_cntrs)+1):nrow(dm)]

    mapping <- apply(dm, 2, which.min)
    vapply(df[[target_col]], function(i) mapping[i], integer(1))
  }

  #--- 1. scale
  data_scaled <- scale(data[, scale_cols])

  #--- 2. optimal k for k-means
  optimal_k <- 6
  if (clus_method %in% c("kmeans","all")) {
    res <- NbClust(data_scaled, distance="euclidean",
                   min.nc=2, max.nc=10,
                   method="kmeans", index="silhouette")
    optimal_k <- res$Best.nc[1]
    message("Optimal k (silhouette): ", optimal_k)
  }

  #--- 3. clustering
  clusters <- list()
  if (clus_method %in% c("kmeans","all")) {
    km <- kmeans(data_scaled, centers=optimal_k, nstart=10)
    data$kmeans_opt <- km$cluster; clusters$kmeans <- km
  }
  if (clus_method %in% c("pam","all")) {
    pm <- pam(data_scaled, k=optimal_k)
    data$pam_opt <- pm$clustering; clusters$pam <- pm
    data$pam_opt_aligned <- align_clusters(data, "kmeans_opt", "pam_opt")
  }
  if (clus_method %in% c("hclust","all")) {
    hc <- hcut(data_scaled, k=optimal_k, hc_method="complete")
    data$hc_opt <- hc$cluster; clusters$hclust <- hc
    data$hc_opt_aligned <- align_clusters(data, "kmeans_opt", "hc_opt")
  }
  if (clus_method %in% c("gmm","all")) {
    gm <- Mclust(data_scaled, G=optimal_k)
    data$gmm_opt <- gm$classification; clusters$gmm <- gm
    data$gmm_opt_aligned <- align_clusters(data, "kmeans_opt", "gmm_opt")
  }

  #--- 4. plotting
  if (plot) {
    plt <- list()
    mkr <- function(col, title) {
      ggplot(data) +
        geom_tile(aes_string(x = x_col, y = y_col, fill = col)) +
        scale_fill_brewer(palette="Set1") +
        labs(title=title,
             x="Longitude", y="Latitude",
             fill="Cluster") +
        theme_minimal()
    }
    if ("kmeans_opt" %in% names(data))
      plt$kmeans <- mkr("factor(kmeans_opt)", "K-means")
    if ("pam_opt_aligned" %in% names(data))
      plt$pam    <- mkr("factor(pam_opt_aligned)", "PAM")
    if ("hc_opt_aligned" %in% names(data))
      plt$hclust <- mkr("factor(hc_opt_aligned)", "Hierarchical")
    if ("gmm_opt_aligned" %in% names(data))
      plt$gmm    <- mkr("factor(gmm_opt_aligned)", "GMM")
    print((plt$kmeans + plt$pam)/(plt$hclust + plt$gmm))
  }

  #--- 5. interpolation
  # Step 5: Interpolation
  interpolated <- NULL
  clusters_rcl <- NULL

  if (interp %in% c("NN", "both")) {
    # build an sf points object & get its bounding box
    pts_sf <- sf::st_as_sf(data, coords = c(x_col, y_col), crs = 4326)
    bb     <- sf::st_bbox(pts_sf)

    # create a Terra extent and empty raster
    ex     <- terra::ext(bb$xmin, bb$xmax, bb$ymin, bb$ymax)
    r_base <- terra::rast(ex, resolution = resolution, crs = "EPSG:4326")

    # pre-compute site coords
    coords <- sf::st_coordinates(pts_sf)

    # nearest-neighbor assignment for each clustering
    for (m in c("kmeans_opt","pam_opt_aligned","hc_opt_aligned","gmm_opt_aligned")) {
      nn <- apply(
        as.matrix(terra::crds(r_base)),
        1,
        function(pt) {
          i <- which.min((coords[,1] - pt[1])^2 + (coords[,2] - pt[2])^2)
          data[[m]][i]
        }
      )
      # add a layer named after m
      r_base[[m]] <- nn
    }

    interpolated <- r_base
  }

  if (interp %in% c("Tps", "both")) {
    # build a Terra vector and extent
    pts_vec <- terra::vect(data, geom = c(x_col, y_col), crs = "EPSG:4326")
    ex2     <- terra::ext(pts_vec)
    grid    <- terra::rast(ex2, resolution = resolution, crs = terra::crs(pts_vec))

    # thin‐plate spline models
    tps_models <- list(
      kmeans_tps = fields::Tps(data[,c(x_col,y_col)], data$kmeans_opt),
      pam_tps    = fields::Tps(data[,c(x_col,y_col)], data$pam_opt_aligned),
      hc_tps     = fields::Tps(data[,c(x_col,y_col)], data$hc_opt_aligned),
      gmm_tps    = fields::Tps(data[,c(x_col,y_col)], data$gmm_opt_aligned)
    )

    # interpolate each one over the grid
    preds <- lapply(tps_models, function(mod) terra::interpolate(grid, mod))
    stack <- do.call(c, preds)
    names(stack) <- names(preds)

    # example reclassification (adjust breaks as needed)
    # matrix: from, to, becomes
    rcl_mat <- matrix(c(
      0, 1.5, 1,
      1.5, 2.5, 2,
      2.5, 3.5, 3,
      3.5, 4.5, 4,
      4.5, 5.5, 5,
      5.5, 6.5, 6
    ), ncol = 3, byrow = TRUE)

    clusters_rcl <- lapply(stack, classify, rcl = rcl_mat)
    clusters_rcl <- do.call(c, clusters_rcl)
  }

  return(list(
    data         = data,
    interpolated = interpolated,
    clusters     = clusters_rcl
  ))
}
