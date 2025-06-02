#' Map Bioregion Clusters and Spatially Interpolate Assignments
#'
#' @description
#' Runs several unsupervised clustering algorithms on multivariate
#' environmental data, *aligns* their cluster labels to a common
#' k-means reference, visualises the resulting partitions, and produces
#' gridded bioregion surfaces by nearest-neighbour (NN) and/or
#' thin-plate-spline (TPS) interpolation.
#'
#' @details
#' * **Clustering**   The optimal number of clusters *k* for k-means is
#'   selected with **NbClust** using the silhouette criterion
#'   (`min.nc = 2`, `max.nc = 10`).
#' * **Label alignment**   Because different algorithms label clusters
#'   arbitrarily, labels are remapped to match k-means centroids via a
#'   nearest-centroid assignment.
#' * **Interpolation**   For each algorithm the point-wise clusters are
#'   projected to a raster grid:
#'   – **NN**: a computationally cheap 1-nearest neighbour assignment
#'   – **TPS**: smooth thin-plate splines fitted with **fields**.
#'   Predicted fractional surfaces are optionally re-classified to
#'   discrete cluster IDs.
#'
#' @param data A data frame containing longitude / latitude columns and the
#'   variables in `scale_cols`.
#' @param scale_cols Character vector naming columns to *z*-scale before
#'   clustering.
#' @param clus_method Character; one or more of `"kmeans"`, `"pam"`,
#'   `"hclust"`, `"gmm"`, or `"all"` (default) for all four.
#' @param plot Logical; if `TRUE` (default) a 2×2 patchwork of cluster
#'   maps is printed.
#' @param interp Character; interpolation type: `"NN"`, `"Tps"`, or
#'   `"both"` (default).
#' @param x_col,y_col Names of the longitude and latitude columns
#'   (defaults `"x"`, `"y"`).
#' @param resolution Numeric. Grid cell size (map units) used to build the
#'   interpolation raster (default `0.5`).
#'
#' @return A list with three components
#' \describe{
#'   \item{data}{Original `data` with added cluster columns
#'     (`kmeans_opt`, `pam_opt_aligned`, `hc_opt_aligned`,
#'     `gmm_opt_aligned`).}
#'   \item{interpolated}{A `SpatRaster` stack of continuous NN or TPS
#'     predictions (one layer per algorithm).}
#'   \item{clusters}{A `SpatRaster` stack of re-classified discrete
#'     clusters (TPS only, `NULL` otherwise).}
#' }
#'
#' @importFrom NbClust NbClust
#' @importFrom cluster pam
#' @importFrom factoextra hcut
#' @importFrom stats kmeans dist
#' @importFrom mclust Mclust
#' @importFrom sf st_as_sf st_bbox st_coordinates
#' @importFrom terra rast ext crds interpolate classify vect crs
#' @importFrom fields Tps
#' @importFrom ggplot2 ggplot geom_tile aes_string scale_fill_brewer labs
#'   theme_minimal
#' @import patchwork
#'
#' @examples
#' \dontrun{
#' ## Synthetic example (≈ 10 s)
#' set.seed(1)
#' n  <- 150
#' df <- data.frame(
#'   x  = runif(n, 23, 24),               # longitude
#'   y  = runif(n, -34, -33),             # latitude
#'   t1 = rnorm(n),                       # environmental predictors
#'   t2 = rnorm(n)
#' )
#'
#' res <- map_bioreg(
#'   data        = df,
#'   scale_cols  = c("t1", "t2"),
#'   clus_method = "kmeans",  # keep it light
#'   plot        = FALSE,
#'   interp      = "NN",
#'   resolution  = 0.2
#' )
#'
#' names(res)
#' head(res$data)
#' res$interpolated            # view with terra::plot(res$interpolated)
#' }
#'
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
