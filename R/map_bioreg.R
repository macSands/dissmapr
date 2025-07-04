#' Map Bioregion Clusters and Spatially Interpolate Assignments
#'
#' @description
#' Runs several unsupervised clustering algorithms on multivariate
#' environmental data, aligns their cluster labels to a common
#' k-means reference, visualises the resulting partitions (optionally),
#' and produces gridded bioregion surfaces by nearest-neighbour (NN)
#' and/or thin-plate-spline (TPS) interpolation.
#'
#' @param data         Data frame containing longitude/latitude plus the
#'                     variables in `scale_cols`.
#' @param scale_cols   Character vector of column names to z-scale before
#'                     clustering.
#' @param clus_method  Which algorithm(s) to run: one or more of
#'                     `"kmeans"`, `"pam"`, `"hclust"`, `"gmm"`, or `"all"`.
#' @param show_plot    Logical; if `TRUE` (default), draw a 2×2 patchwork
#'                     of the point-map partitions. Set `FALSE` to skip.
#' @param interp       Character; interpolation type: `"NN"`, `"Tps"`, or
#'                     `"both"` (default).
#' @param x_col,y_col  Names of the longitude and latitude columns
#'                     (defaults `"x"`, `"y"`).
#' @param resolution   Numeric; grid cell size for the interpolation raster.
#' @param bndy_fc      Optional \code{sf} or \code{SpatVector} polygon to overlay.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{data}}{Original `data` with factor cluster columns
#'         (`kmeans_opt`, `pam_opt`, `pam_opt_aligned`,
#'         `hc_opt`, `hc_opt_aligned`, `gmm_opt`, `gmm_opt_aligned`),
#'         all sharing levels 1…k.}
#'   \item{\code{interpolated}}{A `SpatRaster` of continuous NN or TPS
#'         predictions (one layer per algorithm).}
#'   \item{\code{clusters}}{A `SpatRaster` of re-classified discrete
#'         clusters (TPS only), with the same levels 1…k.}
#' }
#'
#' @importFrom NbClust NbClust
#' @importFrom cluster pam
#' @importFrom factoextra hcut
#' @importFrom stats kmeans dist
#' @importFrom mclust Mclust mclustBIC
#' @importFrom sf st_as_sf st_bbox st_coordinates
#' @importFrom terra rast ext crds interpolate classify vect crs
#' @importFrom fields Tps
#' @importFrom ggplot2 ggplot geom_tile aes_string scale_fill_brewer labs theme_minimal
#' @import patchwork
#' @export
map_bioreg = function(data,
                       scale_cols,
                       clus_method = "all",
                       show_plot   = TRUE,
                       interp      = "both",
                       x_col       = "x",
                       y_col       = "y",
                       resolution  = 0.5,
                       bndy_fc     = NULL) {
  #--- sanity
  if (!all(c(x_col, y_col, scale_cols) %in% names(data))) {
    stop("`data` must contain: ", paste(c(x_col,y_col,scale_cols), collapse=", "))
  }

  #--- helper to align labels
  align_clusters = function(df, ref_col, tgt_col) {
    ref_ctr = aggregate(df[,c(x_col,y_col)], by=list(cluster=df[[ref_col]]), FUN=mean)
    tgt_ctr = aggregate(df[,c(x_col,y_col)], by=list(cluster=df[[tgt_col]]), FUN=mean)
    dm      = as.matrix(dist(rbind(ref_ctr[,-1], tgt_ctr[,-1])))
    dm      = dm[1:nrow(ref_ctr), (nrow(ref_ctr)+1):ncol(dm)]
    map_idx = apply(dm, 2, which.min)
    factor(map_idx[df[[tgt_col]]], levels = seq_len(nrow(ref_ctr)))
  }

  #--- 1. scale
  data_scaled = scale(data[, scale_cols])

  #--- 2. choose k via silhouette on k-means
  if (!("kmeans" %in% clus_method || "all" %in% clus_method)) {
    stop("`clus_method` must include `kmeans` for reference")
  }
  nc = NbClust(data_scaled, distance="euclidean",
                min.nc=2, max.nc=10,
                method="kmeans", index="silhouette")
  k  = nc$Best.nc[1]

  #--- 3. clustering w/ fixed factor levels 1:k
  km = kmeans(data_scaled, centers=k, nstart=10)
  data$kmeans_opt = factor(km$cluster, levels = 1:k)

  if ("pam" %in% clus_method || "all" %in% clus_method) {
    pm = pam(data_scaled, k = k)
    data$pam_opt         = factor(pm$clustering, levels = 1:k)
    data$pam_opt_aligned = align_clusters(data, "kmeans_opt", "pam_opt")
  }
  if ("hclust" %in% clus_method || "all" %in% clus_method) {
    hc = hcut(data_scaled, k = k, hc_method = "complete")
    data$hc_opt         = factor(hc$cluster, levels = 1:k)
    data$hc_opt_aligned = align_clusters(data, "kmeans_opt", "hc_opt")
  }
  if ("gmm" %in% clus_method || "all" %in% clus_method) {
    gm = mclust::Mclust(data_scaled, G = k)
    data$gmm_opt         = factor(gm$classification, levels = 1:k)
    data$gmm_opt_aligned = align_clusters(data, "kmeans_opt", "gmm_opt")
  }

  #--- 4. optional plot
  if (show_plot) {

    ## helper that builds one map, optionally adds a boundary -----------------
    mk = function(fill_col, title) {

      p = ggplot(data) +
        geom_tile(aes_string(x = x_col, y = y_col, fill = fill_col)) +
        scale_fill_brewer(palette = "Set3", drop = FALSE) +
        labs(title = title, x = "Longitude", y = "Latitude",
             fill = "Cluster") +
        theme_minimal()

      # ---- add the boundary if one was supplied -----------------------------
      if (!is.null(bndy_fc)) {
        ## convert sp -> sf on-the-fly if necessary
        if (inherits(bndy_fc, "Spatial")) bndy_fc = sf::st_as_sf(bndy_fc)

        p = p +
          geom_sf(data = bndy_fc,            # outline only
                  fill  = NA,
                  colour = "black",
                  linewidth = .4)
      }
      p
    }

    ## build the four possible panels ----------------------------------------
    p1 = mk("kmeans_opt"      , "K-means")
    p2 = if ("pam_opt_aligned" %in% names(data)) mk("pam_opt_aligned", "PAM")
    p3 = if ("hc_opt_aligned"  %in% names(data)) mk("hc_opt_aligned" , "Hierarchical")
    p4 = if ("gmm_opt_aligned" %in% names(data)) mk("gmm_opt_aligned", "GMM")

    ## print the figure (keeps old layout)
    print((p1 + p2) / (p3 + p4))
  }

  #--- 5. prepare grid
  pts  = sf::st_as_sf(data, coords = c(x_col, y_col), crs = 4326)
  bb   = sf::st_bbox(pts)
  base = rast(ext(bb$xmin, bb$xmax, bb$ymin, bb$ymax),
               resolution = resolution,
               crs        = "EPSG:4326")

  #--- 6. NN interpolation
  interpolated = NULL
  if (interp %in% c("NN","both")) {
    crd = sf::st_coordinates(pts)
    for (col in c("kmeans_opt","pam_opt_aligned","hc_opt_aligned","gmm_opt_aligned")) {
      if (!col %in% names(data)) next
      vals = apply(crds(base), 1, function(pt) {
        i = which.min((crd[,1]-pt[1])^2 + (crd[,2]-pt[2])^2)
        as.integer(as.character(data[[col]][i]))
      })
      base[[col]] = vals
    }
    interpolated = base
  }

  #--- 7. TPS interpolation + dynamic reclass
  clusters_rcl = NULL
  if (interp %in% c("Tps","both")) {
    grid_r = rast(ext(pts), resolution=resolution, crs=crs(base))
    tps_mods = list(
      kmeans = Tps(data[,c(x_col,y_col)], as.integer(data$kmeans_opt)),
      pam    = if("pam_opt_aligned"%in%names(data)) Tps(data[,c(x_col,y_col)], as.integer(data$pam_opt_aligned)),
      hc     = if("hc_opt_aligned"%in%names(data)) Tps(data[,c(x_col,y_col)], as.integer(data$hc_opt_aligned)),
      gmm    = if("gmm_opt_aligned"%in%names(data)) Tps(data[,c(x_col,y_col)], as.integer(data$gmm_opt_aligned))
    )
    preds = lapply(tps_mods, function(m) if(!is.null(m)) terra::interpolate(grid_r,m) else NULL)
    valid = names(preds)[!vapply(preds, is.null, logical(1))]
    stack = do.call(c, preds[valid])
    names(stack) = valid
    rcl = cbind(
      from    = seq(0.5, k-0.5, by=1),
      to      = seq(1.5, k+0.5, by=1),
      becomes = 1:k
    )
    clusters_rcl = do.call(c, lapply(stack, classify, rcl=rcl))
  }

  list(
    data         = data,
    interpolated = interpolated,
    clusters     = clusters_rcl
  )
}
