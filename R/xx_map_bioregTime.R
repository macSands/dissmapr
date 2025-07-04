library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cluster)      # pam()
library(factoextra)   # hcut()
library(mclust)       # Mclust()
library(NbClust)
library(terra)
library(fields)
library(patchwork)
library(clue)

# -------------------------------------------------------------------------
# cluster_all_scenarios --------------------------------------------------
# -------------------------------------------------------------------------
cluster_all_scenarios <- function(data,                 # list OR single df
                                  scale_cols,
                                  method      = c("kmeans", "pam",
                                                  "hclust", "gmm", "all"),
                                  k_override  = NULL,
                                  x_col       = "x",
                                  y_col       = "y",
                                  interpolate = c("none", "nn", "tps", "both"),
                                  res         = 0.5,
                                  crs         = "EPSG:4326") {

  ## ───────────────────────────── helpers ──────────────────────────────── ##
  align_to_ref <- function(ref, tgt) {
    if (length(unique(tgt)) != length(unique(ref)))
      tgt <- factor(tgt, levels = sort(unique(ref)))  # force same k
    M    <- table(ref, tgt)
    cost <- max(M) - M
    map  <- clue::solve_LSAP(cost)                    # Hungarian
    lut  <- setNames(as.integer(colnames(M)[map]),
                     rownames(M))                     # ref → matched tgt
    factor(lut[as.character(tgt)], levels = sort(unique(ref)))
  }

  mode_row <- function(r) {
    r <- na.omit(r)
    r[which.max(tabulate(match(r, r)))]
  }

  fill_nn_full <- function(df, x, y, cls, res, crs) {
    v  <- terra::vect(df, geom = c(x, y), crs = crs)
    r  <- terra::rast(terra::ext(v), resolution = res, crs = crs)
    terra::rasterize(v, r, field = cls, touches = TRUE)
  }

  fill_tps_full <- function(df, x, y, cls, res, crs) {
    v  <- terra::vect(df, geom = c(x, y), crs = crs)
    r  <- terra::rast(terra::ext(v), resolution = res, crs = crs)
    m  <- fields::Tps(df[, c(x, y)], as.integer(df[[cls]]))
    # terra::classify(terra::interpolate(r, m),
    #                 rcl = cbind(from = seq(0.5, max(df[[cls]]) - 0.5, 1),
    #                             to   = seq(1.5, max(df[[cls]]) + 0.5, 1),
    #                             becomes = 1:max(df[[cls]])))
    rcl <- data.frame(
      from    = 0.5 + 0:(cls-1),
      to      = 1.5 + 0:(cls-1),
      becomes = 1:cls
    )
    terra::classify(terra::interpolate(r, m), rcl = rcl)
  }
  ## ───────────────────────────── data prep ────────────────────────────── ##
  if (inherits(data, "data.frame")) {
    df_list <- list(current = data)
  } else if (is.list(data) && all(vapply(data, is.data.frame, TRUE))) {
    df_list <- data
    if (is.null(names(df_list)))
      names(df_list) <- paste0("set", seq_along(df_list))
  } else stop("`data` must be a data.frame or a list of data.frames")

  method      <- match.arg(method, several.ok = TRUE)
  if ("all" %in% method) method <- c("kmeans", "pam", "hclust", "gmm")
  interpolate <- match.arg(interpolate)

  dat  <- dplyr::bind_rows(df_list, .id = "scenario")
  stopifnot(all(c(x_col, y_col, scale_cols) %in% names(dat)))
  zmat <- scale(dat[, scale_cols])

  k <- if (!is.null(k_override)) k_override else
    NbClust::NbClust(zmat, "euclidean", 2, 10,
                     "kmeans", "silhouette")$Best.nc[1]

  ## ────────────────────────── run algorithms ─────────────────────────── ##
  if ("kmeans" %in% method)
    dat$kmeans <- kmeans(zmat, k, nstart = 10)$cluster
  if ("pam" %in% method)
    dat$pam    <- cluster::pam(zmat, k)$clustering
  if ("hclust" %in% method)
    dat$hclust <- factoextra::hcut(zmat, k, hc_method = "complete")$cluster
  if ("gmm" %in% method)
    dat$gmm    <- mclust::Mclust(zmat, G = k)$classification

  ## ─────────── align every non-kmeans column to the kmeans reference ─── ##
  algo_cols <- intersect(c("kmeans", "pam", "hclust", "gmm"), names(dat))
  for (col in setdiff(algo_cols, "kmeans"))
    dat[[col]] <- align_to_ref(dat$kmeans, dat[[col]])

  ## modal label across chosen algorithms
  dat$cluster_mode <- factor(apply(dat[algo_cols], 1, mode_row),
                             levels = 1:k)

  ## colour palette
  pal <- if (k <= 12) RColorBrewer::brewer.pal(k, "Set3") else
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(k)

  one_scn  <- length(unique(dat$scenario)) == 1
  mult_alg <- length(algo_cols) > 1

  plots <- if (one_scn && mult_alg) {           # single df, multiple algos
    lapply(algo_cols, function(col)
      ggplot2::ggplot(dat) +
        ggplot2::geom_tile(ggplot2::aes_string(
          x = x_col, y = y_col, fill = paste0("factor(", col, ")"))) +
        ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
        ggplot2::labs(title = paste("Algorithm:", col), fill = "Cluster") +
        ggplot2::theme_minimal())
  } else {                                     # default: facet by scenario
    dat |>
      dplyr::group_split(scenario) |>
      lapply(function(d)
        ggplot2::ggplot(d) +
          ggplot2::geom_tile(ggplot2::aes_string(
            x = x_col, y = y_col, fill = "cluster_mode")) +
          ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
          ggplot2::labs(title = unique(d$scenario),
                        subtitle = paste("Mode of:",
                                         paste(algo_cols, collapse = "/")),
                        fill = "Cluster") +
          ggplot2::theme_minimal())
  }
  print(patchwork::wrap_plots(plots))

  ## ───────────────────────── rasters (optional) ───────────────────────── ##
  make_stack <- function(fun) {
    lapply(split(dat, dat$scenario), function(d) {
      lays <- lapply(algo_cols, function(col)
        fun(d, x_col, y_col, col, res, crs))
      names(lays) <- algo_cols
      do.call(c, lays)                    # SpatRaster stack
    })
  }
  res_nn  <- res_tps <- NULL
  if (interpolate %in% c("nn",   "both")) res_nn  <- make_stack(fill_nn_full)
  if (interpolate %in% c("tps",  "both")) res_tps <- make_stack(fill_tps_full)

  ## ─────────────────────────── return ----------------------------------- ##
  list(
    nn      = res_nn,
    tps     = res_tps,
    table   = dat,
    plots   = plots,
    methods = algo_cols)
}

# -------------------------------------------------------------------------
# example call -------------------------------------------------------------
# -------------------------------------------------------------------------
## -------------------------------- single data frame
single_out <- cluster_all_scenarios(
  data        = predictors_df,
  scale_cols  = c("pred_zetaExp", "centroid_lon", "centroid_lat"),
  method      = c("kmeans","pam"),
  k_override  = 5,
  x_col       = "centroid_lon",
  y_col       = "centroid_lat",
  interpolate = "nn")

## -------------------------------- list of data frames
multi_out <- cluster_all_scenarios(
  data        = by_scn,      # named list: current, 2030, 2040, 2050
  scale_cols  = c("pred_zetaExp", "centroid_lon", "centroid_lat"),
  method      = "all",
  k_override  = 6,
  x_col       = "centroid_lon",
  y_col       = "centroid_lat",
  interpolate = "both")


multi_out$nn   # list of NN-filled SpatRaster stacks
multi_out$tps  # list of TPS-smoothed stacks
multi_out$table
multi_out$plots
