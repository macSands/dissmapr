#' Map Bioregion Clusters – single year or time series
#'
#' @param data A data frame **or** a *named list* of data frames
#'             (one per scenario when `multi = TRUE`).
#' @param multi Logical; `FALSE` = current behaviour, `TRUE` = treat `data`
#'              as multiple scenarios and enforce consistent cluster IDs.
#' @param ref   Name of the reference scenario in `data` used to pick *k*
#'              and label centroids (ignored if `multi = FALSE`).
#' @inheritParams map_bioreg
#'
#' @return If `multi = FALSE` a single list (as before).
#'         If `multi = TRUE` a named list of such lists, one per scenario.
#' @export
map_bioregTime <- function(data,
                           scale_cols,
                           clus_method = "all",
                           multi       = FALSE,
                           ref         = "current",
                           k_override  = NULL,
                           show_plot   = TRUE,
                           interp      = "both",
                           x_col       = "x",
                           y_col       = "y",
                           resolution  = 0.5,
                           bndy_fc     = NULL) {

  ## ---------------- helper to align one method -----------------------------
  align_one <- function(df, method_col, ref_cent) {
    tgt_cent <- aggregate(df[, c(x_col, y_col)],
                          by = list(cluster = df[[method_col]]), mean)
    dm <- as.matrix(dist(rbind(ref_cent, tgt_cent[,-1])))
    dm <- dm[seq_len(nrow(ref_cent)),
             -(seq_len(nrow(ref_cent)))]
    map <- apply(dm, 2, which.min)
    factor(map[df[[method_col]]],
           levels = seq_len(nrow(ref_cent)))
  }

  ## ---------------- function to run one scenario ---------------------------
  run_scenario <- function(df, k_fix = NULL, ref_list = NULL) {

    z <- scale(df[, scale_cols])
    k <- if (is.null(k_fix)) {
      NbClust(z, "euclidean", 2, 10, "kmeans", "silhouette")$Best.nc[1]
    } else k_fix

    out <- list()

    ## ---- k-means -----------------------------------------------------------
    if ("kmeans" %in% clus_method || "all" %in% clus_method) {
      km <- kmeans(z, centers = k, nstart = 10)
      df$kmeans_opt <- factor(km$cluster, levels = 1:k)
      out$kmeans_cent <- aggregate(df[, c(x_col, y_col)],
                                   by = list(cluster = df$kmeans_opt), mean)
    }

    ## ---- PAM ---------------------------------------------------------------
    if ("pam" %in% clus_method || "all" %in% clus_method) {
      pm <- cluster::pam(z, k)
      df$pam_opt <- factor(pm$clustering, levels = 1:k)
      df$pam_opt_aligned <- if (is.null(ref_list$pam_cent))
        df$pam_opt
      else
        align_one(df, "pam_opt", ref_list$pam_cent)
      out$pam_cent <- aggregate(df[, c(x_col, y_col)],
                                by = list(cluster = df$pam_opt_aligned), mean)
    }

    ## ---- hierarchical ------------------------------------------------------
    if ("hclust" %in% clus_method || "all" %in% clus_method) {
      hc <- factoextra::hcut(z, k, hc_method = "complete")
      df$hc_opt <- factor(hc$cluster, levels = 1:k)
      df$hc_opt_aligned <- if (is.null(ref_list$hc_cent))
        df$hc_opt
      else
        align_one(df, "hc_opt", ref_list$hc_cent)
      out$hc_cent <- aggregate(df[, c(x_col, y_col)],
                               by = list(cluster = df$hc_opt_aligned), mean)
    }

    ## ---- GMM ---------------------------------------------------------------
    if ("gmm" %in% clus_method || "all" %in% clus_method) {
      gm <- mclust::Mclust(z, G = k)
      df$gmm_opt <- factor(gm$classification, levels = 1:k)
      df$gmm_opt_aligned <- if (is.null(ref_list$gmm_cent))
        df$gmm_opt
      else
        align_one(df, "gmm_opt", ref_list$gmm_cent)
      out$gmm_cent <- aggregate(df[, c(x_col, y_col)],
                                by = list(cluster = df$gmm_opt_aligned), mean)
    }

    ## ---- quick plot (only existing columns) --------------------------------
    if (show_plot) {
      have <- grep("_opt", names(df), value = TRUE)
      mk <- function(col) {
        ggplot(df) +
          geom_tile(aes_string(x = x_col, y = y_col, fill = col)) +
          scale_fill_manual(values = palette_set3(k_fix %||% k), drop = FALSE) +
          labs(title = gsub("_.*", "", col), fill = "Cluster") +
          theme_minimal()
      }
      print(patchwork::wrap_plots(lapply(have, mk)))
    }

    out$data <- df
    out
  }

  palette_set3 <- function(k) {
    if (k <= 12) RColorBrewer::brewer.pal(k, "Set3")
    else         colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(k)
  }

  ## ---------------- single snapshot ----------------------------------------
  if (!multi) {
    return(run_scenario(data,
                        k_fix = k_override,
                        ref_list = list(),
                        palette = palette_set3(k_override %||% 10)))
  }

  ## ---------------- multi-year ---------------------------------------------
  if (!is.list(data))
    stop("With multi = TRUE, `data` must be a named list")

  ## reference ----------------------------------------------------------------
  ref_out <- run_scenario(data[[ref]],
                          k_fix = k_override,
                          ref_list = list())
  k_fix   <- k_override %||% nlevels(ref_out$data$kmeans_opt)

  ## others -------------------------------------------------------------------
  res <- purrr::imap(data, \(df, nm)
                     if (nm == ref) ref_out
                     else run_scenario(df, k_fix, ref_list = ref_out))

  return(res)
}
