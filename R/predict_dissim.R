#' Predict Dissimilarity Using Environmental Data
#'
#' @param block_sp    Data frame containing geographic and environmental data; must include `grid_id`, and the columns named by `x_col` & `y_col`.
#' @param sbe_scaled  Scaled environmental predictors (same number of rows as `block_sp`).
#' @param zeta_model  Model object from `zetadiv::Zeta.msgdm`.
#' @param mean_rich   Data frame with mean richness values (`order`, `value`).
#' @param mean_turn   Data frame with mean turnover values (`order`, `value`).
#' @param sbs_xy      Data frame with the same sites in the same order as `sbe_scaled`; must include the columns named by `x_col` & `y_col`.
#' @param rsa         (Optional) an `sf` or `SpatVector` boundary layer to overlay (default: `NULL`).
#' @param x_col       Name of the longitude column in both `block_sp` & `sbs_xy` (default: `"x"`).
#' @param y_col       Name of the latitude  column in both `block_sp` & `sbs_xy` (default: `"y"`).
#'
#' @return A data frame with all predictors plus the dissimilarity predictions.
#' @export
predict_dissim <- function(
    block_sp,
    sbe_scaled,
    zeta_model,
    mean_rich,
    mean_turn,
    sbs_xy,
    rsa       = NULL,
    x_col     = "x",
    y_col     = "y"
) {
  library(geosphere)
  library(dplyr)
  library(zetadiv)
  library(ggplot2)

  #— 1) sanity checks
  req1 <- c("grid_id", x_col, y_col)
  if (!all(req1 %in% names(block_sp))) {
    stop("`block_sp` must have columns: ", paste(req1, collapse = ", "))
  }
  req2 <- c(x_col, y_col)
  if (!all(req2 %in% names(sbs_xy))) {
    stop("`sbs_xy` must have columns: ", paste(req2, collapse = ", "))
  }

  #— 2) compute pairwise distances (km), using user‐specified lon/lat
  dist_o2 <- calculate_pairwise_distances_matrix(
    data     = block_sp,
    x_col    = x_col,
    y_col    = y_col
  )

  #— 3) mean distance per site
  mean_dist <- dist_o2 %>%
    group_by(site_from) %>%
    summarize(distance = mean(value, na.rm = TRUE), .groups = "drop")

  #— 4) bring back coords onto mean_dist
  mean_dist[[x_col]] <- block_sp[[x_col]][
    match(mean_dist$site_from, block_sp$grid_id)
  ]
  mean_dist[[y_col]] <- block_sp[[y_col]][
    match(mean_dist$site_from, block_sp$grid_id)
  ]

  #— 5) build predictor table
  predictors_df <- sbe_scaled %>%
    mutate(distance = mean_dist$distance)

  #— 6) ispline + predict
  splines <- Ispline(predictors_df, order.ispline = 2)
  preds   <- Predict.msgdm(
    zeta_model$model,
    splines$splines,
    reg.type = "ispline",
    type     = "response"
  )

  #— 7) attach richness, turnover, raw & transformed preds
  predictors_df <- predictors_df %>%
    mutate(
      rich_o1           = mean_rich$value[mean_rich$order == 1],
      turn_o2           = mean_turn$value[mean_turn$order == 2],
      pred_zeta         = preds,
      pred_zetaExp      = exp(preds) / (1 + exp(preds)),
      log_pred_zetaExp  = log(pred_zetaExp),
      #—and now the spatial coords
      !!x_col           := sbs_xy[[x_col]],
      !!y_col           := sbs_xy[[y_col]]
    )

  #— 8) plot
  pal <- colorRampPalette(c("blue","green","yellow","orange","red","darkred"))

  p <- ggplot(predictors_df, aes_string(
    x    = x_col,
    y    = y_col,
    fill = "pred_zetaExp"
  )) +
    geom_tile() +
    geom_text(
      aes_string(label = "rich_o1"),
      color     = "yellow",
      size      = 1.5,
      # fontface  = "bold",
      check_overlap = TRUE
    ) +
    scale_fill_gradientn(colors = pal(10)) +
    theme_minimal() +
    labs(
      x    = x_col,
      y    = y_col,
      fill = "Predicted Turnover"
    ) +
    theme(
      panel.grid   = element_blank(),
      panel.border = element_blank()
    )

  if (!is.null(rsa)) {
    p <- p + geom_sf(
      data        = rsa,
      inherit.aes = FALSE,
      fill        = NA,
      color       = "black",
      alpha       = 0.5
    )
  }

  print(p)
  invisible(predictors_df)
}
