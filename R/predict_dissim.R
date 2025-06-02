#' Predict Pairwise Compositional Turnover (ζ-dissimilarity)
#'
#' @description
#' Generates spatial predictions of order-2 compositional turnover
#' (\eqn{\zeta_2}) from a fitted **zetadiv** multi-site generalised
#' dissimilarity model (`Zeta.msgdm`).  See
#' [calculate_pairwise_distances_matrix()] for the distance routine that this
#' function calls internally.
#'
#' @param block_sp  Data frame of site metadata; **must** contain a unique
#'   \code{grid_id} column and coordinate columns given by \code{x_col},
#'   \code{y_col}.
#' @param sbe_scaled Data frame or matrix of *scaled* environmental predictors
#'   (row order identical to \code{block_sp}).
#' @param zeta_model A fitted model returned by
#'   \code{zetadiv::Zeta.msgdm()}.
#' @param mean_rich Data frame with columns \code{order} and \code{value}
#'   giving the mean species richness (order = 1) for annotation.
#' @param mean_turn Data frame with columns \code{order} and \code{value}
#'   giving the mean turnover (order = 2) for annotation.
#' @param sbs_xy Data frame with site coordinates in the same order as
#'   \code{sbe_scaled}; must contain \code{x_col}, \code{y_col}.
#' @param rsa Optional \code{sf} or \code{SpatVector} polygon to draw on the
#'   figure.
#' @param x_col,y_col Names of the longitude / latitude columns
#'   (defaults \code{"x"}, \code{"y"}).
#'
#' @return
#' A data frame containing all predictors plus:
#' \describe{
#'   \item{distance}{Mean great-circle distance (km) to all other sites.}
#'   \item{rich_o1}{Mean richness (order 1) from \code{mean_rich}.}
#'   \item{turn_o2}{Mean turnover (order 2) from \code{mean_turn}.}
#'   \item{pred_zeta}{Linear predictor on the logit scale.}
#'   \item{pred_zetaExp}{Predicted turnover on the 0–1 scale.}
#'   \item{log_pred_zetaExp}{Natural-log of \code{pred_zetaExp}.}
#' }
#' A \code{ggplot2} heat-map is printed for visual inspection.
#'
#' @seealso
#' \code{\link{calculate_pairwise_distances_matrix}}
#'
#' @importFrom dplyr group_by summarize mutate %>%
#' @importFrom geosphere distHaversine
#' @importFrom zetadiv Zeta.msgdm Ispline Predict.msgdm
#' @importFrom ggplot2 ggplot aes_string geom_tile geom_text
#'   scale_fill_gradientn labs theme_minimal theme
#'
#' @examples
#' \dontrun{
#' ##----------------------------------------------------------
#' ## Toy example with 15 sites, two predictors, 60 species
#' ##----------------------------------------------------------
#' library(zetadiv)
#' set.seed(123)
#'
#' n_sites <- 15
#' block_sp <- data.frame(
#'   grid_id = sprintf("g%02d", 1:n_sites),
#'   x       = runif(n_sites, 22, 24),
#'   y       = runif(n_sites, -34, -33)
#' )
#'
#' env_raw <- data.frame(
#'   temp = rnorm(n_sites, 20, 2),
#'   rain = rnorm(n_sites, 600, 40)
#' )
#' sbe_scaled <- scale(env_raw) |> as.data.frame()
#'
#' ## Simulate a simple community matrix
#' comm <- matrix(rbinom(n_sites * 60, 1, 0.15),
#'                nrow = n_sites,
#'                dimnames = list(block_sp$grid_id, paste0("sp", 1:60)))
#'
#' ## Fit a Zeta.msgdm model (quick settings for example)
#' z_mod <- zetadiv::Zeta.msgdm(
#'   data.spec = comm,
#'   data.env  = env_raw,
#'   xy        = block_sp[, c("x", "y")],
#'   order     = 2,
#'   sam       = 100
#' )
#'
#' mean_rich <- data.frame(order = 1, value = mean(rowSums(comm)))
#' mean_turn <- data.frame(order = 2, value = mean(z_mod$zeta.val))
#'
#' res <- predict_dissim(
#'   block_sp   = block_sp,
#'   sbe_scaled = sbe_scaled,
#'   zeta_model = z_mod,
#'   mean_rich  = mean_rich,
#'   mean_turn  = mean_turn,
#'   sbs_xy     = block_sp
#' )
#'
#' head(res)
#' }
#'
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
