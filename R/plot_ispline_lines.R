#' Plot ispline partial effects with quantile and start-point markers
#'
#' @description
#' Given a tidy data frame of ispline basis outputs (from \code{\link{run_ispline_models}}),
#' this function identifies the ispline column matching a specified covariate,
#' draws thinner ispline curves for each \code{zOrder}, and overlays:
#' \itemize{
#'   \item Small symbols at user-defined quantiles of the raw covariate.
#'   \item A larger symbol at each curve's starting point (minimum covariate value).
#' }
#'
#' @param ispline_data A data frame as returned by \code{\link{run_ispline_models}},
#'   containing raw covariates, their `_is` spline bases, and a \code{zOrder} column.
#' @param x_var A string; the name of the raw covariate to plot (e.g. \dQuote{dist}).
#'   The function will look for a matching spline column named \code{<x_var>_is}.
#' @param orders Character vector of \code{zOrder} levels (in desired legend/order).
#'   Default is \code{unique(ispline_data$zOrder)}.
#' @param cols   Character vector of colours, one per order.
#'   Default uses \code{scales::hue_pal()}.
#' @param shapes Integer vector of plotting symbols (pch codes), one per order.
#'   Default is \code{15:}(...).
#' @param probs  Numeric vector of probabilities (between 0 and 1) at which to place
#'   small quantile markers. Default \code{c(0, .25, .5, .75, 1)}.
#' @param line_size    Numeric; line width for the spline curves. Default 0.5.
#' @param point_size   Numeric; size of the quantile markers. Default 1.5.
#' @param start_size   Numeric; size of the big start-point markers. Default 3.
#' @param start_stroke Numeric; stroke width for the big start markers. Default 0.
#'
#' @return A \pkg{ggplot2} object showing:
#' \itemize{
#'   \item \code{geom_line()} for each curve.
#'   \item \code{geom_point()} at the specified quantiles.
#'   \item A larger \code{geom_point()} at each curve's minimum.
#' }
#'
#' @details
#' 1. Detects all columns ending in \code{_is} and fuzzy-matches \code{x_var} against them.
#' 2. Extracts the raw covariate name by stripping \code{_is}.
#' 3. Computes one "start" point per \code{zOrder} at the minimum raw covariate.
#' 4. Computes quantile-closest points per \code{zOrder} at the user's \code{probs}.
#'
#' @examples
#' \dontrun{
#' # load sample data
#' data(bird.spec.fine)
#' data(bird.env.fine)
#'
#' # prepare inputs
#' xy.bird    <- bird.spec.fine[,1:2]
#' spp.bird   <- bird.spec.fine[,3:102]
#' envir.bird <- bird.env.fine[,3:9]
#'
#' # fit & gather ispline tables
#' ispline_tabs_all <- run_ispline_models(
#'   spp_df    = spp.bird,
#'   env_df    = envir.bird,
#'   xy_df     = xy.bird,
#'   orders    = 2:6,
#'   sam       = 100,
#'   normalize = "Jaccard",
#'   reg_type  = "ispline"
#' )
#'
#' # Line plot for "dist"
#' plot_ispline_lines(
#'   ispline_data = ispline_tabs_all,
#'   x_var        = "dist",            # will match "dist_is"
#'   orders       = paste("Order", 2:6),
#'   cols         = c('green','cyan','purple','blue','black'),
#'   shapes       = c(15,16,17,18,19)
#' )
#' }
#'
#' @seealso
#' \code{\link{run_ispline_models}}, \code{\link[ggplot2]{ggplot}}
#'
#' @export
plot_ispline_lines <- function(ispline_data,
                               x_var,
                               orders    = NULL,
                               cols      = NULL,
                               shapes    = NULL,
                               probs     = c(0, .25, .5, .75, 1),
                               line_size    = 0.5,
                               point_size   = 1.5,
                               start_size   = 3,
                               start_stroke = 0) {
  # find all *_is columns
  spline_cols <- grep("_is$", names(ispline_data), value = TRUE)
  # match x_var to a single spline column
  matches <- spline_cols[grepl(x_var, spline_cols, ignore.case = TRUE)]
  if (length(matches) != 1) {
    stop("`x_var = '", x_var, "'` must match exactly one of: ",
         paste(spline_cols, collapse = ", "))
  }
  spline_col <- matches
  raw_col    <- sub("_is$", "", spline_col)
  if (!raw_col %in% names(ispline_data)) {
    stop("Raw covariate '", raw_col, "' not found.")
  }

  # defaults
  if (is.null(orders)) orders <- unique(ispline_data$zOrder)
  if (is.null(cols))   cols   <- scales::hue_pal()(length(orders))
  if (is.null(shapes)) shapes <- seq(15, length.out = length(orders))

  # compute start & quantile points
  start_pts <- ispline_data %>%
    dplyr::group_by(zOrder) %>%
    dplyr::slice_min(.data[[raw_col]], n = 1) %>%
    dplyr::ungroup()

  quant_pts <- ispline_data %>%
    dplyr::group_split(zOrder) %>%
    purrr::map_dfr(function(df) {
      z   <- unique(df$zOrder)
      qs  <- quantile(df[[raw_col]], probs)
      idx <- sapply(qs, function(q) which.min(abs(df[[raw_col]] - q)))
      tibble::tibble(
        !!raw_col    := df[[raw_col]][idx],
        !!spline_col := df[[spline_col]][idx],
        zOrder       = z
      )
    })

  # plot
  ggplot2::ggplot(ispline_data,
                  ggplot2::aes(x = .data[[raw_col]],
                               y = .data[[spline_col]],
                               color = zOrder)) +
    ggplot2::geom_line(size = line_size) +
    ggplot2::geom_point(data   = quant_pts,
                        ggplot2::aes(shape = zOrder),
                        size   = point_size,
                        stroke = 1) +
    ggplot2::geom_point(data   = start_pts,
                        ggplot2::aes(x = .data[[raw_col]],
                                     y = .data[[spline_col]],
                                     shape = zOrder),
                        size   = start_size,
                        stroke = start_stroke) +
    ggplot2::scale_color_manual(values = cols,   labels = orders) +
    ggplot2::scale_shape_manual(values = shapes, labels = orders) +
    ggplot2::labs(x     = raw_col,
                  y     = spline_col,
                  title = paste("Ispline effect of", raw_col, "by zOrder")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title    = element_blank(),
                   legend.position = "right")
}
