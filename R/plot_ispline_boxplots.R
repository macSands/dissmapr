#' Plot facetted boxplots for all ispline basis columns
#'
#' @description
#' Creates a multi-facet boxplot of every ispline basis column in your data,
#' grouped by a specified order factor (e.g. zeta orders). Each spline term
#' (columns ending in \code{_is}) gets its own facet, and uses a
#' color-blind-friendly Viridis palette.
#'
#' @param ispline_data   A data frame as returned by \code{\link{run_ispline_models}},
#'                       containing raw covariates, their spline bases (suffixed \code{_is}),
#'                       and an order column.
#' @param ispline_suffix A string suffix identifying your spline columns.
#'                       Default is \dQuote{_is}.
#' @param order_col      The name of the grouping column (e.g. \dQuote{zOrder}).
#' @param palette        One of the Viridis options (\dQuote{viridis}, \dQuote{magma},
#'                       \dQuote{plasma}, \dQuote{cividis}, etc.). Default \dQuote{viridis}.
#' @param direction      Integer 1 or -1 to control palette direction. Default -1 (reversed).
#' @param ncol           Number of columns in the facet wrap. Default 3.
#' @param outlier_size   Size of the outlier points. Default 0.5.
#'
#' @return A \pkg{ggplot2} object showing one boxplot per spline term.
#'
#' @details
#' 1. Automatically detects all columns whose names end with \code{ispline_suffix}.
#' 2. Pivots the data to long format for \pkg{ggplot2}.
#' 3. Facets a boxplot for each spline variable with independent scales.
#'
#' @examples
#' \dontrun{
#' # load example data
#' data(bird.spec.fine)
#' data(bird.env.fine)
#'
#' # prepare inputs
#' xy.bird    <- bird.spec.fine[,1:2]
#' spp.bird   <- bird.spec.fine[,3:102]
#' envir.bird <- bird.env.fine[,3:9]
#'
#' # Fit & gather ispline tables
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
#' # Facetted boxplots of all *_is columns
#' plot_ispline_boxplots(
#'   ispline_data   = ispline_tabs_all,
#'   ispline_suffix = "_is",
#'   order_col      = "zOrder",
#'   palette        = "viridis",
#'   direction      = -1,
#'   ncol           = 3
#' )
#' }
#'
#' @seealso \code{\link{run_ispline_models}}, \code{\link{plot_ispline_lines}}
#' @export
plot_ispline_boxplots <- function(ispline_data,
                                  ispline_suffix = "_is",
                                  order_col      = "zOrder",
                                  palette        = "viridis",
                                  direction      = -1,
                                  ncol           = 3,
                                  outlier_size   = 0.5) {
  # detect spline columns
  spline_cols <- grep(paste0(ispline_suffix, "$"),
                      names(ispline_data),
                      value = TRUE)

  # pivot to long
  df_long <- ispline_data %>%
    dplyr::select(all_of(c(order_col, spline_cols))) %>%
    tidyr::pivot_longer(
      cols      = all_of(spline_cols),
      names_to  = "variable",
      values_to = "value"
    )

  # plot
  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = .data[[order_col]],
                               y = value,
                               fill = .data[[order_col]])) +
    ggplot2::geom_boxplot(outlier.size = outlier_size) +
    ggplot2::facet_wrap(~ variable, scales = "free", ncol = ncol) +
    ggplot2::scale_fill_viridis_d(option    = palette,
                                  direction = direction,
                                  name      = "Order\n(intensity ^)") +
    ggplot2::labs(
      x     = "Order",
      y     = "Value",
      title = "Distribution of ispline terms by Order"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text       = ggplot2::element_text(face = "bold"),
      legend.position  = "right"
    )
}
