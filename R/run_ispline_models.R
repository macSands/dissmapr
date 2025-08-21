#' Run multiple Zeta.msgdm ispline models and return both models and combined ispline table
#'
#' @description
#' Fits \code{Zeta.msgdm} models of type \dQuote{ispline} for a series of zeta-orders,
#' extracts the raw environmental covariates (plus distance) and their ispline bases,
#' and returns both the list of fitted models and one tidy data frame combining all orders.
#'
#' @param spp_df        A data frame or matrix of species incidence/abundance.
#' @param env_df        A data frame of environmental covariates.
#' @param xy_df         A two-column data frame or matrix of site coordinates.
#' @param orders        Integer vector of zeta orders to fit (e.g. 2:6).
#' @param sam           Integer; number of random samples per order (passed to \code{Zeta.msgdm}).
#' @param distance.type Character; distance metric for \code{Zeta.msgdm} (default \dQuote{Euclidean}).
#' @param normalize     Character; normalization method for \code{Zeta.msgdm} (default \dQuote{Jaccard}).
#' @param reg_type      Character; regression type for \code{Zeta.msgdm} (default \dQuote{ispline}).
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{zeta_gdm_list}}{A list of the fitted \code{Zeta.msgdm()} objects, named by order.}
#'   \item{\code{ispline_table}}{A tibble with one row per sample, containing all
#'     original covariates (including \code{distance}), the ispline bases suffixed
#'     \code{_is}, and a \code{zOrder} column.}
#' }
#'
#' @examples
#' \dontrun{
#' data(bird.spec.fine); data(bird.env.fine)
#' xy   <- bird.spec.fine[,1:2]
#' spp  <- bird.spec.fine[,3:102]
#' env  <- bird.env.fine[,3:9]
#'
#' out <- run_ispline_models(
#'   spp_df        = spp,
#'   env_df        = env,
#'   xy_df         = xy,
#'   orders        = 2:6,
#'   sam           = 100,
#'   normalize     = "Jaccard",
#'   reg_type      = "ispline"
#' )
#' names(out)
#' head(out$ispline_table)
#' }
#'
#' @seealso
#' \code{\link[zetadiv]{Zeta.msgdm}}, \code{\link[zetadiv]{Return.ispline}}
#' @export
run_ispline_models <- function(
    spp_df,
    env_df,
    xy_df,
    orders         = 2:6,
    sam            = 100,
    distance.type  = "Euclidean",
    normalize      = "Jaccard",
    reg_type       = "ispline"
) {
  # 1. Fit one Zeta.msgdm model per order
  zeta_gdm_list <- purrr::map(orders, function(ord) {
    zetadiv::Zeta.msgdm(
      spp_df,
      env_df,
      xy_df,
      sam           = sam,
      order         = ord,
      distance.type = distance.type,
      normalize     = normalize,
      reg.type      = reg_type
    )
  })
  names(zeta_gdm_list) <- paste0("Order", orders)

  # 2. Extract ispline tables and bind into one tibble
  ispline_table <- purrr::map2_dfr(
    zeta_gdm_list, names(zeta_gdm_list),
    function(zeta_obj, ord_name) {
      out     <- zetadiv::Return.ispline(zeta_obj, data.env = env_df, distance = TRUE)
      env_tab <- out$env
      iso_tab <- out$Ispline
      colnames(iso_tab) <- paste0(colnames(iso_tab), "_is")
      dplyr::bind_cols(env_tab, iso_tab) %>%
        dplyr::mutate(zOrder = ord_name)
    }
  )

  # 3. Return both the model list and the combined table
  list(
    zeta_gdm_list = zeta_gdm_list,
    ispline_table = ispline_table
  )
}
