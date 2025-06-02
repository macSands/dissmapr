#' Remove Highly Correlated Predictors
#'
#' @description
#' Detects pairs (or groups) of strongly collinear predictors and eliminates
#' the minimum subset necessary to keep every absolute pairwise correlation
#' below a user‐defined threshold.  Correlations are computed with
#' \code{stats::cor()}; the variables to discard are chosen via the
#' \code{caret::findCorrelation()} algorithm.  An optional heat-map of the
#' correlation matrix is produced with **corrplot** for rapid inspection.
#'
#' @details
#' * Non‐numeric columns are silently dropped prior to correlation
#'   calculation.
#' * When \code{cols} is supplied (numeric or character), only those columns are
#'   tested; otherwise all numeric columns in \code{data} are used.
#' * The names of removed and retained variables are printed to the console for
#'   transparency.
#'
#' @param data A data frame with candidate predictor variables.
#' @param cols Optional numeric or character vector specifying columns to
#'   consider; defaults to **all** numeric columns.
#' @param threshold Numeric in \eqn{[0,1]} specifying the absolute correlation
#'   cut-off (default \code{0.7}).
#' @param plot Logical; if \code{TRUE} (default) a correlation heat-map with
#'   coefficients is drawn.
#'
#' @return A data frame containing the original rows but only the subset of
#' predictor columns whose absolute pairwise correlations are \eqn{<}
#' \code{threshold}.
#'
#' @seealso
#' \code{\link[corrplot]{corrplot}}, \code{\link[caret]{findCorrelation}}
#'
#' @importFrom stats cor
#' @importFrom caret findCorrelation
#' @importFrom corrplot corrplot
#'
#' @examples
#' set.seed(99)
#' n  <- 200
#' df <- data.frame(
#'   a = rnorm(n),
#'   b = rnorm(n),
#'   c = rnorm(n) * 0.8 + rnorm(n) * 0.2,  # moderately corr. with 'a'
#'   d = rnorm(n) * 0.9 + rnorm(n) * 0.1   # moderately corr. with 'b'
#' )
#'
#' ## Remove predictors with |r| ≥ 0.75
#' df_reduced <- rm_correlated(df, threshold = 0.75, plot = FALSE)
#' names(df_reduced)
#'
#' ## Visualise the correlation structure & removals at a stricter threshold
#' rm_correlated(df, threshold = 0.6, plot = TRUE)
#'
#' @export
rm_correlated <- function(data, cols = NULL, threshold = 0.7, plot = TRUE) {
  library(caret)
  library(corrplot)

  # Select only specified columns or default to all columns in the data
  if (is.null(cols)) {
    vars <- data
  } else {
    vars <- data[, cols]
  }

  # Compute the correlation matrix
  cor_matrix <- cor(vars, use = "pairwise.complete.obs")

  # Optionally plot the correlation matrix
  if (plot) {
    corrplot(cor_matrix, method = "color", tl.cex = 0.6, tl.col = "black",
             addCoef.col = "black", number.cex = 0.4)
  }

  # Identify highly correlated variables
  highlyCorrelated <- findCorrelation(cor_matrix, cutoff = threshold, names = TRUE)

  # Remove highly correlated variables
  vars_reduced <- vars[, !names(vars) %in% highlyCorrelated]

  # Output the results
  cat("Variables removed due to high correlation:\n")
  print(highlyCorrelated)

  cat("\nVariables retained:\n")
  print(names(vars_reduced))

  # Return the reduced dataset
  return(vars_reduced)
}
