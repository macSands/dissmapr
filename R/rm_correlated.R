#' Remove Highly Correlated Variables
#'
#' This function identifies and removes variables with high correlation from a dataset.
#' It calculates the correlation matrix, visualizes it (optional), and removes variables
#' exceeding a specified correlation threshold.
#'
#' @param data Data frame containing the variables to analyze.
#' @param cols Numeric or character vector. Columns to consider for correlation analysis. Defaults to all columns.
#' @param threshold Numeric. Correlation threshold above which variables are considered highly correlated (default: 0.7).
#' @param plot Logical. Whether to plot the correlation matrix (default: TRUE).
#'
#' @return Data frame with highly correlated variables removed.
#' @export
#'
#' @examples
#' # Example dataset
#' data <- data.frame(
#'   var1 = rnorm(100),
#'   var2 = rnorm(100),
#'   var3 = rnorm(100) + 0.9 * rnorm(100)
#' )
#'
#' # Remove highly correlated variables
#' reduced_data <- rm_correlated(data, threshold = 0.8, plot = TRUE)
#'
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
