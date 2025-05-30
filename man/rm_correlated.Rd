% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rm_correlated.R
\name{rm_correlated}
\alias{rm_correlated}
\title{Remove Highly Correlated Variables}
\usage{
rm_correlated(data, cols = NULL, threshold = 0.7, plot = TRUE)
}
\arguments{
\item{data}{Data frame containing the variables to analyze.}

\item{cols}{Numeric or character vector. Columns to consider for correlation analysis. Defaults to all columns.}

\item{threshold}{Numeric. Correlation threshold above which variables are considered highly correlated (default: 0.7).}

\item{plot}{Logical. Whether to plot the correlation matrix (default: TRUE).}
}
\value{
Data frame with highly correlated variables removed.
}
\description{
This function identifies and removes variables with high correlation from a dataset.
It calculates the correlation matrix, visualizes it (optional), and removes variables
exceeding a specified correlation threshold.
}
\examples{
# Example dataset
data <- data.frame(
  var1 = rnorm(100),
  var2 = rnorm(100),
  var3 = rnorm(100) + 0.9 * rnorm(100)
)

# Remove highly correlated variables
reduced_data <- rm_correlated(data, threshold = 0.8, plot = TRUE)

}
