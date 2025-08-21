# Add a tiny helper for the test suite
dfify <- function(x) if (is.data.frame(x)) x else as.data.frame(x)
