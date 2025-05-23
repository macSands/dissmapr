#' Predict Dissimilarity Using Environmental Data
#'
#' This function predicts turnover or zeta diversity values using models from `zetadiv::Zeta.msgdm` and new environmental data.
#'
#' @param block_sp Data frame containing geographic and environmental data.
#' @param sbe_scaled Scaled environmental predictors.
#' @param zeta_model Model object from `zetadiv::Zeta.msgdm`.
#' @param mean_rich Data frame with mean richness values for calibration.
#' @param mean_turn Data frame with mean turnover values for calibration.
#' @param sbs_xy Data frame with spatial coordinates.
#' @param rsa Spatial object for overlaying boundaries on plots.
#'
#' @return A data frame with predicted dissimilarity and other metrics.
#' @export
#'
#' @examples
#' # Example usage:
#' result <- predict_dissim(block_sp, sbe_scaled, zeta_model, mean_rich, mean_turn, sbs_xy, rsa)
#' head(result)
#'
predict_dissim <- function(block_sp,
                           sbe_scaled,
                           zeta_model,
                           mean_rich,
                           mean_turn,
                           sbs_xy,
                           rsa) {
  library(geosphere)
  library(dplyr)
  library(zetadiv)
  library(ggplot2)

  # Calculate geographic distance between all sites
  dist_o2 <- calculate_pairwise_distances_matrix(block_sp)

  # Summarize distance - mean distance per site
  mean_dist_o2 <- dist_o2 %>%
    group_by(site_from) %>%
    summarize(value = mean(value, na.rm = TRUE))

  # Add lon/lat coordinates back to mean distances
  mean_dist_o2$x <- block_sp$x[match(mean_dist_o2$site_from, block_sp$grid_id)]
  mean_dist_o2$y <- block_sp$y[match(mean_dist_o2$site_from, block_sp$grid_id)]

  # Combine data in a new predictors data frame
  predictors_df <- cbind(sbe_scaled, distance = mean_dist_o2$value)

  # Generate I-splines for predictors
  predictors_ispline <- zetadiv::Ispline(predictors_df, order.ispline = 2)

  # Predict turnover using the model
  predict_ispline <- zetadiv::Predict.msgdm(
    zeta_model$model,
    predictors_ispline$splines,
    reg.type = 'ispline',
    type = "response"
  )

  # Add species richness and turnover values to predictors_df
  predictors_df$rich_o1 <- as.vector(mean_rich$value[mean_rich$order == 1])
  predictors_df$turn_o2 <- as.vector(mean_turn$value[mean_turn$order == 2])
  predictors_df$pred_zeta <- predict_ispline
  predictors_df$pred_zetaExp <- exp(predict_ispline) / (1 + exp(predict_ispline))
  predictors_df$x <- sbs_xy$x
  predictors_df$y <- sbs_xy$y

  # Plot predicted turnover
  color_palette <- colorRampPalette(c("blue", "green", "yellow", "orange", "red", "darkred"))

  print(ggplot() +
    geom_tile(data = predictors_df, aes(x = x, y = y, fill = log(pred_zetaExp))) +
    scale_fill_gradientn(colors = color_palette(10)) +
    geom_sf(data = rsa, fill = NA, color = "black", alpha = 0.5) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude", fill = "Predicted Turnover") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank()
    ))

  # Return predictors_df with predictions
  return(predictors_df)
}
