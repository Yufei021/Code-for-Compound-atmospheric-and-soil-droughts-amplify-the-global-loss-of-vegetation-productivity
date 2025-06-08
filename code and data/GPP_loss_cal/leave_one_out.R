leave_one_out_cv <- function(gpp_years, gpp_values) {
  errors <- numeric(length(gpp_values))
  
  for (i in 1:length(gpp_values)) {
    # Set aside the i-th data point for verification
    train_years <- gpp_years[-i]
    train_values <- gpp_values[-i]
    
    # Use the remaining data for fitting
    fit <- tryCatch({
      smooth.spline(train_years, train_values)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(fit)) {
      # Predict the left-out data points
      validation_year <- gpp_years[i]
      predicted_value <- predict(fit, validation_year)$y
      actual_value <- gpp_values[i]
      errors[i] <- actual_value - predicted_value
    }
  }
  
  return(errors)
}
