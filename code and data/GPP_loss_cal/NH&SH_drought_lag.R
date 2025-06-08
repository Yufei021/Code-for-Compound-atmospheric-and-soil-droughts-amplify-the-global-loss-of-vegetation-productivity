# Calculate gpp_loss with lagged effect (i.e. gpp has been declining after the drought event)
# Output is NH_GPP_loss_lagged01.tif
# Load data and resize --------------------------------------------------------
# Load package
library(ncdf4)
library(raster)
library(reticulate)

# Setting up the Python environment
use_python("D:/anaconda3/python.exe", required = TRUE)

# Import the NumPy module
np <- import("numpy")

# Load .npy file
drought_events <- np$load("D:/paper01/VPD&SM(M03)/NH_only_vpd.npy")

# Resize the dimension of drought_eventsyears <- 37  
months_per_year <- 12

drought_events_reshaped <- array(0, dim = c(months_per_year, years, 180, 720))

# Reorganize the original drought_events by year and month
for (year in 1:years) {
  for (month in 1:months_per_year) {
    index_in_drought_events <- (year - 1) * 12 + month  
    drought_events_reshaped[month, year, , ] <- drought_events[index_in_drought_events, , ]
  }
}
dim(drought_events_reshaped)

# Load GPP
gpp_path <- "D:/GPP_data/GIMMS_NIRv/tif/WGS_1_tif_res_detr_deseason/"

original_gpp_array <- array(0, dim = c(444, 360, 720))

month_counter <- 1

# Loop through each month of each year
for (year in 1982:2018) {
  for (month in 1:12) {
    gpp_file <- sprintf("%s%04d_%02d.tif", gpp_path, year, month)
    
    if (file.exists(gpp_file)) {
      gpp_raster <- raster(gpp_file)
      
      # Check if the number of rows and columns of the raster data is 360 x 720
      if (nrow(gpp_raster) == 360 && ncol(gpp_raster) == 720) {
        # Extract all rows and columns of data
        gpp_data <- as.matrix(gpp_raster)
      } else {
        stop("GPP data dimensions do not meet expectations (360x720)")
      }
      
      # Store the data for that month in original_gpp_array
      original_gpp_array[month_counter, , ] <- gpp_data
      
      # Increment month counter
      month_counter <- month_counter + 1
    } else {
      warning(sprintf("File does not exist: %s", gpp_file))
    }
  }
}
# Output the dimensions of the 3D array to confirm
dim(original_gpp_array) 

gpp_array <- original_gpp_array[, 1:360, 1:720]
# Output the extracted dimensions to confirm
dim(gpp_array)

# Resize the dimensions of gpp_array
years <- 37  
months_per_year <- 12

gpp_array_reshaped <- array(0, dim = c(months_per_year, years, 360, 720))

# Rearrange by year and month
for (year in 1:years) {
  for (month in 1:months_per_year) {
    index_in_gpp_array <- (year - 1) * 12 + month  
    gpp_array_reshaped[month, year, , ] <- gpp_array[index_in_gpp_array, , ]
  }
}
dim(gpp_array_reshaped)


# Create a new array to store labels
drought_events_reshaped02 <- drought_events_reshaped

# Iterate through each pixel to detect drought events
for (lat in 1:360) {
  for (lon in 1:720) {
    # Find the start and end indices of drought events
    drought_sequence <- drought_events_reshaped[, , lat, lon]
    drought_start_indices <- which(drought_sequence == 1, arr.ind = TRUE)
    
    if (nrow(drought_start_indices) > 1) {  # Proceed only if there are multiple drought events
      for (i in 1:(nrow(drought_start_indices) - 1)) {  # Ensure no out-of-bounds error
        start <- drought_start_indices[i, ]
        end <- drought_start_indices[i + 1, ]
        
        start_month <- start[1]  # Starting month
        start_year <- start[2]   # Starting year
        next_start_month <- end[1]  # Starting month of next drought event
        next_start_year <- end[2]   # Starting year of next drought event
        
        # Calculate the time difference (in months) between two drought events
        start_index <- (start_year - 1) * 12 + start_month
        next_start_index <- (next_start_year - 1) * 12 + next_start_month
        
        if (start_index + 1 < next_start_index) {
          
          gpp_values <- vector()
          for (index in start_index:(next_start_index - 1)) {
            month <- (index - 1) %% 12 + 1
            year <- floor((index - 1) / 12) + 1
            gpp_values <- c(gpp_values, gpp_array_reshaped[month, year, lat, lon])
          }
          
          # Limit to a maximum of 12 months of data
          gpp_values <- gpp_values[1:min(12, length(gpp_values))]
          
          # Check for valid GPP data and identify continuously declining months
          if (length(gpp_values) > 1 && sum(!is.na(gpp_values)) == length(gpp_values)) {
            # Try to identify the start of a persistent decline
            found_decline_start <- FALSE
            decline_start_index <- 1
            
            for (start_point in 1:(length(gpp_values) - 1)) {
              if (!found_decline_start) {
                # Try to fit a model from each potential starting point to the end
                for (end_point in (start_point + 1):length(gpp_values)) {
                  time_points <- seq(start_point, end_point)
                  test_fit <- lm(gpp_values[start_point:end_point] ~ time_points)
                  
                  # Check if there's an overall downward trend
                  if (nrow(summary(test_fit)$coefficients) >= 2 && summary(test_fit)$coefficients[2, 1] < 0) {
                    decline_start_index <- start_point
                    found_decline_start <- TRUE
                    break
                  }
                }
              } else {
                break
              }
            }
            
            # If the start of a continuous decline is found
            if (found_decline_start) {
              # Mark all months from the start point to the end of data as 2
              for (index in (start_index + decline_start_index - 1):(start_index + length(gpp_values))) {
                month <- (index - 1) %% 12 + 1
                year <- floor((index - 1) / 12) + 1
                drought_events_reshaped02[month, year, lat, lon] <- 2
              }
            }
          }
          
          # Check valid GPP data and perform segmented linear regression to detect decline then rise
          if (length(gpp_values) > 1 && sum(!is.na(gpp_values)) == length(gpp_values)) {
            # Try to identify the start of a persistent decline
            found_decline_start <- FALSE
            found_subsequent_rise <- FALSE
            decline_start_index <- 1
            
            for (start_point in 1:(length(gpp_values) - 1)) {
              if (!found_decline_start) {
                # Try to fit a model from each potential starting point to the end
                for (end_point in (start_point + 1):length(gpp_values)) {
                  time_points <- seq(start_point, end_point)
                  test_fit <- lm(gpp_values[start_point:end_point] ~ time_points)
                  
                  # Check if there's an overall downward trend
                  if (nrow(summary(test_fit)$coefficients) >= 2 && summary(test_fit)$coefficients[2, 1] < 0) {
                    decline_start_index <- start_point
                    found_decline_start <- TRUE
                    break
                  }
                }
              } else {
                break
              }
            }
            
            # Check if a subsequent rise exists
            if (found_decline_start) {
              # Try fitting a linear model from end of decline to the end
              if (decline_start_index < length(gpp_values)) {
                rise_start_point <- decline_start_index + 1
                time_points <- seq(rise_start_point, length(gpp_values))
                rise_fit <- lm(gpp_values[rise_start_point:length(gpp_values)] ~ time_points)
                if (nrow(summary(rise_fit)$coefficients) >= 2 && summary(rise_fit)$coefficients[2, 1] > 0) {
                  found_subsequent_rise <- TRUE
                }
              }
              
              # Mark as "decline then rise"
              if (found_subsequent_rise) {
                for (index in (start_index + decline_start_index - 1):(start_index + length(gpp_values))) {
                  month <- (index - 1) %% 12 + 1
                  year <- floor((index - 1) / 12) + 1
                  drought_events_reshaped02[month, year, lat, lon] <- 3
                }
              }
            }
          }
        }
      }
    }
  }
}


# Export the array to .R file -----------------------------------
# Define output filename
output_file <- "G:/VPD&SM_CDE(M03)/NH_only_sm_reshaped.RData"

# Save the array to an .RData file
save(drought_events_reshaped02, file = output_file)

# Notify that saving is complete
cat("Array saved to file:", output_file, "\n")







