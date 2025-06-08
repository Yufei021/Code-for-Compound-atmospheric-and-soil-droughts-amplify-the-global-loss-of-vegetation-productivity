# Calculate GPP loss with lag effect (i.e., GPP continues to decline after drought events)
# Use smooth.spline function with leave-one-out cross-validation
# Running the code: leave_one_out.R at first

# Load data and reshape dimensions --------------------------------------------------------
# Load libraries
library(ncdf4)
library(raster)
library(reticulate)

# Load .RData file (the results of running the NH&SH_drought_lag.R code)
load("G:/paper01/TEST/only_events_reshaped.RData")

# Read GPP data
# Path to GPP data
gpp_path <- "G:/paper01/GPP_data/GIMMS_NIRv/res/tif/WGS_1_tif_res_detr_deseason/"

# Initialize original_gpp_array with dimensions [444, 360, 720]
original_gpp_array <- array(0, dim = c(444, 360, 720))

# Define a counter to track the total number of months
month_counter <- 1

# Loop through each year and month
for (year in 1982:2018) {
  for (month in 1:12) {
    # Construct file name, e.g., "NIRv.GPP.198201.v1.nc"
    gpp_file <- sprintf("%s%04d_%02d.tif", gpp_path, year, month)
    
    # Read GPP data for the month (if the file exists)
    if (file.exists(gpp_file)) {
      gpp_raster <- raster(gpp_file)
      
      # Check if raster has 360 rows and 720 columns
      if (nrow(gpp_raster) == 360 && ncol(gpp_raster) == 720) {
        # Extract all data as a matrix
        gpp_data <- as.matrix(gpp_raster)
      } else {
        stop("GPP data dimensions do not match expected 360x720")
      }
      
      # Store data into original_gpp_array (3D array: [total_months, lat, lon])
      original_gpp_array[month_counter, , ] <- gpp_data
      
      # Increment month counter
      month_counter <- month_counter + 1
    } else {
      warning(sprintf("File not found: %s", gpp_file))
    }
  }
}
# Output the dimensions of the 3D array to confirm
dim(original_gpp_array) 

gpp_array <- original_gpp_array[, 1:360, 1:720]
# Output dimensions of extracted data to confirm
dim(gpp_array)

# Reshape gpp_array dimensions
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

# Fit normal-condition GPP values and calculate GPP Loss -------------------------------------------
library(tcltk)

# Initialize 3D GPP loss array
gpp_loss_array <- array(0, dim = c(12, 360, 720))  # 12 months, 360 rows, 720 columns

# Create a total difference matrix
total_diff <- array(0, dim = c(360, 720))  # To store total difference for each pixel

# Initialize successful_fits list
successful_fits <- list()

# Total number of steps: 12 months * 360 rows * 720 columns
total_steps <- 12 * 360 * 720

# Create progress bar
pb <- tkProgressBar(title = "Processing", min = 0, max = total_steps, width = 300)

# Current step counter
step_count <- 0

# Pre-compute non-drought and drought year indices to avoid repeated computation
non_drought_indices <- drought_events_reshaped02 == 0  # Indices for non-drought years
drought_indices <- drought_events_reshaped02 > 0        # Indices for drought years

# Use expand.grid to generate lat-lon combinations
lat_lon_combinations <- expand.grid(lat = 1:360, lon = 1:720)

# Initialize gpp_differences list
gpp_differences <- list()

# Use lapply to iterate over each month
for (month in 1:12) {
  month_diff <- array(0, dim = c(37, 360, 720))  # Array to store GPP differences for drought years per pixel
  
  # Use lapply to iterate over each latitude-longitude combination
  results <- lapply(1:nrow(lat_lon_combinations), function(index) {
    lat <- lat_lon_combinations$lat[index]
    lon <- lat_lon_combinations$lon[index]
    
    tryCatch({
      # Update progress bar
      step_count <<- step_count + 1
      setTkProgressBar(pb, step_count, label = paste("Processing", step_count, "of", total_steps))
      
      # Get GPP values for this month during non-drought years
      non_drought_gpp_values <- gpp_array_reshaped[month, non_drought_indices[month, , lat, lon], lat, lon]
      non_drought_years <- 1982 + which(non_drought_indices[month, , lat, lon]) - 1
      
      # Remove NA and Inf values from non-drought data
      valid_indices <- which(!is.na(non_drought_gpp_values) & is.finite(non_drought_gpp_values))
      non_drought_gpp_values <- non_drought_gpp_values[valid_indices]
      non_drought_years <- non_drought_years[valid_indices]
      
      num_points <- length(non_drought_gpp_values)
      
      # Skip this pixel if data points are too few or not enough years
      if (num_points < 3 || length(unique(non_drought_years)) < 4) {
        return(NULL)
      }
      
      # Perform leave-one-out cross-validation to calculate errors
      errors <- leave_one_out_cv(non_drought_years, non_drought_gpp_values)
      mean_error <- mean(errors, na.rm = TRUE)
      
      # Calculate GPP differences for drought years (values > 0)
      drought_gpp_values <- gpp_array_reshaped[month, drought_indices[month, , lat, lon], lat, lon]
      drought_years <- 1982 + which(drought_indices[month, , lat, lon]) - 1
      
      # Fit a smoothing spline model using non-drought GPP values
      spline_model <- smooth.spline(non_drought_years, non_drought_gpp_values)
      
      # Compute GPP difference for each drought year and ensure numeric vector is returned
      gpp_differences_local <- unlist(sapply(drought_years, function(year) {
        fitted_gpp <- predict(spline_model, year)$y
        actual_gpp <- gpp_array_reshaped[month, year - 1982 + 1, lat, lon]
        actual_gpp - fitted_gpp
      }))
      
      # Store GPP difference for each drought year into month_diff
      for (i in seq_along(drought_years)) {
        year <- drought_years[i]
        if (i <= length(gpp_differences_local)) {
          month_diff[year - 1982 + 1, lat, lon] <- gpp_differences_local[i]
        }
      }
      
      # Check if lengths match
      if (length(drought_years) != length(gpp_differences_local)) {
        warning("Lengths of drought_years and gpp_differences do not match!")
      }
      
      # Accumulate the total difference matrix
      total_diff[lat, lon] <<- total_diff[lat, lon] + sum(month_diff[, lat, lon], na.rm = TRUE)
      
      # Update gpp_loss_array for this month with the differences from month_diff
      gpp_loss_array[month, lat, lon] <<- gpp_loss_array[month, lat, lon] + sum(month_diff[, lat, lon], na.rm = TRUE)
      
      # Print GPP loss for each pixel and month
      cat("Month:", month, "Lat:", lat, "Lon:", lon, 
          "GPP Loss:", gpp_loss_array[month, lat, lon], "\n")
      
      return(NULL)
    }, error = function(e) {
      return(NULL) 
    })
  })
}

# Close the progress bar
close(pb)


# Output results ----------------------------------------------------------------------

# Define output file name
output_file <- "G:/paper01/TEST/CDE_loss.RData"

# Save array to .RData file
save(gpp_loss_array, file = output_file)

# Print confirmation message
cat("Array has been saved to file:", output_file, "\n")

# Define file path
input_file <- "G:/VPD&SM_CDE(M03)/GIMMS_NIRv/CDE_chara/VPD_reshaped.RData"




