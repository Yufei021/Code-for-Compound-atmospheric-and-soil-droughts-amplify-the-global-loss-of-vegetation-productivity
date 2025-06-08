library(ncdf4)
library(raster)
library(reticulate)

# CDE_loss.RData contains the monthly GPP_loss data (from January to December) over the study period
load("G:/paper01/VPD&SM_CDE(M03)/GIMMS_NIRv/GPP_loss_monthly/CDE_loss.RData")
load("G:/paper01/VPD&SM_CDE(M03)/GIMMS_NIRv/CDE_chara/CDE_reshaped.RData")

drought_events_reshaped00 <- drought_events_reshaped02

# Read GPP data ------------------------------------------
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
      
      # Check if the raster dimensions are 360 x 720
      if (nrow(gpp_raster) == 360 && ncol(gpp_raster) == 720) {
        # Extract all rows and columns
        gpp_data <- as.matrix(gpp_raster)
      } else {
        stop("GPP data dimensions do not match expected size (360x720)")
      }
      
      # Store the monthly data into original_gpp_array (3D array: [total_months, lat, lon])
      original_gpp_array[month_counter, , ] <- gpp_data
      
      # Increment the month counter
      month_counter <- month_counter + 1
    } else {
      warning(sprintf("File does not exist: %s", gpp_file))
    }
  }
}
# Output the dimensions of the 3D array for confirmation
dim(original_gpp_array) 

gpp_array <- original_gpp_array[, 1:360, 1:720]
# Output the dimensions of the extracted data for confirmation
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

# Calculate recovery time ---------------------------------------------

library(tcltk)

# Create a progress bar
pb <- tkProgressBar(title = "Progress", min = 0, max = 360 * 720, width = 300)

# Create an array to store recovery time
recovery_time <- array(NA, dim = dim(drought_events_reshaped02))

# Create a 2D matrix to store average recovery time
average_recovery_time <- matrix(NA, nrow = 360, ncol = 720)
total_recovery_time <- matrix(NA, nrow = 360, ncol = 720)
# Initialize a 2D matrix to store the count of invalid recovery events
invalid_recovery_count <- matrix(0, nrow = 360, ncol = 720)

# Loop through each pixel
for (lat in 1:360) {
  for (lon in 1:720) {
    drought_sequence <- drought_events_reshaped02[, , lat, lon]
    gpp_sequence <- gpp_array_reshaped[, , lat, lon]
    
    # Extract drought event indices
    drought_indices <- which(drought_sequence == 1, arr.ind = TRUE)
    
    # Add year information
    if (ncol(drought_indices) == 1) {
      years <- rep(1:dim(drought_sequence)[2], each = 12)
      drought_indices <- cbind(drought_indices, years[drought_indices[, 1]])
    }
    
    if (nrow(drought_indices) > 0) {
      
      # Fix grouping logic: consider temporal continuity
      drought_events <- split(
        drought_indices,
        cumsum(
          c(TRUE, drought_indices[-1, 2] != drought_indices[-nrow(drought_indices), 2] |  # Across years
              drought_indices[-1, 1] != drought_indices[-nrow(drought_indices), 1] + 1)    # Non-consecutive months
        )
      )
      
      total_recovery_time_lat_lon <- 0  # Total recovery time for current pixel
      total_drought_events_lat_lon <- 0  # Total number of drought events for current pixel
      
      for (event in drought_events) {
        if (is.null(event) || nrow(as.matrix(event)) < 1) next
        
        # Force matrix structure
        event <- matrix(event, ncol = 2, byrow = FALSE)
        
        if (ncol(event) < 2 || nrow(event) < 1) next
        
        # Get the last month and year of the drought event
        end_month <- event[nrow(event), 1]
        end_year <- event[nrow(event), 2]
        
        # Calculate the index of the ending month
        end_index <- (end_year - 1) * 12 + end_month
        
        # Search for the month when GPP returns to normal level
        recovery_found <- FALSE
        year <- end_year
        month <- end_month
        while (!recovery_found) {
          # Move to the next month
          month <- month + 1
          if (month > 12) {
            month <- 1
            year <- year + 1
          }
          # Check if index exceeds data range
          if (year > ncol(gpp_sequence) || month > nrow(gpp_sequence)) {
            break
          }
          
          # Get the normal GPP value
          normal_gpp <- gpp_loss_array01[month, year, lat, lon]
          
          # Check if GPP exceeds the normal value
          if (!is.na(gpp_sequence[month, year]) && !is.na(normal_gpp) && 
              gpp_sequence[month, year] > normal_gpp) {
            recovery_found <- TRUE
            
            # Check whether it remains above for at least 3/2 months (Yao's criterion)
            valid_recovery <- TRUE
            for (i in 1:2) {
              # Calculate future month and year
              future_month <- (month + i - 1) %% 12 + 1
              future_year <- year + (month + i - 1) %/% 12
              
              # Check validity of future dates
              if (future_month > nrow(gpp_sequence) || future_year > ncol(gpp_sequence) ||
                  is.na(gpp_sequence[future_month, future_year]) ||
                  is.na(gpp_loss_array01[future_month, future_year, lat, lon]) ||
                  gpp_sequence[future_month, future_year] <= gpp_loss_array01[future_month, future_year, lat, lon]) {
                valid_recovery <- FALSE
                break
              }
            }
            
            # If it satisfies the 2-month criterion, record recovery time
            if (valid_recovery) {
              if (is.na(recovery_time[month, year, lat, lon])) {
                recovery_time[month, year, lat, lon] <- ((year - 1) * 12 + month) - end_index
                
                # Print recovery time info
                print(sprintf("Lat: %d, Lon: %d, Recovery Time Found: %d months (valid for 2 months)",
                              lat, lon,
                              recovery_time[month, year, lat, lon]))
              }
              break
            }
          }
          if (recovery_found) break
        }
        
        # If recovery was not found, mark as NA
        if (!recovery_found) {
          recovery_time[end_month, end_year, lat, lon] <- NA
        }
        
        # # Count drought events
        # total_drought_events_lat_lon <- total_drought_events_lat_lon + 1
      }
      
      # Calculate total recovery time for this pixel
      valid_recovery_times <- recovery_time[, , lat, lon]
      
      # Filter valid recovery times (≤ 12 months)
      valid_recovery_times <- valid_recovery_times[!is.na(valid_recovery_times) & valid_recovery_times <= 12]
      
      # Filter invalid recovery times (> 12 months)
      invalid_recovery_times <- recovery_time[, , lat, lon]
      invalid_recovery_times <- invalid_recovery_times[!is.na(invalid_recovery_times) & invalid_recovery_times > 12]
      # Record the number of invalid recoveries in the 2D array
      invalid_recovery_count[lat, lon] <- length(invalid_recovery_times)
      
      # Count valid drought events
      total_valid_drought_events <- length(valid_recovery_times)
      
      # Sum recovery time
      total_recovery_time_lat_lon <- sum(valid_recovery_times, na.rm = TRUE)
      
      # Calculate average recovery time per drought event
      if (total_valid_drought_events > 0) {
        average_recovery_time_lat_lon <- total_recovery_time_lat_lon / total_valid_drought_events
      } else {
        average_recovery_time_lat_lon <- NA  # Set NA if no valid drought events
      }
      
      # Store the pixel's average recovery time
      total_recovery_time[lat, lon] <- average_recovery_time_lat_lon
    }
    
    # Update progress bar
    progress <- (lat - 1) * 720 + lon
    setTkProgressBar(pb, progress, label = sprintf("Progress: %.2f%%", progress / (360 * 720) * 100))
  }
}

# Close the progress bar
close(pb)

# Plotting --------------------------------------------

# # Check the number of regions marked as 1 in drought_events_reshaped02
# sum(drought_events_reshaped02 == 1, na.rm = TRUE)

# Suppose we want to plot the GPP data for the first month
event_for_month <- total_recovery_time  # Extract data for the first layer
# Remove values equal to 0 by replacing them with NA
event_for_month[event_for_month == 0] <- NA
# event_for_month_matrix <- matrix(event_for_month, nrow = 180, ncol = 360)  # Adjust dimensions as needed

# Convert to a raster object
event_raster <- raster(event_for_month)

# Set spatial resolution of the raster (you may need to adjust origin and intervals)
extent(event_raster) <- c(-180.0499999999, 179.950000001, -90, 90.05)  # Set longitude and latitude range

# Set coordinate reference system to match the map
crs(event_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

# Plot the raster data
plot(event_raster, main = "GPP for Month 1", xlab = "Longitude", ylab = "Latitude")

output_file <- "G:/VPD&SM_CDE(M03)/GIMMS_NIRv/recovery_time/12/v2_Chai/CDE_invalid02.tif"  # Specify output file path
writeRaster(event_raster, filename = output_file, format = "GTiff", overwrite = TRUE)




