### author: Puja Pande
### date: 01 October 2025
### description: this script does all the spatial wrangling required for the degree heating weeks attribute from NOAA

#### load libraries ####
library(terra)
library(lubridate)
library(dplyr)
library(pbapply)


#### creating the rasterbrick ####
root_dir <- "~/Library/CloudStorage/Dropbox/Coracle/NOAA DHW Data/"
rasters <- list.files(
  root_dir,
  pattern = "ct5km_dhw-max_v3.1_.*\\.nc$",
  full.names = TRUE,
  recursive = TRUE
)
rasters <- sort(rasters)
rasterbrick_dhw <- rast(rasters, subds = "degree_heating_week")

# check
rasterbrick_dhw
nlyr(rasterbrick_dhw)
# Define start and end dates
start_date <- as.Date("1985-04-01")
end_date   <- as.Date("2023-10-01")
# Generate a monthly sequence
date_seq <- seq(start_date, end_date, by = "month")
# Create names like DHW4-1985, DHW5-1985, ... , DHW10-2023
names(rasterbrick_dhw) <- paste0("DHW", format(date_seq, "%m-%Y"))


#### reading in the panel data to extract DHW values or closest non-NA DHW values ####
full_panel_data <- read.csv("~/Library/CloudStorage/Dropbox/Coracle/Cleaned GCBD Datasets/GBCD_full_data_cleaned.csv")

#### extracting DHW values for each row in the panel data ####
# Add year-month column
full_panel_data$ym <- format(
  ymd(paste(full_panel_data$Date_Year, full_panel_data$Date_Month, "01", sep = "-")),
  "%Y-%m"
)

# Raster dates
raster_dhw_dates <- as.Date(time(rasterbrick_dhw))
raster_dhw_ym <- format(raster_dhw_dates, "%Y-%m")

# Map each row's date to layer index
idx <- match(full_panel_data$ym, raster_dhw_ym)

# Extract DHW values directly at sites (with NA handling)
full_panel_data$dhw <- pbsapply(1:nrow(full_panel_data), function(i) {
  layer_idx <- idx[i]
  if (is.na(layer_idx)) return(NA)
  
  # Create point for this site
  point <- terra::vect(cbind(full_panel_data$Longitude[i], full_panel_data$Latitude[i]), 
                       crs = terra::crs(rasterbrick_dhw))
  
  # Extract value at point
  val <- terra::extract(rasterbrick_dhw[[layer_idx]], point)[1, 2]  # [row, column after ID]
  
  # If NA, extract from buffer and take nearest non-NA
  if (is.na(val)) {
    buff <- terra::buffer(point, width = 10000)  # 10km buffer
    vals <- terra::extract(rasterbrick_dhw[[layer_idx]], buff)[, 2]  # All values except ID column
    non_na_vals <- vals[!is.na(vals)]
    if (length(non_na_vals) > 0) {
      val <- non_na_vals[1]
    }
  }
  
  return(val)
})

#### creating DHW lags from 0 - 5 months ####
## efficient lag function - batch extraction
# 1. Extract full time series (with NAs)
# 1. Extract full time series for all sites in one call
all_vals <- terra::extract(
  rasterbrick_dhw,
  full_panel_data[, c("Longitude_Degrees", "Latitude_Degrees")]
)

# Drop the ID column
all_vals <- as.matrix(all_vals[, -1])

# 2. Create lagged variables with progress bar
lags <- -1:5
lagged_df <- pblapply(lags, function(L) {
  shifted_idx <- idx - L
  valid <- !is.na(shifted_idx) & shifted_idx > 0 & shifted_idx <= ncol(all_vals)
  vals <- rep(NA_real_, length(idx))
  vals[valid] <- all_vals[cbind(seq_along(idx)[valid], shifted_idx[valid])]
  vals
})

# 3. Attach lagged columns to full_panel_data and fill NAs using buffer
for (j in seq_along(lags)) {
  col_name <- paste0("dhw_lag", lags[j])
  full_panel_data[[col_name]] <- lagged_df[[j]]
  
  # Fill remaining NAs with buffer extraction
  na_indices <- which(is.na(full_panel_data[[col_name]]))
  if (length(na_indices) > 0) {
    full_panel_data[[col_name]][na_indices] <- pbsapply(na_indices, function(i) {
      layer_idx <- idx[i] - lags[j]
      if (is.na(layer_idx) || layer_idx < 1 || layer_idx > nlyr(rasterbrick_dhw)) return(NA)
      
      point <- terra::vect(
        cbind(full_panel_data$Longitude_Degrees[i], full_panel_data$Latitude_Degrees[i]),
        crs = terra::crs(rasterbrick_dhw)
      )
      
      buff <- terra::buffer(point, width = 10000)
      vals <- terra::extract(rasterbrick_dhw[[layer_idx]], buff)[, 2]
      non_na_vals <- vals[!is.na(vals)]
      if (length(non_na_vals) > 0) return(non_na_vals[1])
      return(NA)
    })
  }
}


## saving this panel data
write.csv(full_panel_data, "Panel_Data_AllDHW.csv")

