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

writeRaster(rasterbrick_dhw, filename = "rasterbrick_dhw.tif")

#### reading in the panel data to extract DHW values or closest non-NA DHW values ####
## doing this for the entire GBCD dataset and the all_bleaching_events dataset
full_panel_data <- read.csv("~/Library/CloudStorage/Dropbox/Coracle/Cleaned GCBD Datasets/GBCD_full_data_cleaned.csv")
full_panel_data <- read.csv("~/Library/CloudStorage/Dropbox/Coracle/Cleaned GCBD Datasets/All_Bleaching_Events_Data.csv")

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
write.csv(full_panel_data, "All_Bleaching_Events_Data_AllDHW.csv")




## doing this for mermaid data ##
# =============================================================================
# PART 1: Load original rasterbrick (1985-2023)
# =============================================================================

cat("=============================================================================\n")
cat("PART 1: Extracting DHW from original data (1985-2023)\n")
cat("=============================================================================\n\n")

# Load Mermaid data
mermaid_data <- readRDS("MermaidBleachingData.RDS")
cat("Mermaid data loaded:", nrow(mermaid_data), "observations\n")
cat("Date range:", min(mermaid_data$year), "-", max(mermaid_data$year), "\n\n")

# Create year-month column
mermaid_data$ym <- sprintf("%04d-%02d", mermaid_data$year, mermaid_data$month)

# Load original rasterbrick (UPDATE THIS PATH to your original data)
cat("Loading original rasterbrick...\n")
rasterbrick_dhw_original <- rast("rasterbrick_dhw.tif")  # Or wherever your original file is

cat("Original rasterbrick loaded:\n")
cat("  Layers:", nlyr(rasterbrick_dhw_original), "\n")
cat("  Layer names (first 5):", head(names(rasterbrick_dhw_original), 5), "\n\n")

# Extract dates from original rasterbrick
# Assuming names are like "DHW04-1985", "DHW05-1985", etc.
# If you know the original dates you used
start_date <- as.Date("1985-04-01")
end_date   <- as.Date("2023-10-01")
raster_dates_original <- seq(start_date, end_date, by = "month")
raster_dhw_ym_original <- format(raster_dates_original, "%Y-%m")
cat("Original DHW date range:", range(raster_dhw_ym_original), "\n\n")

# Map each row's date to layer index
idx_original <- match(mermaid_data$ym, raster_dhw_ym_original)

# Check how many observations have matching data
cat("Observations matching original DHW data:", sum(!is.na(idx_original)), "/", nrow(mermaid_data), "\n\n")

# Extract DHW values from original data
cat("Extracting DHW values from original data...\n")
mermaid_data$dhw_original <- pbsapply(1:nrow(mermaid_data), function(i) {
  layer_idx <- idx_original[i]
  
  if (is.na(layer_idx)) return(NA)
  
  lon <- mermaid_data$longitude[i]
  lat <- mermaid_data$latitude[i]
  
  if (is.na(lon) | is.na(lat)) return(NA)
  
  point <- terra::vect(cbind(lon, lat), crs = terra::crs(rasterbrick_dhw_original))
  val <- terra::extract(rasterbrick_dhw_original[[layer_idx]], point)[1, 2]
  
  # If NA, use 10km buffer
  if (is.na(val)) {
    buff <- terra::buffer(point, width = 10000)
    vals <- terra::extract(rasterbrick_dhw_original[[layer_idx]], buff)[, 2]
    non_na_vals <- vals[!is.na(vals)]
    if (length(non_na_vals) > 0) {
      val <- non_na_vals[1]
    }
  }
  
  return(val)
})

cat("\nPart 1 complete!\n")
cat("  Observations with DHW:", sum(!is.na(mermaid_data$dhw_original)), "\n")
cat("  DHW range:", range(mermaid_data$dhw_original, na.rm = TRUE), "°C-weeks\n\n")

# =============================================================================
# PART 2: Load new calculated files (2023-2025)
# =============================================================================

cat("=============================================================================\n")
cat("PART 2: Extracting DHW from new calculated data (2023-2025)\n")
cat("=============================================================================\n\n")

root_dir <- "~/Library/CloudStorage/Dropbox/Coracle/NOAA DHW Data/DHW 2023 - 2025_Calculated"
rasters_mermaid <- list.files(
  root_dir,
  pattern = "ct5km_dhw-max_v3.1_.*\\.nc$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(rasters_mermaid) == 0) {
  cat("No new DHW files found. Skipping Part 2.\n\n")
  mermaid_data$dhw_new <- NA
} else {
  cat("Found", length(rasters_mermaid), "new DHW files\n")
  
  rasters_mermaid <- sort(rasters_mermaid)
  rasterbrick_dhw_new <- rast(rasters_mermaid, subds = "degree_heating_week")
  
  # Extract dates from filenames
  extract_date_from_filename <- function(filename) {
    date_str <- gsub(".*_(\\d{6})\\.nc$", "\\1", basename(filename))
    year <- as.numeric(substr(date_str, 1, 4))
    month <- as.numeric(substr(date_str, 5, 6))
    return(as.Date(paste(year, month, "01", sep = "-")))
  }
  
  raster_dates_new <- sapply(rasters_mermaid, extract_date_from_filename)
  raster_dates_new <- as.Date(raster_dates_new, origin = "1970-01-01")
  names(rasterbrick_dhw_new) <- paste0("DHW", format(raster_dates_new, "%m-%Y"))
  raster_dhw_ym_new <- format(raster_dates_new, "%Y-%m")
  
  cat("New DHW date range:", range(raster_dhw_ym_new), "\n\n")
  
  # Map each row's date to layer index
  idx_new <- match(mermaid_data$ym, raster_dhw_ym_new)
  
  cat("Observations matching new DHW data:", sum(!is.na(idx_new)), "/", nrow(mermaid_data), "\n\n")
  
  if (sum(!is.na(idx_new)) > 0) {
    cat("Extracting DHW values from new data...\n")
    mermaid_data$dhw_new <- pbsapply(1:nrow(mermaid_data), function(i) {
      layer_idx <- idx_new[i]
      
      if (is.na(layer_idx)) return(NA)
      
      lon <- mermaid_data$longitude[i]
      lat <- mermaid_data$latitude[i]
      
      if (is.na(lon) | is.na(lat)) return(NA)
      
      point <- terra::vect(cbind(lon, lat), crs = terra::crs(rasterbrick_dhw_new))
      val <- terra::extract(rasterbrick_dhw_new[[layer_idx]], point)[1, 2]
      
      # If NA, use 10km buffer
      if (is.na(val)) {
        buff <- terra::buffer(point, width = 10000)
        vals <- terra::extract(rasterbrick_dhw_new[[layer_idx]], buff)[, 2]
        non_na_vals <- vals[!is.na(vals)]
        if (length(non_na_vals) > 0) {
          val <- non_na_vals[1]
        }
      }
      
      return(val)
    })
    
    cat("\nPart 2 complete!\n")
    cat("  Observations with DHW:", sum(!is.na(mermaid_data$dhw_new)), "\n")
    cat("  DHW range:", range(mermaid_data$dhw_new, na.rm = TRUE), "°C-weeks\n\n")
  } else {
    cat("No observations match new DHW data dates.\n\n")
    mermaid_data$dhw_new <- NA
  }
}

# =============================================================================
# PART 3: Combine DHW values
# =============================================================================

cat("=============================================================================\n")
cat("PART 3: Combining DHW values\n")
cat("=============================================================================\n\n")

# Combine: use original data where available, fill in with new data where needed
mermaid_data$dhw <- ifelse(!is.na(mermaid_data$dhw_original), 
                           mermaid_data$dhw_original, 
                           mermaid_data$dhw_new)
