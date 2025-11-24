## author: Puja Pande
## date: 21 November 2025
## description: joining mermaid data and adding dhw values

#### load libraries ####
library(tidyr)
library(dplyr)
library(terra)
library(pbapply)

#### load data ####
mermaid_data <- readRDS("MermaidBleachingData.RDS")
rasterbrick_dhw <- rast("rasterbrick_dhw.tif")

#### extracting DHW values or closest non-NA DHW values ####
# Raster dates
raster_dhw_dates <- as.Date(time(rasterbrick_dhw))
raster_dhw_ym <- format(raster_dhw_dates, "%Y-%m")

# Create year-month column in mermaid_data
mermaid_data$ym <- format(mermaid_data$date, "%Y-%m")

# Map each row's date to layer index
idx <- match(mermaid_data$ym, raster_dhw_ym)

# Extract DHW values directly at sites (with NA handling)
mermaid_data$dhw <- pbsapply(1:nrow(mermaid_data), function(i) {
  layer_idx <- idx[i]
  if (is.na(layer_idx)) return(NA)
  
  # Create point for this site
  point <- terra::vect(cbind(mermaid_data$longitude[i], mermaid_data$latitude[i]), 
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

#### doing lags ####
#### creating DHW lags from 0 - 5 months ####
## efficient lag function - batch extraction
# 1. Extract full time series for all sites in one call
all_vals <- terra::extract(
  rasterbrick_dhw,
  mermaid_data[, c("longitude", "latitude")]
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

# 3. Attach lagged columns to mermaid_data and fill NAs using buffer
for (j in seq_along(lags)) {
  col_name <- paste0("dhw_lag", lags[j])
  mermaid_data[[col_name]] <- lagged_df[[j]]
  
  # Fill remaining NAs with buffer extraction
  na_indices <- which(is.na(mermaid_data[[col_name]]))
  if (length(na_indices) > 0) {
    mermaid_data[[col_name]][na_indices] <- pbsapply(na_indices, function(i) {
      layer_idx <- idx[i] - lags[j]
      if (is.na(layer_idx) || layer_idx < 1 || layer_idx > nlyr(rasterbrick_dhw)) return(NA)
      
      point <- terra::vect(
        cbind(mermaid_data$longitude[i], mermaid_data$latitude[i]),
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

mermaid_data <- mermaid_data %>% rename(`dhw_lag.1` = `dhw_lag-1`)

#### merging with GCBD dataset ####
gcbd_data <- read.csv("All_Bleaching_Events_Data_AllDHW.csv")

head(mermaid_data)
head(gcbd_data)


# Create a lookup table of unique site locations with their ecoregion and realm
location_lookup <- gcbd_data %>%
  select(Latitude_Degrees, Longitude_Degrees, Ecoregion_Name, Realm_Name) %>%
  distinct()

# Join mermaid_data with location_lookup based on coordinates
# Using a spatial join approach with tolerance for slight coordinate differences
mermaid_data_enriched <- mermaid_data %>%
  left_join(
    location_lookup,
    by = c("latitude" = "Latitude_Degrees", "longitude" = "Longitude_Degrees")
  )

# If exact matches don't work well, use a nearest neighbor approach:
library(geosphere)

# Alternative: Find nearest gcbd location for each mermaid site
mermaid_data_enriched <- mermaid_data %>%
  rowwise() %>%
  mutate(
    nearest_idx = which.min(
      distHaversine(
        cbind(longitude, latitude),
        cbind(location_lookup$Longitude_Degrees, location_lookup$Latitude_Degrees)
      )
    ),
    Ecoregion_Name = location_lookup$Ecoregion_Name[nearest_idx],
    Realm_Name = location_lookup$Realm_Name[nearest_idx]
  ) %>%
  select(-nearest_idx) %>%
  ungroup()


















