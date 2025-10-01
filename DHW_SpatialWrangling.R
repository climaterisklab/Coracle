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
## saving this panel data
write.csv(full_panel_data, "Panel_Data_DHW.csv")


