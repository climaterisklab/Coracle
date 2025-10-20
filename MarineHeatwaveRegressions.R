### author: Puja Pande
### date: 06 October 2025
### GMST - marine heatwave regressions

#### load libraries ####
library(terra)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

#### load gmst and DHW data ####
gmst_data <- read.csv("~/Library/CloudStorage/Dropbox/Coracle/ERA5_GMT.csv") %>% select(-X) %>%
  mutate(date = as.Date(paste0(Year, "-06-15")))
# all the data 
dhw_data <- read.csv("Panel_Data_AllDHW.csv") %>% select(-X)

#### quick plot ####
plot(gmst_data$Year, gmst_data$GMT, type = "l", main = "Observed GMST")

#### unique location and rasterbrick ####
# Get unique lat/long pairs
unique_locations <- dhw_data %>%
  distinct(Longitude_Degrees, Latitude_Degrees, .keep_all = TRUE) %>%
  select(Site_ID, Longitude_Degrees, Latitude_Degrees) %>%
  rename(lon = Longitude_Degrees, lat = Latitude_Degrees)

# Read the rasterbrick
rasterbrick_dhw <- rast("rasterbrick_dhw.tif")

#### pre-extracting all the data ####
# Extract all points at once (much faster than looping)
coords_matrix <- as.matrix(unique_locations[, c("lon", "lat")])
all_extracted <- terra::extract(rasterbrick_dhw, coords_matrix)

# Create a list of time series for each location
raster_dates <- as.Date(time(rasterbrick_dhw))
ts_list <- list()

for (i in seq_len(nrow(unique_locations))) {
  vals <- as.numeric(all_extracted[i, -1])  # Drop ID column
  
  ts_list[[i]] <- data.frame(
    location_id = i,
    lon = unique_locations$lon[i],
    lat = unique_locations$lat[i],
    date = raster_dates[seq_along(vals)],
    dhw = vals
  )
}

#### Month/Year Aggregation ####
# Create storage for results
model_results <- list()
model_list <- list()
error_log <- data.frame()

# Bind all ts data once for faster processing
all_ts_data <- bind_rows(ts_list)

# Add month and year columns once
all_ts_data <- all_ts_data %>%
  mutate(
    month = month(date),
    year = year(date)
  )

#### running poisson regression model for every pixel ####
cat("Starting model fitting...\n")
start_time <- Sys.time()

model_results <- list()
result_count <- 0

# Loop over each location
for (loc_id in seq_len(nrow(unique_locations))) {
  
  if (loc_id %% 1000 == 0) {
    cat("Processing location", loc_id, "of", nrow(unique_locations), "\n")
  }
  
  # Get this location's data
  loc_data <- all_ts_data %>% filter(location_id == loc_id)
  
  if (nrow(loc_data) == 0) next
  
  lon <- loc_data$lon[1]
  lat <- loc_data$lat[1]
  
  # Loop over months
  for (month_num in 1:12) {
    
    ts_month <- loc_data %>%
      filter(month == month_num)
    
    if (nrow(ts_month) == 0) next
    
    # Aggregate to yearly
    ts_yearly <- ts_month %>%
      group_by(year) %>%
      summarise(mean_dhw = mean(dhw, na.rm = TRUE), .groups = 'drop')
    
    # Join with GMST
    merged_df <- ts_yearly %>%
      inner_join(gmst_data, by = c("year" = "Year"))
    
    if (nrow(merged_df) < 2) next
    
    # Remove NAs
    merged_df <- merged_df %>% drop_na(mean_dhw, GMT)
    
    if (nrow(merged_df) < 2) next
    
    # Handle zeros
    merged_df$mean_dhw <- pmax(merged_df$mean_dhw, 0.001)
    
    # Fit model with error handling
    tryCatch({
      lm_fit <- glm(mean_dhw ~ GMT, family = "poisson", data = merged_df)
      
      coef_summary <- summary(lm_fit)$coefficients
      
      result_count <- result_count + 1
      
      result_row <- data.frame(
        location_id = loc_id,
        longitude = lon,
        latitude = lat,
        month = month_num,
        month_name = month.name[month_num],
        n_obs = nrow(merged_df),
        mean_dhw = mean(merged_df$mean_dhw, na.rm = TRUE),
        alpha = coef_summary[1, 1],        # intercept (constant term)
        alpha_se = coef_summary[1, 2],     # intercept standard error
        beta = coef_summary[2, 1],         # slope for GMT
        beta_se = coef_summary[2, 2],      # slope SE
        p_value = coef_summary[2, 4],      # p-value for GMT
        aic = AIC(lm_fit)
      )
      
      model_results[[result_count]] <- result_row
      
    }, error = function(e) {
      # Silently skip failed models
    })
  }
}

end_time <- Sys.time()
cat("Model fitting complete in", 
    round(as.numeric(difftime(end_time, start_time, units = "mins")), 2), 
    "minutes\n")

#### combine results ####
final_poisson_results <- bind_rows(model_results)

dhw_unique <- dhw_data[, 1:3] %>%
  group_by(Longitude_Degrees, Latitude_Degrees) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

final_poisson_results <- left_join(
  final_poisson_results,
  dhw_unique[, 1:3],
  join_by(longitude == Longitude_Degrees, latitude == Latitude_Degrees)
)


# Save results to CSV
write.csv(final_poisson_results, "GMST_DHW_Poisson_Results.csv", row.names = FALSE)
final_poisson_results <- read.csv("GMST_DHW_Poisson_Results.csv")



#### quick plots ####
## Histogram of beta coefficients
hist(final_poisson_results$beta, breaks = 50, main = "Histogram of Beta Coefficients", xlab = "Beta Coefficient")
## Scatter plot of beta vs mean DHW
plot(final_poisson_results$mean_dhw, final_poisson_results$beta, 
     main = "Beta vs Mean DHW", xlab = "Mean DHW", ylab = "Beta Coefficient")
## Boxplot of beta by month
boxplot(beta ~ month_name, data = final_poisson_results, 
        main = "Beta Coefficients by Month", xlab = "Month", ylab = "Beta Coefficient", 
        las = 2)

## more important plots 
# Beta coefficients on world map
### Beta Coefficients on Global Map ----------------------------------------

# Load required libraries
library(terra)
library(dplyr)
library(leaflet)
library(viridisLite)

# Create raster list
raster_list <- list()

# Define global grid (0.25° x 0.25°) based on all coordinate ranges
global_ext <- terra::ext(
  range(final_poisson_results$longitude, na.rm = TRUE),
  range(final_poisson_results$latitude, na.rm = TRUE)
)
world_grid <- terra::rast(global_ext, resolution = 0.25, crs = "EPSG:4326")

# Rasterize monthly beta values
for (m in unique(final_poisson_results$month_name)) {
  month_data <- final_poisson_results %>% filter(month_name == m)
  if (nrow(month_data) == 0) next
  
  points_spat <- terra::vect(
    month_data[, c("longitude", "latitude", "beta")],
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  )
  
  beta_raster <- terra::rasterize(points_spat, world_grid, field = "beta", fun = mean)
  raster_list[[m]] <- beta_raster
}

# Combine raster values for color scaling
all_values <- unlist(lapply(raster_list, terra::values))
val_range <- range(all_values, na.rm = TRUE)

# Ensure color scale is symmetric around 0
max_abs <- max(abs(val_range), na.rm = TRUE)
domain <- c(-max_abs, max_abs)

# Define diverging color palette (dark blue → neutral → dark red)
pal <- colorNumeric(
  palette = colorRampPalette(c("darkblue", "white", "darkred"))(256),  # neutral dark gray for 0
  domain = domain,
  na.color = "transparent"
)

pal <- colorNumeric(
  palette = colorRampPalette(c("#08306B", "#252525", "#99000D"))(256),  # deep navy → charcoal → dark crimson
  domain = domain,
  na.color = "transparent"
)


# Initialize leaflet map with dark basemap
leaf <- leaflet() %>%
  #addProviderTiles("CartoDB.DarkMatterNoLabels") %>%
  addProviderTiles("CartoDB.Positron")

# Add each month’s raster as a separate selectable layer
for (m in names(raster_list)) {
  leaf <- leaf %>%
    addRasterImage(
      raster_list[[m]],
      colors = pal,
      opacity = 0.8,
      group = m
    )
}

# Add layer control and legend
leaf <- leaf %>%
  addLayersControl(
    baseGroups = names(raster_list),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addLegend(
    position = "bottomright",
    pal = pal,
    values = all_values,
    title = expression(beta)
  )

# Display map
leaf




