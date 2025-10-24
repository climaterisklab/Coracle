### author: Puja Pande
### date: 06 October 2025
### GMST - marine heatwave regressions

#### load libraries ####
library(terra)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(DescTools)
library(spatstat)


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
      
      # ===== GOODNESS OF FIT METRICS =====
      
      # 1. Deviance-based metrics
      null_dev <- lm_fit$null.deviance
      resid_dev <- lm_fit$deviance
      df_null <- lm_fit$df.null
      df_resid <- lm_fit$df.residual
      
      # Deviance goodness-of-fit p-value
      deviance_pval <- pchisq(resid_dev, df_resid, lower.tail = FALSE)
      
      # Pseudo R-squared (McFadden)
      pseudo_r2 <- 1 - (resid_dev / null_dev)
      #pseudo_r2_mcfadden <- PseudoR2(lm_fit, which = "McFadden")
      #pseudo_r2_nagelkerke <- PseudoR2(lm_fit, which = "Nagelkerke")
      
      # 2. Dispersion parameter
      dispersion <- resid_dev / df_resid
      
      # 3. Pearson chi-square test
      pearson_resid <- residuals(lm_fit, type = "pearson")
      pearson_chisq <- sum(pearson_resid^2)
      pearson_pval <- pchisq(pearson_chisq, df_resid, lower.tail = FALSE)
      
      # 4. AIC and BIC
      #model_aic <- AIC(lm_fit)
      #model_bic <- BIC(lm_fit)
      
      # 5. Log-likelihood
      #log_lik <- logLik(lm_fit)[1]
      
      result_count <- result_count + 1
      
      result_row <- data.frame(
        location_id = loc_id,
        longitude = lon,
        latitude = lat,
        month = month_num,
        month_name = month.name[month_num],
        n_obs = nrow(merged_df),
        mean_dhw = mean(merged_df$mean_dhw, na.rm = TRUE),
        
        # Model coefficients
        alpha = coef_summary[1, 1],
        alpha_se = coef_summary[1, 2],
        beta = coef_summary[2, 1],
        beta_se = coef_summary[2, 2],
        p_value = coef_summary[2, 4],
        
        # Goodness of fit metrics
        #aic = model_aic,
        #bic = model_bic,
        #log_likelihood = log_lik,
        null_deviance = null_dev,
        residual_deviance = resid_dev,
        deviance_pval = deviance_pval,
        pseudo_r2 = pseudo_r2,
        #spatstat_pseudoR2 = spatstat_pseudoR2,
        #pseudo_r2_mcfadden = pseudo_r2_mcfadden, 
        #pseudo_r2_nagelkerke = pseudo_r2_nagelkerke, 
        dispersion = dispersion,
        pearson_chisq = pearson_chisq,
        pearson_pval = pearson_pval,
        df_residual = df_resid
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
# Save results
write.csv(final_poisson_results, "poisson_model_results_with_gof.csv", row.names = FALSE)
final_poisson_results <- read.csv("poisson_model_results_with_gof.csv")



#### dhw, GMT plots ####
# Extract the 10 location_ids corresponding to your site IDs
ten_location_ids <- final_poisson_results %>%
  filter(Site_ID %in% c(8, 9498, 9469, 8521, 6270, 249, 3058, 3567, 7338, 8027)) %>%
  distinct(location_id) %>%
  pull(location_id)

# Create an empty list for storing results
plot_data_list <- list()

for (loc_id in ten_location_ids) {
  
  loc_data <- all_ts_data %>% filter(location_id == loc_id)
  if (nrow(loc_data) == 0) next
  
  lon <- loc_data$lon[1]
  lat <- loc_data$lat[1]
  
  # Aggregate to yearly means
  ts_yearly <- loc_data %>%
    group_by(year) %>%
    summarise(mean_dhw = mean(dhw, na.rm = TRUE), .groups = "drop")
  
  merged_df <- ts_yearly %>%
    inner_join(gmst_data, by = c("year" = "Year")) %>%
    drop_na(mean_dhw, GMT)
  
  if (nrow(merged_df) < 2) next
  
  merged_df <- merged_df %>%
    mutate(
      mean_dhw = pmax(mean_dhw, 0.001),
      location_id = factor(loc_id),
      lon = lon,
      lat = lat
    )
  
  plot_data_list[[as.character(loc_id)]] <- merged_df
}
# All_Bleaching_Events_Data_AllDHW <- read.csv("All_Bleaching_Events_Data_AllDHW.csv")
ten_sites_info <- All_Bleaching_Events_Data_AllDHW %>%
  filter(Site_ID %in% c(8, 9498, 9469, 8521, 6270, 249, 3058, 3567, 7338, 8027)) %>%
  distinct(Site_ID, Ecoregion_Name)

# Merge ecoregion names into your plot data
plot_data_named <- plot_data %>%
  left_join(ten_sites_info, by = "location_id")

# Combine all 10 sites
plot_data <- bind_rows(plot_data_list)

# Verify that location_id exists
str(plot_data$location_id)

# Plot
ggplot(plot_data, aes(x = GMT, y = mean_dhw)) +
  geom_point(color = "steelblue", size = 2) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = TRUE, color = "darkred") +
  facet_wrap(~ location_id, scales = "free_y") +
  labs(
    title = "Relationship between Global Mean Surface Temperature (GMT) and DHW",
    subtitle = "Ten selected reef locations",
    x = "Global Mean Surface Temperature (°C)",
    y = "Degree Heating Weeks (DHW)"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))





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
# 
# pal <- colorNumeric(
#   palette = colorRampPalette(c("#08306B", "#252525", "#99000D"))(256),  # deep navy → charcoal → dark crimson
#   domain = domain,
#   na.color = "transparent"
# )


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




#### looking at model fits ####
# Filter for insignificant points (p-value > 0.05)
insig_data <- final_poisson_results %>% 
  filter(p_value > 0.05)

# Create raster list
raster_list_insig <- list()

# Define global grid (0.25° x 0.25°) based on all coordinate ranges
global_ext <- terra::ext(
  range(final_poisson_results$longitude, na.rm = TRUE),
  range(final_poisson_results$latitude, na.rm = TRUE)
)
world_grid <- terra::rast(global_ext, resolution = 0.25, crs = "EPSG:4326")

# Rasterize monthly insignificant counts
for (m in unique(insig_data$month_name)) {
  month_data <- insig_data %>% filter(month_name == m)
  if (nrow(month_data) == 0) next
  
  points_spat <- terra::vect(
    month_data[, c("longitude", "latitude")],
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  )
  
  # Add a count field (1 per insignificant observation)
  points_spat$insig <- 1
  
  # Rasterize by summing counts of insignificant sites per cell
  insig_raster <- terra::rasterize(points_spat, world_grid, field = "insig", fun = sum)
  raster_list_insig[[m]] <- insig_raster
}

# Combine raster values for color scaling
all_values_insig <- unlist(lapply(raster_list_insig, terra::values))
val_range_insig <- range(all_values_insig, na.rm = TRUE)

# Define color palette (light yellow → dark red)
pal_insig <- colorNumeric(
  palette = colorRampPalette(c("black", "yellow"))(256),
  domain = val_range_insig,
  na.color = "transparent"
)

# Initialize leaflet map (light base map for contrast)
leaf_insig <- leaflet() %>%
  addProviderTiles("CartoDB.Positron")

# Add each month’s raster as a separate selectable layer
for (m in names(raster_list_insig)) {
  leaf_insig <- leaf_insig %>%
    addRasterImage(
      raster_list_insig[[m]],
      colors = pal_insig,
      opacity = 0.8,
      group = m
    )
}

# Add layer control and legend
leaf_insig <- leaf_insig %>%
  addLayersControl(
    baseGroups = names(raster_list_insig),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addLegend(
    position = "bottomright",
    pal = pal_insig,
    values = all_values_insig,
    title = "Count of Insignificant Cells"
  )

# Display map
leaf_insig


# Deviance residuals vs fitted values
plot(fitted(lm_fit), residuals(lm_fit, type = "deviance"),
     xlab = "Fitted values", ylab = "Deviance residuals",
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "red", lty = 2)




#### chi-squared plots ####
#### monthly Pearson χ² maps ####

library(terra)
library(dplyr)
library(leaflet)
library(viridisLite)

# Create a list to store monthly rasters
raster_list_chisq <- list()

# Define global grid (same as before)
global_ext <- terra::ext(
  range(final_poisson_results$longitude, na.rm = TRUE),
  range(final_poisson_results$latitude, na.rm = TRUE)
)
world_grid <- terra::rast(global_ext, resolution = 0.25, crs = "EPSG:4326")

# Rasterize mean Pearson χ² per grid cell for each month
for (m in unique(final_poisson_results$month_name)) {
  month_data <- final_poisson_results %>% filter(month_name == m)
  if (nrow(month_data) == 0) next
  
  points_spat <- terra::vect(
    month_data[, c("longitude", "latitude", "pearson_chisq")],
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  )
  
  chisq_raster <- terra::rasterize(points_spat, world_grid, field = "pearson_chisq", fun = mean)
  raster_list_chisq[[m]] <- chisq_raster
}

# Combine raster values for color scaling
all_values_chisq <- unlist(lapply(raster_list_chisq, terra::values))
val_range_chisq <- range(all_values_chisq, na.rm = TRUE)

# Define a diverging palette (blue → white → red)
pal_chisq <- colorNumeric(
  palette = colorRampPalette(c("darkblue", "white", "darkred"))(256),
  domain = val_range_chisq,
  na.color = "transparent"
)

# Initialize leaflet map
leaf_chisq <- leaflet() %>%
  addProviderTiles("CartoDB.Positron")

# Add monthly raster layers
for (m in names(raster_list_chisq)) {
  leaf_chisq <- leaf_chisq %>%
    addRasterImage(
      raster_list_chisq[[m]],
      colors = pal_chisq,
      opacity = 0.8,
      group = m
    )
}

# Add layer control and legend
leaf_chisq <- leaf_chisq %>%
  addLayersControl(
    baseGroups = names(raster_list_chisq),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addLegend(
    position = "bottomright",
    pal = pal_chisq,
    values = all_values_chisq,
    title = expression("Mean Pearson χ"^2)
  )

# Display the interactive map
leaf_chisq




#### PSEUDO R² PLOTS ####
library(terra)
library(dplyr)
library(leaflet)
library(viridisLite)

#### OPTION 1: McFadden's Pseudo R² Maps ####

# Create a list to store monthly rasters
raster_list_r2 <- list()

# Define global grid (same as before)
global_ext <- terra::ext(
  range(final_poisson_results$longitude, na.rm = TRUE),
  range(final_poisson_results$latitude, na.rm = TRUE)
)

world_grid <- terra::rast(global_ext, resolution = 0.25, crs = "EPSG:4326")

# Rasterize mean Pseudo R² per grid cell for each month
for (m in unique(final_poisson_results$month_name)) {
  month_data <- final_poisson_results %>% filter(month_name == m)
  if (nrow(month_data) == 0) next
  
  points_spat <- terra::vect(
    month_data[, c("longitude", "latitude", "pseudo_r2")],
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  )
  
  r2_raster <- terra::rasterize(points_spat, world_grid, 
                                field = "pseudo_r2", fun = mean)
  raster_list_r2[[m]] <- r2_raster
}

# Combine raster values for color scaling
all_values_r2 <- unlist(lapply(raster_list_r2, terra::values))
val_range_r2 <- range(all_values_r2, na.rm = TRUE)

# Define color palette (white → yellow → red for R²)
# Use a sequential palette since R² is bounded [0,1]
pal_r2 <- colorNumeric(
  palette = colorRampPalette(c("darkblue", "white", "darkred"))(256),
  domain = c(min(val_range_r2), max(val_range_r2, 0.5)),  # Cap at 0.5 or max value
  na.color = "transparent"
)

# Initialize leaflet map
leaf_r2 <- leaflet() %>%
  addProviderTiles("CartoDB.Positron")

# Add monthly raster layers
for (m in names(raster_list_r2)) {
  leaf_r2 <- leaf_r2 %>%
    addRasterImage(
      raster_list_r2[[m]],
      colors = pal_r2,
      opacity = 0.8,
      group = m
    )
}

# Add layer control and legend
leaf_r2 <- leaf_r2 %>%
  addLayersControl(
    baseGroups = names(raster_list_r2),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addLegend(
    position = "bottomright",
    pal = pal_r2,
    values = all_values_r2,
    title = "McFadden's Pseudo R²",
    labFormat = labelFormat(digits = 3)
  )

# Display the interactive map
leaf_r2




