## author: Puja Pande
## date: 15 Jan 2026
## description: addition of MPA and MPA interaction 

#### Adding MPA 
library(sf)

# Read all files
poly_files <- c(
  "~/Library/CloudStorage/Dropbox/Coracle/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_0/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp",
  "~/Library/CloudStorage/Dropbox/Coracle/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_1/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp",
  "~/Library/CloudStorage/Dropbox/Coracle/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_2/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp"
)

point_files <- c(
  "~/Library/CloudStorage/Dropbox/Coracle/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_0/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp",
  "~/Library/CloudStorage/Dropbox/Coracle/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_1/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp",
  "~/Library/CloudStorage/Dropbox/Coracle/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_2/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp"
)

# Merge polygons
poly_list <- lapply(poly_files, st_read)
merged_poly <- do.call(rbind, poly_list)

# Merge points
point_list <- lapply(point_files, st_read)
merged_points <- do.call(rbind, point_list)

# Get all unique column names
all_cols <- union(names(merged_poly), names(merged_points))

# Add missing columns to each dataset
for (col in all_cols) {
  if (!col %in% names(merged_poly)) {
    merged_poly[[col]] <- NA
  }
  if (!col %in% names(merged_points)) {
    merged_points[[col]] <- NA
  }
}

# Reorder columns to match
merged_poly <- merged_poly[, all_cols]
merged_points <- merged_points[, all_cols]

# Now combine
all_marine <- rbind(merged_poly, merged_points)

#st_write(all_marine, "merged_marine_all.gpkg", delete_dsn = TRUE)



#### adding mpa to data 
merged <- read.csv("Final_Relevant_scripts/Merged_Mermaid_Panel_Bleaching_Data.csv")

# Convert your data to sf object (points)
merged_sf <- st_as_sf(merged, 
                      coords = c("Longitude_Degrees", "Latitude_Degrees"),  # Adjust column names
                      crs = 4326)  # WGS84

# Fix invalid geometries in MPA data
all_marine <- st_make_valid(all_marine)

# If still fails, switch off spherical geometry (s2)
sf_use_s2(FALSE)
merged_with_mpa <- st_join(merged_sf, all_marine, left = TRUE)
sf_use_s2(TRUE)  # Turn back on after

# Add binary MPA indicator
merged_with_mpa$in_mpa <- !is.na(merged_with_mpa$SITE_ID)

# Convert back to dataframe
merged_final <- st_drop_geometry(merged_with_mpa)

# Check results
table(merged_final$in_mpa)

#### running models with MPA interaction
library(fixest)

#### model ####
model_negative_binomial <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                    data = merged_final)
model_negative_binomial_mpa <- fenegbin(Percent_Bleached ~ log1p(dhw) + log1p(dhw)*in_mpa | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                        data = merged_final)
summary(model_negative_binomial)
summary(model_negative_binomial_mpa)
## mpa effect not significant
