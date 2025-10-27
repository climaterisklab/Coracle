### author: Puja Pande
### date: 21 October 2025
### description: EEZ intersections with shallow water tropical reefs

# ============================================================================
# Download and merge reef mask tiles from Google Drive (R version)
# Requires: install.packages(c("googledrive", "terra", "sf"))
# ============================================================================

library(googledrive)
library(terra)

# ============================================================================
# CONFIGURATION
# ============================================================================

FILE_PREFIX <- "global_reef_mask_30m_styled"  # Files starting with this
DOWNLOAD_FOLDER <- "./reef_tiles"             # Local download folder
OUTPUT_FILE <- "merged_reef_mask.tif"         # Final merged output
GDRIVE_FOLDER <- "GEE_Exports"                # Google Drive folder (optional)

# ============================================================================
# OPTION 1: DOWNLOAD FROM GOOGLE DRIVE AND MERGE
# ============================================================================

download_and_merge_from_gdrive <- function() {
  
  # Authenticate with Google Drive (will open browser first time)
  cat("Authenticating with Google Drive...\n")
  drive_auth()
  
  # Create download folder if it doesn't exist
  if (!dir.exists(DOWNLOAD_FOLDER)) {
    dir.create(DOWNLOAD_FOLDER, recursive = TRUE)
  }
  
  # Search for files matching the pattern
  cat(paste0("Searching for files starting with '", FILE_PREFIX, "'...\n"))
  
  # Build search query
  if (!is.null(GDRIVE_FOLDER) && GDRIVE_FOLDER != "") {
    # Search within specific folder
    folder <- drive_find(pattern = paste0("^", GDRIVE_FOLDER, "$"), 
                         type = "folder")
    
    if (nrow(folder) > 0) {
      folder_id <- folder$id[1]
      files <- drive_ls(path = as_id(folder_id), 
                        pattern = paste0("^", FILE_PREFIX))
    } else {
      cat("Warning: Folder not found. Searching all files.\n")
      files <- drive_find(pattern = paste0("^", FILE_PREFIX))
    }
  } else {
    # Search all files
    files <- drive_find(pattern = paste0("^", FILE_PREFIX))
  }
  
  if (nrow(files) == 0) {
    stop(paste0("No files found matching '", FILE_PREFIX, "'"))
  }
  
  cat(paste0("Found ", nrow(files), " files to download\n\n"))
  
  # Download each file
  downloaded_files <- c()
  
  for (i in 1:nrow(files)) {
    file_name <- files$name[i]
    file_path <- file.path(DOWNLOAD_FOLDER, file_name)
    
    cat(paste0("Downloading ", i, "/", nrow(files), ": ", file_name, "\n"))
    
    drive_download(
      file = as_id(files$id[i]),
      path = file_path,
      overwrite = TRUE
    )
    
    downloaded_files <- c(downloaded_files, file_path)
  }
  
  cat("\nAll files downloaded!\n\n")
  return(downloaded_files)
}
# 
# # ============================================================================
# # OPTION 2: MERGE LOCAL FILES (if already downloaded)
# # ============================================================================
# 
# get_local_tiles <- function(folder_path, file_pattern) {
#   
#   # Find all matching .tif files
#   search_pattern <- paste0(file_pattern, ".*\\.tif$")
#   all_files <- list.files(folder_path, 
#                           pattern = search_pattern, 
#                           full.names = TRUE)
#   
#   if (length(all_files) == 0) {
#     stop(paste0("No files found in '", folder_path, 
#                 "' matching pattern '", file_pattern, "'"))
#   }
#   
#   cat(paste0("Found ", length(all_files), " local tiles\n"))
#   return(all_files)
# }
# 
# # ============================================================================
# # MERGE RASTER TILES
# # ============================================================================

merge_tiles <- function(tile_files, output_file) {
  
  if (length(tile_files) == 0) {
    stop("No tiles to merge!")
  }
  
  cat(paste0("\nMerging ", length(tile_files), " tiles...\n"))
  
  # Read all rasters into a list
  cat("Reading raster tiles...\n")
  raster_list <- lapply(tile_files, rast)
  
  # Create a SpatRasterCollection
  cat("Creating raster collection...\n")
  rast_collection <- sprc(raster_list)
  
  # Merge/mosaic all tiles
  cat("Merging tiles (this may take a while)...\n")
  merged_raster <- mosaic(rast_collection)
  
  # Write the merged raster
  cat(paste0("Writing merged raster to: ", output_file, "\n"))
  writeRaster(merged_raster, 
              output_file, 
              overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))
  
  cat("\nSuccess! Merged raster saved.\n\n")
  
  # Print raster info
  cat("Raster info:\n")
  print(merged_raster)
  
  return(merged_raster)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cat("=" , rep("=", 69), "\n", sep = "")
cat("REEF MASK TILE MERGER\n")
cat("=" , rep("=", 69), "\n\n", sep = "")

# Try downloading from Google Drive first
tryCatch({
  
  cat("OPTION 1: Downloading from Google Drive...\n")
  cat("=" , rep("=", 69), "\n", sep = "")
  
  tile_files <- download_and_merge_from_gdrive()
  merged <- merge_tiles(tile_files, OUTPUT_FILE)
  
}, error = function(e) {
  
  cat("\nError with Google Drive download:", e$message, "\n")
  cat("\nTrying OPTION 2: Merging local files...\n")
  cat("=" , rep("=", 69), "\n", sep = "")
  
  # Try merging local files instead
  tryCatch({
    
    tile_files <- get_local_tiles(DOWNLOAD_FOLDER, FILE_PREFIX)
    merged <- merge_tiles(tile_files, OUTPUT_FILE)
    
  }, error = function(e2) {
    cat("\nError:", e2$message, "\n")
    cat("\nMake sure you have either:\n")
    cat("1. Authenticated with Google Drive (googledrive package), OR\n")
    cat("2. Downloaded the tiles manually to:", DOWNLOAD_FOLDER, "\n")
  })
})

# ============================================================================
# ALTERNATIVE: SIMPLE MANUAL MERGE (if you have files locally)
# ============================================================================

# If the above doesn't work, uncomment this simpler version:
# 
# library(terra)
# 
# # Set your folder containing the tiles
# tile_folder <- "./reef_tiles"
# 
# # Get all .tif files
# tile_files <- list.files(tile_folder, 
#                          pattern = "global_reef_mask_30m_styled.*\\.tif$",
#                          full.names = TRUE)
# 
# # Read and merge
# rasters <- lapply(tile_files, rast)
# collection <- sprc(rasters)
# merged <- mosaic(collection)
# 
# # Save
# writeRaster(merged, "merged_reef_mask.tif", overwrite = TRUE)
# 
# cat("Done! Merged raster saved to: merged_reef_mask.tif\n")










#### intersecting exclusive economic zones with shallow water tropical reefs ####
## done in r 
library(terra)
# Load shapefiles
shape1 <- vect("~/Library/CloudStorage/Dropbox/Coracle/Global_reef_maps/UNEP-WCMC/14_001_WCMC008_CoralReefs2018_v4_1/01_Data/WCMC008_CoralReef2018_Py_v4_1.shp")
shape2 <- vect("~/Library/CloudStorage/Dropbox/Coracle/World_EEZ_v12_20231025/eez_v12.shp")
# Ensure same CRS
shape2 <- project(shape2, crs(shape1))
# Intersect
intersected <- terra::intersect(shape1, shape2)
class(intersected)
writeVector(intersected, "coral_reefs_in_eez.shp", overwrite = TRUE)

## done using QGIS - to make sure everything is correct
intersected_qgis <- vect("Coral_Reefs_in_EEZ_QGIS.shp")

names(intersected_qgis)
head(intersected_qgis)

unique(intersected_qgis$SOVEREIGN1)
unique(intersected_qgis$SOVEREIGN2)
unique(intersected_qgis$SOVEREIGN3)


## plots
library(dplyr)

countries_with_reefs <- intersected_qgis %>%
  as.data.frame() %>%
  select(SOVEREIGN1, SOVEREIGN2, SOVEREIGN3) %>%
  tidyr::pivot_longer(cols = everything(), values_to = "country") %>%
  filter(!is.na(country), country != "") %>%
  distinct(country) %>%
  arrange(country)

reef_counts <- intersected_qgis %>%
  as.data.frame() %>%
  select(SOVEREIGN1, SOVEREIGN2, SOVEREIGN3) %>%
  tidyr::pivot_longer(cols = everything(), values_to = "country") %>%
  filter(!is.na(country), country != "") %>%
  group_by(country) %>%
  summarise(n_reefs = n()) %>%
  arrange(desc(n_reefs))

library(terra)

# Reproject for accurate area calculations
intersected_eqarea <- project(intersected_qgis, "ESRI:54009")

# Add area (km²)
intersected_eqarea$area_km2 <- expanse(intersected_eqarea, unit = "km")

reef_area <- intersected_eqarea %>%
  as.data.frame() %>%
  select(SOVEREIGN1, SOVEREIGN2, SOVEREIGN3, area_km2) %>%
  tidyr::pivot_longer(cols = c(SOVEREIGN1, SOVEREIGN2, SOVEREIGN3), values_to = "country") %>%
  filter(!is.na(country), country != "") %>%
  group_by(country) %>%
  summarise(total_area_km2 = sum(area_km2, na.rm = TRUE)) %>%
  arrange(desc(total_area_km2))

reef_area



library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

# Load world shapefile (country polygons)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Make sure reef_counts$country matches naming in world$name
# You might need to adjust or harmonize names, e.g. "United States" vs "United States of America"
reef_counts <- reef_counts %>%
  mutate(country = case_when(
    country == "United States" ~ "United States of America",
    country == "Democratic Republic of the Congo" ~ "Congo, The Democratic Republic of the",
    TRUE ~ country
  ))

# Merge reef counts with world shapefile
world_reef <- left_join(world, reef_counts, by = c("name" = "country"))

# Plot
ggplot(data = world_reef) +
  geom_sf(aes(fill = n_reefs), color = "grey50", size = 0.1) +
  scale_fill_viridis_c(
    option = "plasma",
    trans = "log10",
    na.value = "lightgrey",
    name = "Reef count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    title = "Distribution of Shallow-Water Tropical Reefs by Country (EEZ Intersections)",
    caption = "Data: UNEP-WCMC Coral Reefs (2018) & Marine Regions EEZs"
  )





