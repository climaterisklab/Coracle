### author: Puja Pande
### date: 21 October 2025
### description: EEZ intersections with shallow water tropical reefs

#### libraries ####
library(googledrive)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)


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


#### did intersection in QGIS ####
# Read the CSV
eez_coral <- read.csv("/Users/pujapande/Documents/Coracle/zonal_statisticas_coral_eez.csv")

colnames(eez_coral)
head(eez_coral)
summary(eez_coral)

# Filter only EEZs with coral and pivot to get all sovereigns
eez_by_sovereign <- eez_coral %>%
  filter(coral_count > 0) %>%
  select(SOVEREIGN1, SOVEREIGN2, coral_count) %>%
  pivot_longer(cols = c(SOVEREIGN1, SOVEREIGN2), 
               names_to = "sovereign_type", 
               values_to = "sovereign") %>%
  filter(!is.na(sovereign), sovereign != "") %>%
  group_by(sovereign) %>%
  summarise(total_coral_pixels = sum(coral_count, na.rm = TRUE)) %>%
  arrange(desc(total_coral_pixels))

# View top countries
print("Top 20 sovereign states:")
print(head(eez_by_sovereign, 20))

# Top 20 bar chart
top20 <- head(eez_by_sovereign, 20)

ggplot(top20, aes(x = reorder(sovereign, total_coral_pixels), y = total_coral_pixels)) +
  geom_col(fill = "coral") +
  coord_flip() +
  labs(title = "Top 20 Sovereign States by Coral Reef Coverage",
       x = "Sovereign State",
       y = "Total Coral Pixels") +
  theme_minimal() +
  theme(text = element_text(size = 12))


# Aggregate coral counts by sovereign state
reef_counts <- eez_coral %>%
  filter(coral_count > 0) %>%
  select(SOVEREIGN1, SOVEREIGN2, coral_count) %>%
  pivot_longer(cols = c(SOVEREIGN1, SOVEREIGN2), 
               names_to = "sovereign_type", 
               values_to = "country") %>%
  filter(!is.na(country), country != "") %>%
  group_by(country) %>%
  summarise(total_coral_pixels = sum(coral_count, na.rm = TRUE),
            n_reefs = n()) %>%
  arrange(desc(total_coral_pixels))

# Load world shapefile
world <- ne_countries(scale = "medium", returnclass = "sf")

# Harmonize country names
reef_counts_clean <- reef_counts %>%
  mutate(country = case_when(
    country == "United States" ~ "United States of America",
    country == "United Kingdom" ~ "United Kingdom",
    country == "France" ~ "France",
    country == "Netherlands" ~ "Netherlands",
    country == "Australia" ~ "Australia",
    country == "New Zealand" ~ "New Zealand",
    country == "Norway" ~ "Norway",
    TRUE ~ country
  ))

# Merge reef counts with world shapefile
world_reef <- left_join(world, reef_counts_clean, by = c("name" = "country"))

# Create the map with total coral pixels (NO scientific notation)
ggplot(data = world_reef) +
  geom_sf(aes(fill = total_coral_pixels), color = "grey50", size = 0.1) +
  scale_fill_viridis_c(
    option = "plasma",
    trans = "log10",
    na.value = "lightgrey",
    name = "Coral Pixels",
    breaks = c(10, 100, 1000, 10000, 100000),
    labels = comma  # This formats numbers with commas
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    title = "Distribution of Coral Reefs by Country (EEZ)",
    caption = "Data: Coral reef raster clipped to EEZ boundaries"
  )
