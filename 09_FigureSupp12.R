#### author: Puja Pande
#### description: builds Supplementary Figure S12 — global map of panel bleaching
####              survey sites, coloured/sized by number of surveys per site, over
####              a bathymetry basemap

######## panel bleaching sites 
All_Bleaching_Events_Data_AllDHW <- read.csv("path/to/Coracle/Datasets/Final_Panel_Data_allDHW_01July.csv")

library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
library(terra)
library(marmap)

robin_pac <- "+proj=robin +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Shares figure_1_updated.R's on-disk cache so the bathymetry raster and
# world polygons are only ever built once across scripts, rather than
# re-downloading/reprojecting (memory-heavy, and has crashed Positron's ark
# kernel before).
cache_dir <- "path/to/Coracle/Positron Test/map_cache"
if (!dir.exists(cache_dir)) dir.create(cache_dir)
world_cache    <- file.path(cache_dir, "world_robin.rds")
bathy_cache    <- file.path(cache_dir, "bathy_robin.tif")
bathy_df_cache <- file.path(cache_dir, "bathy_df_robin.rds")

Sys.setenv(GDAL_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
terra::terraOptions(threads = 1, memfrac = 0.3)

if (file.exists(world_cache)) {
  world <- readRDS(world_cache)
} else {
  world <- ne_countries(scale = "medium", returnclass = "sf") %>%
    st_transform(robin_pac)
  saveRDS(world, world_cache)
}

if (file.exists(bathy_df_cache)) {
  bathy_df <- readRDS(bathy_df_cache)
} else {
  message(
    "No cached bathy_df found. This step (raster -> data.frame conversion) ",
    "has crashed Positron's ark kernel before. Strongly recommended: run this ",
    "block once via `Rscript this_file.R` in a terminal instead of inside ",
    "Positron, then just re-source this script in Positron afterwards — it ",
    "will pick up the cache and skip straight past this block."
  )

  if (file.exists(bathy_cache)) {
    bathy_robin <- rast(bathy_cache)
  } else {
    bathy_raw <- getNOAA.bathy(
      lon1 = -180, lon2 = 180,
      lat1 = -90,  lat2 = 90,
      resolution = 20, keep = TRUE
    )

    bathy_xyz <- as.xyz(bathy_raw) %>%
      as.data.frame() %>%
      setNames(c("lon", "lat", "depth"))
    rm(bathy_raw); gc()

    bathy_terra <- rast(bathy_xyz, type = "xyz", crs = "EPSG:4326")
    rm(bathy_xyz); gc()
    bathy_terra[bathy_terra > 0] <- NA
    bathy_robin <- project(bathy_terra, robin_pac, method = "bilinear")
    rm(bathy_terra); gc()

    writeRaster(bathy_robin, bathy_cache, overwrite = TRUE)
  }

  writeRaster(bathy_robin, file.path(cache_dir, "bathy_tmp_conv.tif"), overwrite = TRUE)
  bathy_df <- as.data.frame(rast(file.path(cache_dir, "bathy_tmp_conv.tif")), xy = TRUE)
  file.remove(file.path(cache_dir, "bathy_tmp_conv.tif"))

  rm(bathy_robin); gc()

  saveRDS(bathy_df, bathy_df_cache)
}

ocean_palette <- colorRampPalette(c("#AED6F1"))(100)

lat_graticule <- st_graticule(lat = seq(-60, 60, 20), lon = seq(-180, 180, 20), ndiscr = 200) %>%
  filter(type == "N") %>%
  st_transform(robin_pac)

ref_lats   <- c(-23.4368, 0, 23.4368)
ref_labels <- c("23.5°S", "0°", "23.5°N")

ref_graticule <- st_graticule(lat = ref_lats, lon = seq(-180, 180, 20), ndiscr = 200) %>%
  filter(type == "N") %>%
  st_transform(robin_pac)

ref_label_df <- data.frame(lat = ref_lats, label = ref_labels)
ref_label_df$y <- vapply(ref_lats, function(l) {
  pt <- st_sfc(st_point(c(0, l)), crs = 4326) %>% st_transform(robin_pac)
  st_coordinates(pt)[2]
}, numeric(1))
ref_label_df$x <- -16700000

base_map_layers <- list(
  geom_raster(data = bathy_df, aes(x = x, y = y, fill = depth)),
  scale_fill_gradientn(colours = ocean_palette, limits = c(-6000, 0),
                       na.value = "#AED6F1", guide = "none"),
  geom_sf(data = world, fill = "grey85", colour = "grey60", linewidth = 0.2),
  #geom_sf(data = lat_graticule, colour = "grey50", linewidth = 0.15, alpha = 0.4),
  geom_sf(data = ref_graticule, colour = "grey55", linewidth = 0.35,
          linetype = "dotted", alpha = 0.8),
  geom_text(data = ref_label_df, aes(x = x, y = y, label = label),
            hjust = 0, vjust = -0.4, size = 6 / .pt, colour = "grey20",
            family = "Arial", inherit.aes = FALSE)
)

map_theme <- theme_void(base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "#AED6F1", colour = NA),
    panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.4),
    legend.position   = "bottom",
    legend.title      = element_text(size = 6),
    legend.text       = element_text(size = 6)
  )

map_coord <- coord_sf(crs = robin_pac, xlim = c(-17000000, 17000000),
                       ylim = c(-4000000, 4000000), expand = FALSE)

#### ── SITE SETS ──────────────────────────────────────────────────────────
combined_raw <- read.csv("/path/to/Coracle/Datasets - 01 July 2026/Final_Combined_Data_all_levels_01July2026.csv")

# Keep all MERMAID rows (Bleaching_Level is always NA there — MERMAID doesn't
# distinguish population/colony level) plus only GCBD rows explicitly marked
# "Population" — this drops both GCBD "Colony"-level rows and the ~6.7k GCBD
# rows with unspecified Bleaching_Level.
combined_raw <- combined_raw %>%
  filter(Source == "Mermaid" | (Source == "GCBD" & Bleaching_Level == "Population"))

all_sites_sf <- combined_raw %>%
  distinct(Site_ID, Latitude_Degrees, Longitude_Degrees) %>%
  st_as_sf(coords = c("Longitude_Degrees", "Latitude_Degrees"), crs = 4326) %>%
  st_transform(robin_pac)

panel_site_ids <- All_Bleaching_Events_Data_AllDHW %>%
  group_by(Site_ID) %>%
  filter(n() > 1) %>%
  pull(Site_ID) %>%
  unique()

panel_sites_sf <- All_Bleaching_Events_Data_AllDHW %>%
  filter(Site_ID %in% panel_site_ids) %>%
  distinct(Site_ID, Latitude_Degrees, Longitude_Degrees) %>%
  st_as_sf(coords = c("Longitude_Degrees", "Latitude_Degrees"), crs = 4326) %>%
  st_transform(robin_pac)

p_sites <- ggplot() +
  base_map_layers +
  geom_sf(data = all_sites_sf, aes(colour = "All sites (GCBD + MERMAID)"),
          size = 1, alpha = 0.6, shape = 16) +
  geom_sf(data = panel_sites_sf, aes(colour = "Panel sites"),
          size = 1.4, alpha = 0.85, shape = 16) +
  scale_colour_manual(
    name   = NULL,
    values = c("All sites (GCBD + MERMAID)" = "black",
               "Panel sites"                 = "#CC0000")
  ) +
  guides(colour = guide_legend(
    title.position = "top",
    title.hjust    = 0.5,
    label.position = "right",
    nrow           = 1,
    override.aes   = list(size = c(2.5, 2.5), alpha = 1)
  )) +
  map_coord +
  map_theme +
  theme(
    legend.direction    = "horizontal",
    legend.position     = "bottom",
    legend.key          = element_blank(),
    legend.key.width    = unit(0.5, "cm"),
    legend.spacing.x    = unit(0.1, "cm"),
    legend.box.spacing  = unit(4, "pt"),
    legend.margin       = margin(t = 0, b = 4),
    plot.margin         = margin(t = 0, r = 2, b = 0, l = 2)
  )

p_sites

ggsave("supp_panel_sites_map.png", p_sites,
       width  = 183,
       height = 100,
       units  = "mm",
       dpi    = 600)

