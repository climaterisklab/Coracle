#### author: Puja Pande
#### description: builds Supplementary Figure S12 — global map of panel bleaching
####              survey sites, coloured/sized by number of surveys per site, over
####              a bathymetry basemap

######## panel bleaching sites 
All_Bleaching_Events_Data_AllDHW <- read.csv("path/to/Coracle/Datasets/Final_Panel_Data_allDHW_01July.csv")

library(ggplot2)
library(sf)
library(rnaturalearth)
library(tidyterra)

robin_pac <- "+proj=robin +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(robin_pac)

panel_sf <- All_Bleaching_Events_Data_AllDHW %>%
  distinct(Site_ID, Latitude_Degrees, Longitude_Degrees) %>%
  st_as_sf(coords = c("Longitude_Degrees", "Latitude_Degrees"), crs = 4326) %>%
  st_transform(robin_pac)

library(marmap)

bathy_raw <- getNOAA.bathy(
  lon1 = -180, lon2 = 180,
  lat1 = -90,  lat2 = 90,
  resolution = 20, keep = TRUE
)

bathy_xyz   <- as.xyz(bathy_raw) %>%
  as.data.frame() %>%
  setNames(c("lon", "lat", "depth"))

bathy_terra <- rast(bathy_xyz, type = "xyz", crs = "EPSG:4326")
bathy_terra[bathy_terra > 0] <- NA
bathy_robin <- project(bathy_terra, robin_pac, method = "bilinear")

ocean_palette <- colorRampPalette(c("#4a9fd4", "#AED6F1"))(100)

# Add n_obs to panel_sf
panel_sf <- All_Bleaching_Events_Data_AllDHW %>%
  group_by(Site_ID) %>%
  mutate(n_obs = n()) %>%
  distinct(Site_ID, Latitude_Degrees, Longitude_Degrees, n_obs) %>%
  st_as_sf(coords = c("Longitude_Degrees", "Latitude_Degrees"), crs = 4326) %>%
  st_transform(robin_pac)

panel_sf <- panel_sf %>%
  mutate(obs_cat = cut(n_obs,
                       breaks = c(0, 10, 20, 30, Inf),
                       labels = c("2–10", "11–20", "21–30", "30+"),
                       right  = TRUE))

p_sites <- ggplot() +
  geom_spatraster(data = bathy_robin) +
  scale_fill_gradientn(colours = ocean_palette, limits = c(-6000, 0),
                       na.value = "#AED6F1", guide = "none") +
  geom_sf(data = world, fill = "grey85", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = panel_sf,
          aes(colour = obs_cat, size = obs_cat),
          alpha = 0.8, shape = 16) +
  scale_colour_manual(
    name   = "No. of surveys",
    values = c("2–10"  = "#F0E442",  # yellow
               "11–20" = "#E69F00",  # orange
               "21–30" = "#D55E00",  # vermillion
               "30+"   = "#7B3500")  # dark brown
  ) +
  scale_size_manual(
    name   = "No. of surveys",
    values = c("2–10" = 3, "11–20" = 5, "21–30" = 7, "30+" = 9)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  coord_sf(crs = robin_pac, xlim = c(-17000000, 17000000),
           ylim = c(-4000000, 4000000), expand = FALSE) +
  theme_void(base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "#AED6F1", colour = NA),
    panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.4),
    legend.position  = "bottom",
    legend.title     = element_text(size = 17),
    legend.text      = element_text(size = 15)
  )

p_sites
