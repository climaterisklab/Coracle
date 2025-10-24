### author: Puja Pande
### date: 21 October 2025
### description: EEZ intersections with shallow water tropical reefs












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





