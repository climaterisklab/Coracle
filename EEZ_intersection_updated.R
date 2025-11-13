library(sf)
library(ggplot2)
library(dplyr)

# Path to your exported QGIS results
reef_eez <- st_read("/Users/pujapande/Documents/Coracle/zonal_stats_final.gpkg")

# Check structure
glimpse(reef_eez)


reef_eez <- reef_eez %>%
  mutate(
    reef_area_km2 = reef_sum * (30 * 30) / 1e6,   # assuming 30 m pixels
    reef_share = reef_sum / reef_count
  )

library(ggplot2)

ggplot(reef_eez) +
  geom_sf(aes(fill = AREA_KM2), color = "grey30", size = 0.1) +
  scale_fill_viridis_c(
    option = "turbo",
    name = "Reef area (km²)",
    trans = "sqrt",
    na.value = "lightgrey"
  ) +
  labs(
    title = "Estimated Coral Reef Area by EEZ",
    subtitle = "Derived from 30 m global reef mask",
    caption = "Data: UNEP-WCMC EEZ v12 + reef mask zonal statistics"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )



reef_eez %>%
  st_drop_geometry() %>%
  arrange(desc(reef_area_km2)) %>%
  slice_head(n = 15) %>%
  ggplot(aes(x = reorder(TERRITORY1, reef_area_km2), y = reef_area_km2)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Reef area (km²)",
    title = "Top 15 EEZs by coral reef area"
  ) +
  theme_minimal()



#### better plots ####
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(scales)

# Assuming reef_eez is already loaded
# reef_eez <- st_read("your_file.shp")

# ============================================================================
# Aggregate by Sovereign State
# ============================================================================

# Filter only EEZs with coral and pivot to get all sovereigns
eez_by_sovereign <- reef_eez %>%
  st_drop_geometry() %>%
  filter(reef_count > 0) %>%
  select(SOVEREIGN1, SOVEREIGN2, SOVEREIGN3, reef_count) %>%
  pivot_longer(cols = c(SOVEREIGN1, SOVEREIGN2, SOVEREIGN3), 
               names_to = "sovereign_type", 
               values_to = "sovereign") %>%
  filter(!is.na(sovereign), sovereign != "") %>%
  group_by(sovereign) %>%
  summarise(
    total_coral_pixels = sum(reef_count, na.rm = TRUE),
    total_coral_area_km2 = sum(reef_count * 900 / 1e6, na.rm = TRUE),  # Convert to km²
    n_eezs = n()
  ) %>%
  arrange(desc(total_coral_area_km2))

# View top countries
cat("\n=== Top 20 Sovereign States by Coral Coverage ===\n")
print(head(eez_by_sovereign, 20))

# ============================================================================
# Top 20 Bar Chart
# ============================================================================

top20 <- head(eez_by_sovereign, 20)

bar_plot <- ggplot(top20, aes(x = reorder(sovereign, total_coral_area_km2), 
                              y = total_coral_area_km2)) +
  geom_col(fill = "coral") +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Top 20 Sovereign States by Coral Reef Coverage",
    x = "Sovereign State",
    y = "Total Coral Area (km²)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 10)
  )

print(bar_plot)

# ============================================================================
# World Map by Sovereign
# ============================================================================

# Load world shapefile
world <- ne_countries(scale = "medium", returnclass = "sf")

# Harmonize country names between datasets
reef_counts_clean <- eez_by_sovereign %>%
  mutate(sovereign = case_when(
    sovereign == "United States" ~ "United States of America",
    sovereign == "United Kingdom" ~ "United Kingdom",
    sovereign == "France" ~ "France",
    sovereign == "Netherlands" ~ "Netherlands",
    sovereign == "Australia" ~ "Australia",
    sovereign == "New Zealand" ~ "New Zealand",
    sovereign == "Norway" ~ "Norway",
    sovereign == "Tanzania" ~ "United Republic of Tanzania",
    sovereign == "Venezuela" ~ "Venezuela",
    sovereign == "Vietnam" ~ "Vietnam",
    sovereign == "Dominican Republic" ~ "Dominican Republic",
    sovereign == "Micronesia" ~ "Federated States of Micronesia",
    TRUE ~ sovereign
  ))

# Merge reef counts with world shapefile
world_reef <- left_join(world, reef_counts_clean, by = c("name" = "sovereign"))

# Create the map
map_plot <- ggplot(data = world_reef) +
  geom_sf(aes(fill = total_coral_area_km2), color = "grey50", size = 0.1) +
  scale_fill_viridis_c(
    option = "plasma",
    trans = "log10",
    na.value = "lightgrey",
    name = "Coral Area (km²)",
    breaks = c(1, 10, 100, 1000, 10000),
    labels = comma
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.caption = element_text(size = 9, color = "grey50")
  ) +
  labs(
    title = "Distribution of Coral Reefs by Country (EEZ)",
    caption = "Data: Allen Coral Atlas 30m resolution clipped to EEZ boundaries"
  )

print(map_plot)

# ============================================================================
# Combined Plot
# ============================================================================

library(patchwork)

combined <- map_plot / bar_plot + 
  plot_layout(heights = c(2, 1.5))

print(combined)

# Save high-resolution version
ggsave(
  "coral_reef_by_sovereign.png",
  combined,
  width = 14,
  height = 12,
  dpi = 300,
  bg = "white"
)

cat("\nPlot saved as: coral_reef_by_sovereign.png\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=== SUMMARY STATISTICS ===\n")
cat(sprintf("Total sovereign states with coral: %d\n", nrow(eez_by_sovereign)))
cat(sprintf("Total coral area globally: %s km²\n", 
            comma(sum(eez_by_sovereign$total_coral_area_km2))))
cat(sprintf("\nTop 5 states contain %.1f%% of global coral area\n",
            100 * sum(head(eez_by_sovereign, 5)$total_coral_area_km2) / 
              sum(eez_by_sovereign$total_coral_area_km2)))


