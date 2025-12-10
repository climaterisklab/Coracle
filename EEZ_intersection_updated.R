library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(scales)

# ============================================================================
# Load Data
# ============================================================================

# Path to your NEW exported QGIS results
reef_eez <- st_read("/Users/pujapande/Documents/Coracle/eez_finally_working_09_12.gpkg")

# Check structure and find the count field name
glimpse(reef_eez)
names(reef_eez)

# IMPORTANT: Check what your count field is called!
# It might be "fid_count", "reef_count", "DN_count", etc.
# Replace "fid_count" below with the actual name

# ============================================================================
# Calculate Reef Areas (if counting discrete polygons)
# ============================================================================

# If you have a count field (number of reef polygons per EEZ):
reef_eez <- reef_eez %>%
  mutate(
    # Rename your count field to a standard name
    reef_polygon_count = fid_count  # CHANGE 'fid_count' to your actual field name!
  )

# ============================================================================
# Map: Reef Polygon Count by EEZ
# ============================================================================

ggplot(reef_eez) +
  geom_sf(aes(fill = reef_polygon_count), color = "grey30", size = 0.1) +
  scale_fill_viridis_c(
    option = "turbo",
    name = "Number of\nReef Areas",
    trans = "log10",
    na.value = "lightgrey",
    breaks = c(1, 10, 100, 1000, 10000),
    labels = comma
  ) +
  labs(
    title = "Number of Discrete Coral Reef Areas by EEZ",
    subtitle = "Derived from Allen Coral Atlas 30m resolution",
    caption = "Data: UNEP-WCMC EEZ v12 + Allen Coral Atlas"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ============================================================================
# Bar Chart: Top 15 EEZs by Reef Count
# ============================================================================

reef_eez %>%
  st_drop_geometry() %>%
  filter(reef_polygon_count > 0) %>%
  arrange(desc(reef_polygon_count)) %>%
  slice_head(n = 15) %>%
  ggplot(aes(x = reorder(TERRITORY1, reef_polygon_count), y = reef_polygon_count)) +
  geom_col(fill = "coral") +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    x = NULL,
    y = "Number of Reef Areas",
    title = "Top 15 EEZs by Number of Coral Reef Areas"
  ) +
  theme_minimal()

# ============================================================================
# Aggregate by Sovereign State
# ============================================================================

eez_by_sovereign <- reef_eez %>%
  st_drop_geometry() %>%
  filter(reef_polygon_count > 0) %>%
  select(SOVEREIGN1, SOVEREIGN2, SOVEREIGN3, reef_polygon_count) %>%
  pivot_longer(
    cols = c(SOVEREIGN1, SOVEREIGN2, SOVEREIGN3), 
    names_to = "sovereign_type", 
    values_to = "sovereign"
  ) %>%
  filter(!is.na(sovereign), sovereign != "") %>%
  group_by(sovereign) %>%
  summarise(
    total_reef_areas = sum(reef_polygon_count, na.rm = TRUE),
    n_eezs = n()
  ) %>%
  arrange(desc(total_reef_areas))

# View top countries
cat("\n=== Top 20 Sovereign States by Number of Reef Areas ===\n")
print(head(eez_by_sovereign, 20))

# ============================================================================
# Top 20 Sovereigns Bar Chart
# ============================================================================

top20 <- head(eez_by_sovereign, 20)

bar_plot <- ggplot(top20, aes(x = reorder(sovereign, total_reef_areas), 
                              y = total_reef_areas)) +
  geom_col(fill = "coral") +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Top 20 Sovereign States by Number of Coral Reef Areas",
    x = "Sovereign State",
    y = "Total Number of Reef Areas"
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

# Harmonize country names
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
  geom_sf(aes(fill = total_reef_areas), color = "grey50", size = 0.1) +
  scale_fill_viridis_c(
    option = "plasma",
    trans = "log10",
    na.value = "lightgrey",
    name = "Number of\nReef Areas",
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
    title = "Distribution of Coral Reef Areas by Country (EEZ)",
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
  "coral_reef_areas_by_sovereign.png",
  combined,
  width = 14,
  height = 12,
  dpi = 300,
  bg = "white"
)

cat("\nPlot saved as: coral_reef_areas_by_sovereign.png\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=== SUMMARY STATISTICS ===\n")
cat(sprintf("Total sovereign states with coral: %d\n", nrow(eez_by_sovereign)))
cat(sprintf("Total discrete reef areas globally: %s\n", 
            comma(sum(eez_by_sovereign$total_reef_areas))))
cat(sprintf("\nTop 5 states contain %.1f%% of all reef areas\n",
            100 * sum(head(eez_by_sovereign, 5)$total_reef_areas) / 
              sum(eez_by_sovereign$total_reef_areas)))