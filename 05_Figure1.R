#### author: Puja Pande
#### description: builds main text Figure 1 (sites map, and model coefficient panels), assembled with patchwork
#### note: large fonts used for screenshot/presentation use; panel labels added manually

library(cowplot)
library(patchwork)
library(ggplot2)
library(scales)
library(extrafont)
library(ncdf4)
library(dplyr)
library(lubridate)
library(fixest)
library(sf)
library(rnaturalearth)
library(terra)
library(marmap)
loadfonts(device = "pdf", quiet = TRUE)

#### Font size constants — Nature guidelines (5-7 pt at final print size) ####
BS  <- 7    # base_size
AT  <- 7    # axis title
ATX <- 6    # axis text
LT  <- 6    # legend text
TAG <- 8    # panel tag (a/b/c), bold
LW  <- 0.5  # thin line
LWM <- 0.8  # medium line

#### ── PANEL DATA ───────────────────────────────────────────────────────────
All_Bleaching_Events_Data_AllDHW <- read.csv("path/to/Coracle/Datasets - 01 July 2026/Final_Panel_Data_allDHW_01July.csv")
All_Bleaching_Events_Data_AllDHW$abs_lat <- abs(All_Bleaching_Events_Data_AllDHW$Latitude_Degrees)
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(Ecoregion_Month = interaction(Ecoregion_Name, Date_Month, drop = TRUE))

#### ── PANEL A: SITES MAP (basemap shared, cached to disk after first build) ────
# Rebuilding the global bathymetry raster (download + reproject) on every run
# is memory-heavy; in a persistent IDE session (unlike a fresh Rscript call)
# that spike stacks on top of whatever else is already loaded and can crash
# the R session. Cache it once and read from disk on subsequent runs.
robin_pac <- "+proj=robin +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Absolute path (not relative to whatever the current session's working
# directory happens to be) so the cache is always found regardless of
# whether this is run via Rscript or sourced inside an IDE session.
cache_dir <- "/Users/pujapande/Library/CloudStorage/Dropbox/Coracle/Positron Test/map_cache"
if (!dir.exists(cache_dir)) dir.create(cache_dir)
world_cache    <- file.path(cache_dir, "world_robin.rds")
bathy_cache    <- file.path(cache_dir, "bathy_robin.tif")
bathy_df_cache <- file.path(cache_dir, "bathy_df_robin.rds")  # cache the DATA FRAME, not just the raster

# Single-threaded GDAL: multi-threaded raster reprojection/conversion has been
# observed to crash the R session under Positron/ark's sandboxed execution
# context — `as.data.frame()` on a SpatRaster segfaults ark specifically, even
# though it's fine under plain Rscript. Set this before ANY terra call.
Sys.setenv(GDAL_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
terra::terraOptions(threads = 1, memfrac = 0.3)

if (file.exists(world_cache)) {
  world <- readRDS(world_cache)
} else {
  world <- ne_countries(scale = "medium", returnclass = "sf") %>%
    st_transform(robin_pac)
  saveRDS(world, world_cache)
}

# ---- bathy_df: this is the crash point. Cache the finished data frame so
# the risky raster->data.frame conversion only ever has to run ONCE, and
# ideally that first run happens via plain `Rscript`, not inside Positron.
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

  # Round-trip through disk instead of a direct in-memory as.data.frame()
  # call — writeRaster + rast() + as.data.frame() on the disk-backed object
  # has been more stable than converting the in-memory SpatRaster directly.
  writeRaster(bathy_robin, file.path(cache_dir, "bathy_tmp_conv.tif"), overwrite = TRUE)
  bathy_df <- as.data.frame(rast(file.path(cache_dir, "bathy_tmp_conv.tif")), xy = TRUE)
  file.remove(file.path(cache_dir, "bathy_tmp_conv.tif"))

  rm(bathy_robin); gc()

  saveRDS(bathy_df, bathy_df_cache)
}

ocean_palette <- colorRampPalette(c("#AED6F1"))(100)

# Latitude graticule (parallels only, no meridians). Parallels are exactly
# straight horizontal lines under the Robinson projection, so the spurious
# chord st_graticule() draws where its longitude sequence crosses the
# antimeridian of this off-centre (lon_0 = 150) projection retraces the same
# horizontal line rather than cutting a visible diagonal across the map —
# safe to use directly without manually splitting at the seam.
lat_graticule <- st_graticule(lat = seq(-60, 60, 20), lon = seq(-180, 180, 20), ndiscr = 200) %>%
  filter(type == "N") %>%
  st_transform(robin_pac)

# Equator + tropics as dashed reference lines, with their latitude labelled
# on the left edge of the map. Parallels are straight horizontal lines under
# this Robinson projection (see note above lat_graticule), so a single
# projected point at any longitude gives the correct y for the label.
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
            hjust = 0, vjust = -0.4, size = ATX / .pt, colour = "grey20",
            family = "Arial", inherit.aes = FALSE)
)

map_theme <- theme_void(base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "#AED6F1", colour = NA),
    panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.4),
    legend.position   = "bottom",
    legend.title      = element_text(size = LT),
    legend.text       = element_text(size = ATX),
    plot.tag          = element_text(size = TAG, face = "bold")
  )

map_coord <- coord_sf(crs = robin_pac, xlim = c(-17000000, 17000000),
                       ylim = c(-4000000, 4000000), expand = FALSE)

# Per-site count of bleaching observations (i.e. surveys), regardless of
# whether bleaching was recorded at each one — a survey reading 0% bleaching
# is still an observation. Built from the raw combined dataset (not the
# panel-regression CSV) so this map's site extent matches Fig. S12's "all
# sites" layer — the panel CSV structurally excludes single-observation
# sites upstream (in Merging_GCBD_MERMAID.R), which was why the two figures
# previously showed different geographic coverage (e.g. Red Sea, northern
# SA / southern Mozambique). Restricted to MERMAID + GCBD-Population rows,
# same filter used for Plot_Supp1_and_12.R's all_sites_sf.
combined_raw <- read.csv("/Users/pujapande/Library/CloudStorage/Dropbox/Coracle/Datasets - 01 July 2026/Final_Combined_Data_all_levels_01July2026.csv") %>%
  filter(Source == "Mermaid" | (Source == "GCBD" & Bleaching_Level == "Population"))

site_bleach_sf <- combined_raw %>%
  group_by(Site_ID) %>%
  mutate(n_bleached = n()) %>%
  distinct(Site_ID, Latitude_Degrees, Longitude_Degrees, n_bleached) %>%
  ungroup() %>%
  st_as_sf(coords = c("Longitude_Degrees", "Latitude_Degrees"), crs = 4326) %>%
  st_transform(robin_pac) %>%
  mutate(bleach_cat = cut(n_bleached,
                          breaks = c(0, 1, 2, 5, 10, Inf),
                          labels = c("1", "2", "3–5", "6–10", "10+"),
                          right  = TRUE))

p_bleach_count <- ggplot() +
  base_map_layers +
  geom_sf(data = site_bleach_sf,
          aes(colour = bleach_cat, size = bleach_cat),
          alpha = 0.8, shape = 16) +
  scale_colour_manual(
    name   = "Number of bleaching observations",
    values = c(
               "1"  = "grey40",
               "2"  = "#F0E442",
               "3–5"  = "#E69F00",
               "6–10" = "#D55E00",
               "10+"  = "#7B3500")
  ) +
  scale_size_manual(
    name   = "Number of bleaching observations",
    values = c("1" = 1, "2" = 1.2, "3–5" = 2, "6–10" = 2.8, "10+" = 3.6),
    guide  = "none"
  ) +
  guides(colour = guide_legend(
    title.position = "top",
    title.hjust    = 0.5,
    label.position = "bottom",
    nrow           = 1,
    override.aes   = list(size = c(1, 1.2, 2, 2.8, 3.6), alpha = 1)
  )) +
  map_coord +
  labs(tag = "a") +
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

#### ── PANEL B ───────────────────────────────────────────────────────────── ####

PLOT_ABS_LAT <- 13
base_path    <- "path/to/Coracle"

modlin      <- read.csv(file.path(base_path, "linear_lat_conley200km_coefs.csv"))
modbin      <- read.csv(file.path(base_path, "FE_binned_DHW_conley200km_coefs.csv"))
modlin_vcov <- readRDS(file.path(base_path, "linear_lat_conley200km_vcov.rds"))
modbin_vcov <- readRDS(file.path(base_path, "FE_binned_DHW_conley200km_vcov.rds"))

stopifnot(
  "modlin_vcov must be 2x2"              = all(dim(modlin_vcov) == c(2, 2)),
  "modlin_vcov row/col names must match" = identical(rownames(modlin_vcov), colnames(modlin_vcov)),
  "modbin_vcov row/col names must match" = identical(rownames(modbin_vcov), colnames(modbin_vcov))
)

b_dhw     <- modlin$x[1]
b_dhw_lat <- modlin$x[2]

parse_bin <- function(term) {
  bin_label <- sub('.*DHW_bin::([^"]+)"\\)', '\\1', term)
  is_lat    <- grepl(":abs_lat", bin_label)
  bin_label <- sub(":abs_lat", "", bin_label)
  list(bin = bin_label, is_lat = is_lat)
}

parsed        <- lapply(modbin$term, parse_bin)
modbin$bin    <- sapply(parsed, `[[`, "bin")
modbin$is_lat <- sapply(parsed, `[[`, "is_lat")

main_terms <- modbin %>% filter(!is_lat) %>% select(bin, estimate) %>% tibble::deframe()
lat_terms  <- modbin %>% filter(is_lat)  %>% select(bin, estimate) %>% tibble::deframe()

get_left_edge <- function(b) as.numeric(sub("\\[([0-9.]+),.*", "\\1", b))

linear_curve <- function(abs_lat, dhw_grid = seq(0, 9, length.out = 200)) {
  x_grid <- cbind(dhw_grid, dhw_grid * abs_lat)
  pred   <- x_grid %*% c(b_dhw, b_dhw_lat)
  se     <- sqrt(rowSums((x_grid %*% modlin_vcov) * x_grid))
  data.frame(dhw = dhw_grid, pred = as.numeric(pred), se = se, abs_lat = abs_lat)
}

binned_estimates <- function(abs_lat) {
  result <- lapply(unique(modbin$bin), function(bin_label) {
    main_name <- paste0('DHW_bin::', bin_label, '")')
    lat_name  <- paste0('DHW_bin::', bin_label, ':abs_lat")')
    pred      <- main_terms[[bin_label]] + lat_terms[[bin_label]] * abs_lat
    var_main  <- modbin_vcov[main_name, main_name]
    var_lat   <- modbin_vcov[lat_name,  lat_name]
    cov_ml    <- modbin_vcov[main_name, lat_name]
    var_pred  <- var_main + abs_lat^2 * var_lat + 2 * abs_lat * cov_ml
    data.frame(dhw = get_left_edge(bin_label), pred = pred,
               se = sqrt(var_pred), abs_lat = abs_lat)
  })
  bind_rows(data.frame(dhw = 0, pred = 0, se = 0, abs_lat = abs_lat),
            bind_rows(result))
}

lat_styles <- c("13°" = "#000000")
lat_labels <- c("13°")

linear_df <- linear_curve(13) %>%
  mutate(lat_label = factor("13°", levels = lat_labels))

bin_df <- binned_estimates(PLOT_ABS_LAT) %>%
  mutate(ci_lo = pred - 1.96 * se, ci_hi = pred + 1.96 * se)

p_binned <- ggplot() +
  geom_ribbon(data = linear_df,
              aes(x = dhw, ymin = pred - 1.96 * se, ymax = pred + 1.96 * se,
                  fill = lat_label), alpha = 0.12) +
  geom_line(data = linear_df,
            aes(x = dhw, y = pred, colour = lat_label), linewidth = LWM) +
  geom_errorbar(data = bin_df,
                aes(x = dhw, ymin = ci_lo, ymax = ci_hi),
                width = 0.15, linewidth = LW, colour = "black") +
  geom_point(data = bin_df, aes(x = dhw, y = pred),
             size = 2, colour = "black", shape = 16) +
  scale_colour_manual(values = lat_styles, labels = lat_labels, name = "Latitude") +
  scale_fill_manual(values   = lat_styles, labels = lat_labels, guide = "none") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 8.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = "Degree heating weeks (°C-weeks)", y = "Coral bleaching (%)", tag = "b") +
  theme_classic(base_size = BS, base_family = "Arial") +
  theme(
    axis.line         = element_blank(),
    panel.border      = element_rect(fill = NA, colour = "black", linewidth = LW),
    axis.ticks        = element_line(linewidth = LW, colour = "black"),
    axis.ticks.length = unit(3, "pt"),
    axis.text         = element_text(size = ATX, colour = "black"),
    axis.title        = element_text(size = AT,  colour = "black"),
    axis.title.y = element_text(size = AT, colour = "black", margin = margin(r = 10)),
    legend.position   = c(0.15, 0.90),
    legend.title      = element_text(size = LT),
    legend.text       = element_text(size = LT),
    legend.key.size   = unit(0.3, "cm"),
    panel.grid        = element_blank(),
    plot.margin       = margin(4, 6, 4, 4)
  )

#### ── PANEL C ───────────────────────────────────────────────────────────── ####

model_linear_lat <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

dhw_max        <- max(All_Bleaching_Events_Data_AllDHW$dhw, na.rm = TRUE)
dhw_seq        <- seq(0, dhw_max, length.out = 300)
median_abs_lat <- median(All_Bleaching_Events_Data_AllDHW$abs_lat, na.rm = TRUE)
y_max          <- 100

# Conley (200 km, spherical) vcov, matching panel B and the reported intervals
vcov_lat <- unclass(vcov_conley(
  model_linear_lat,
  lat      = "Latitude_Degrees",
  lon      = "Longitude_Degrees",
  cutoff   = 200,
  distance = "spherical"
))

beta_dhw     <- coef(model_linear_lat)["dhw"]
beta_dhw_lat <- coef(model_linear_lat)["dhw:abs_lat"]
se_dhw       <- sqrt(vcov_lat["dhw", "dhw"])
se_dhw_lat   <- sqrt(vcov_lat["dhw:abs_lat", "dhw:abs_lat"])
cov_dhw      <- vcov_lat["dhw", "dhw:abs_lat"]

make_linear_curve <- function(lat, label) {
  data.frame(dhw = dhw_seq, lat_group = label) %>%
    mutate(
      pred     = (beta_dhw + beta_dhw_lat * lat) * dhw,
      var_pred = dhw^2 * (se_dhw^2 + lat^2 * se_dhw_lat^2 + 2 * lat * cov_dhw),
      se_pred  = sqrt(pmax(var_pred, 0)),
      ci_lo    = pmax(pred - 1.96 * se_pred, 0),
      ci_hi    = pmin(pred + 1.96 * se_pred, y_max)
    )
}

linear_curves <- bind_rows(
  make_linear_curve(0,              "0°"),
  make_linear_curve(13, "13°"),
  make_linear_curve(23,             "23°")
) %>%
  mutate(lat_group = factor(lat_group, levels = c("0°", "13°", "23°")))

curve_colours <- c("13°" = "#000000", "0°" = "#FF4086", "23°" = "#35A1FF")
curve_lty     <- c("13°" = "solid",   "0°" = "dashed",  "23°" = "dashed")

p_bl_nature <- ggplot() +
  geom_ribbon(data = linear_curves,
              aes(x = dhw, ymin = ci_lo, ymax = ci_hi, fill = lat_group), alpha = 0.12) +
  geom_line(data = linear_curves,
            aes(x = dhw, y = pred, colour = lat_group, linetype = lat_group),
            linewidth = LWM) +
  scale_colour_manual(values   = curve_colours) +
  scale_fill_manual(values     = curve_colours, guide = "none") +
  scale_linetype_manual(values = curve_lty) +
  guides(
    colour   = guide_legend(title = "Latitude",
                            override.aes = list(linewidth = 0.8, fill = NA)),
    linetype = guide_legend(title = "Latitude",
                            override.aes = list(linewidth = 0.8, fill = NA))
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 8.5)) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = "Degree heating weeks (°C-weeks)", y = "Coral bleaching (%)", tag = 'c') +
  theme_classic(base_size = BS, base_family = "Arial") +
  theme(
    axis.line         = element_blank(),
    panel.border      = element_rect(fill = NA, colour = "black", linewidth = LW),
    axis.ticks        = element_line(linewidth = LW, colour = "black"),
    axis.ticks.length = unit(3, "pt"),
    axis.text         = element_text(size = ATX, colour = "black"),
    axis.title        = element_text(size = AT,  colour = "black"),
    legend.text       = element_text(size = LT,  colour = "black"),
    axis.title.y = element_text(size = AT, colour = "black", margin = margin(r = 10)),
    legend.title      = element_text(size = LT,  colour = "black"),
    legend.position   = c(0.15, 0.88),
    legend.key.width  = unit(0.5, "cm"),
    legend.key.height = unit(0.2, "cm"),
    legend.spacing.y  = unit(0.05, "cm"),
    panel.grid        = element_blank(),
    plot.margin       = margin(4, 6, 4, 4)
  )

#### ── ASSEMBLE & SAVE ───────────────────────────────────────────────────── ####
# p_binned/p_bl_nature previously had aspect.ratio = 1, which forced them to
# render as squares narrower than their allotted column — that's what made
# the bottom row occupy less width than panel a's map above it, no matter
# how the spacer widths were tuned. Dropping the fixed aspect ratio lets both
# panels stretch to fill their full column, so the bottom row now spans the
# same total width as panel a without needing manual spacer compensation.
fig <- p_bleach_count / (p_binned | p_bl_nature) +
  plot_layout(heights = c(1, 1))

fig

ggsave("figure1_map_panel_a.png", fig,
       width  = 183,
       height = 150,
       units  = "mm",
       dpi    = 600)



