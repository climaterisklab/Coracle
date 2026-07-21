#### author: Puja Pande
#### description: builds main text Figure 1 (emissions/GMT time series, DHW response
####              curves, and model coefficient panels), assembled with patchwork
#### note: large fonts used for screenshot/presentation use; panel labels added manually

# FIGURE 1 — Large fonts for screenshot/presentation use
# Panel labels will be added manually

library(cowplot)
library(patchwork)
library(ggplot2)
library(scales)
library(extrafont)
library(ncdf4)
library(dplyr)
library(lubridate)
library(fixest)
loadfonts(device = "pdf", quiet = TRUE)

#### Font size constants — large for screen ####
BS  <- 14   # base_size
AT  <- 13   # axis title
ATX <- 11   # axis text
LT  <- 11   # legend text
LW  <- 0.5  # thin line
LWM <- 0.8  # medium line

#### ── DATA ─────────────────────────────────────────────────────────────── ####

emissions <- read.csv("path/to/Coracle/Scripts/cumulative_global_emissions.csv",
                      row.names = 1)
emissions$Year  <- seq(1850, by = 1, length.out = nrow(emissions))
emissions$World <- emissions$World * 3.664 / 1000

nc_gmt   <- nc_open("path/to/Coracle/Scripts/ERA5_GMT.nc")
gmt_vals <- ncvar_get(nc_gmt, "GMT")
gmt_year <- nc_gmt$dim$year$vals
nc_close(nc_gmt)
gmt_df   <- data.frame(year = gmt_year, GMT = gmt_vals)

data <- read.csv("path/to/Coracle/Datasets/Final_Combined_Data_all_levels_01July2026.csv")
data <- data %>%
  filter(!(Source == "GCBD" & is.na(Bleaching_Level))) %>%
  filter(Bleaching_Level == "Population" | is.na(Bleaching_Level))

data$abs_lat <- abs(data$Latitude_Degrees)

left_df <- left_join(emissions, gmt_df, by = c("Year" = "year")) %>%
  filter(Year >= 1985) %>%
  mutate(World = World / 1000)

nc_dhw     <- nc_open("path/to/data/DHWmm_1985-2025.nc")
time_vals  <- ncvar_get(nc_dhw, "time")
time_units <- ncatt_get(nc_dhw, "time", "units")$value
origin     <- as.Date(sub("days since ", "", time_units))
dates      <- origin + time_vals
years      <- year(dates)
dhw_vals   <- ncvar_get(nc_dhw, "DHW")
dhw_vals[is.nan(dhw_vals)] <- NA
nc_close(nc_dhw)

dhw_annual <- data.frame(
  year = unique(years),
  mean_dhw  = sapply(unique(years), function(y) {
    s <- apply(dhw_vals[, years == y], 1, max, na.rm = TRUE)
    mean(s, na.rm = TRUE) }),
  ci_lo_dhw = sapply(unique(years), function(y) {
    s <- apply(dhw_vals[, years == y], 1, max, na.rm = TRUE); n <- sum(!is.na(s))
    mean(s, na.rm = TRUE) - 1.96 * sd(s, na.rm = TRUE) / sqrt(n) }),
  ci_hi_dhw = sapply(unique(years), function(y) {
    s <- apply(dhw_vals[, years == y], 1, max, na.rm = TRUE); n <- sum(!is.na(s))
    mean(s, na.rm = TRUE) + 1.96 * sd(s, na.rm = TRUE) / sqrt(n) })
)

bleaching_obs <- data %>%
  group_by(Date_Year, Site_ID) %>%
  summarise(site_mean = mean(Percent_Bleached, na.rm = TRUE), .groups = "drop") %>%
  group_by(Date_Year) %>%
  summarise(
    mean_bleaching_obs = mean(site_mean, na.rm = TRUE),
    ci_lo_bleaching    = pmax(mean(site_mean) - 1.96 * sd(site_mean) / sqrt(n()), 0),
    ci_hi_bleaching    = mean(site_mean) + 1.96 * sd(site_mean) / sqrt(n()),
    .groups = "drop"
  ) %>% rename(year = Date_Year)

left_df <- left_df %>%
  left_join(dhw_annual,    by = c("Year" = "year")) %>%
  left_join(bleaching_obs, by = c("Year" = "year"))

left_df_1985 <- left_df %>%
  rename(year = Year) %>%
  mutate(
    ci_lo_bleaching = pmax(ci_lo_bleaching, 0),
    ci_hi_bleaching = pmax(ci_hi_bleaching, 0)
  )
left_df_1985$ci_lo_bleaching[37] <- 4.2
left_df_1985$ci_hi_bleaching[37] <- 4.2

gmt_range <- range(left_df_1985$GMT, na.rm = TRUE)
GMT_LO    <- floor(gmt_range[1]   * 2) / 2
GMT_HI    <- ceiling(gmt_range[2] * 2) / 2


#### ── PANEL A DONORS ────────────────────────────────────────────────────── ####

p_spacer <- ggplot() + theme_void()

p_emissions_donor <- ggplot(left_df_1985, aes(x = year, y = World)) +
  geom_line(colour = "#009E73", linewidth = LW) +
  scale_x_continuous(limits = c(1985, 2024), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 0.4), expand = c(0, 0)) +
  labs(y = expression("Cumulative emissions (TtCO"[2]*")"), tag = "a") +
  theme_classic(base_size = BS, base_family = "Arial") +
  theme(
    axis.title.y = element_text(colour = "#009E73", size = AT, vjust = 0),
    axis.text.y  = element_text(colour = "#009E73", size = ATX),
    axis.ticks.y = element_line(colour = "#009E73", linewidth = LW),
    axis.line.y  = element_line(colour = "#009E73", linewidth = LW),
    plot.margin  = margin(0, 0, 0, 0)
  )

p_gmst_donor <- ggplot(left_df_1985, aes(x = year, y = GMT)) +
  geom_line(alpha = 0.5, linetype = "dashed", colour = "grey40", linewidth = LW) +
  scale_x_continuous(limits = c(1985, 2024), expand = c(0, 0)) +
  scale_y_continuous(limits = c(GMT_LO, GMT_HI),
                     breaks = seq(GMT_LO, GMT_HI, by = 0.4), expand = c(0, 0)) +
  labs(y = "GMST (°C)") +
  theme_classic(base_size = BS, base_family = "Arial") +
  theme(
    axis.title.y = element_text(colour = "grey40", size = AT, vjust = 0),
    axis.text.y  = element_text(colour = "grey40", size = ATX),
    axis.ticks.y = element_line(colour = "grey40", linewidth = LW),
    axis.line.y  = element_line(colour = "grey40", linewidth = LW),
    plot.margin  = margin(0, 0, 0, -1000)
  )

p_dhw_donor <- ggplot(left_df_1985, aes(x = year, y = mean_dhw)) +
  geom_line(colour = "#0072B2", linewidth = LW) +
  scale_x_continuous(limits = c(1985, 2024), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2),
                     position = "right", expand = c(0, 0)) +
  labs(y = "Mean DHW across sites (°C-weeks)") +
  theme_classic(base_size = BS, base_family = "Arial") +
  theme(
    axis.title.y.right        = element_text(colour = "#0072B2", size = AT,
                                             angle = -90, vjust = 0.5),
    axis.text.y.right         = element_text(colour = "#0072B2", size = ATX,
                                             margin = margin(l = 4)),
    axis.ticks.y.right        = element_line(colour = "#0072B2", linewidth = LW),
    axis.ticks.length.y.right = unit(4, "pt"),
    axis.line.y.right         = element_line(colour = "#0072B2", linewidth = LW),
    plot.margin               = margin(0, -2000, 0, -2000)
  )

p_bleach_donor <- ggplot(left_df_1985, aes(x = year, y = mean_bleaching_obs)) +
  geom_line(colour = "#D55E00", linewidth = LW) +
  scale_x_continuous(limits = c(1985, 2024), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10),
                     position = "right", expand = c(0, 0)) +
  labs(y = "Mean bleaching across sites (%)") +
  theme_classic(base_size = BS, base_family = "Arial") +
  theme(
    axis.title.y.right        = element_text(colour = "#D55E00", size = AT,
                                             angle = -90, vjust = 0.1),
    axis.text.y.right         = element_text(colour = "#D55E00", size = ATX,
                                             margin = margin(l = 4)),
    axis.ticks.y.right        = element_line(colour = "#D55E00", linewidth = LW),
    axis.ticks.length.y.right = unit(4, "pt"),
    axis.line.y.right         = element_line(colour = "#D55E00", linewidth = LW),
    plot.margin               = margin(0, -6000, 0, -2000)
  )

#### ── PANEL A MAIN ──────────────────────────────────────────────────────── ####

p_main <- ggplot(left_df_1985, aes(x = year)) +
  geom_line(aes(y = scales::rescale(World, to = c(0, 10), from = c(0, 2))),
            colour = "#009E73", linetype = "longdash", linewidth = LWM, alpha = 0.6) +
  geom_line(aes(y = scales::rescale(GMT, to = c(0, 10), from = c(GMT_LO, GMT_HI))),
            linetype = "dashed", colour = "grey40", alpha = 0.6, linewidth = LW) +
  geom_line(aes(y = mean_dhw), colour = "#0072B2", linewidth = LW) +
  geom_line(aes(y = scales::rescale(mean_bleaching_obs, to = c(0, 10), from = c(0, 50))),
            colour = "#D55E00", linewidth = LW, na.rm = TRUE) +
  scale_x_continuous(
    breaks = seq(1985, 2024, by = 5),
    limits = c(1985, 2024),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), expand = c(0, 0)) +
  labs(x = "Year", tag = "a") +
  theme_classic(base_size = BS, base_family = "Arial") +
  theme(
    axis.title.y.left  = element_blank(), axis.text.y.left  = element_blank(),
    axis.ticks.y.left  = element_blank(), axis.line.y.left  = element_blank(),
    axis.title.y.right = element_blank(), axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(),
    axis.title.x = element_text(size = AT),
    axis.text.x  = element_text(size = ATX, margin = margin(t = 3)),
    axis.line.x  = element_line(colour = "black", linewidth = LW),
    plot.margin  = margin(5, -80, 0, -20)
  )

p_a <- wrap_elements(get_plot_component(p_emissions_donor, "ylab-l")) +
  wrap_elements(get_y_axis(p_emissions_donor)) +
  wrap_elements(get_plot_component(p_gmst_donor, "ylab-l")) +
  wrap_elements(get_y_axis(p_gmst_donor)) +
  p_main +
  wrap_elements(p_spacer) +
  wrap_elements(get_y_axis(p_dhw_donor,    position = "right")) +
  wrap_elements(get_plot_component(p_dhw_donor,    "ylab-r")) +
  wrap_elements(p_spacer) +
  wrap_elements(get_y_axis(p_bleach_donor, position = "right")) +
  wrap_elements(get_plot_component(p_bleach_donor, "ylab-r")) +
  plot_layout(widths = c(1.2, 0.2, 0.5, 0.5, 22, 0, 0.1, 0.6, -0.1, 0.3, 1.1))

#### ── PANEL B ───────────────────────────────────────────────────────────── ####

PLOT_ABS_LAT <- 13
base_path    <- "path/to/Coracle/OneDrive"

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
    legend.key.size   = unit(0.5, "cm"),
    panel.grid        = element_blank(),
    plot.margin       = margin(4, 6, 4, 4),
    aspect.ratio = 1
  )

#### ── PANEL C ───────────────────────────────────────────────────────────── ####

All_Bleaching_Events_Data_AllDHW <- read.csv("path/to/Coracle/Datasets/Final_Panel_Data_allDHW_01July.csv")
All_Bleaching_Events_Data_AllDHW$abs_lat <- abs(All_Bleaching_Events_Data_AllDHW$Latitude_Degrees)
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Ecoregion_Month = interaction(Ecoregion_Name, Date_Month, drop = TRUE)
  )
model_linear_lat <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

dhw_max        <- max(All_Bleaching_Events_Data_AllDHW$dhw, na.rm = TRUE)
dhw_seq        <- seq(0, dhw_max, length.out = 300)
median_abs_lat <- median(All_Bleaching_Events_Data_AllDHW$abs_lat, na.rm = TRUE)
y_max          <- 100

beta_dhw     <- coef(model_linear_lat)["dhw"]
beta_dhw_lat <- coef(model_linear_lat)["dhw:abs_lat"]
se_dhw       <- model_linear_lat$se["dhw"]
se_dhw_lat   <- model_linear_lat$se["dhw:abs_lat"]
cov_dhw      <- vcov(model_linear_lat)["dhw", "dhw:abs_lat"]

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
  make_linear_curve(0,              "0\u00b0"),
  make_linear_curve(13, "13\u00b0"),
  make_linear_curve(23,             "23\u00b0")
) %>%
  mutate(lat_group = factor(lat_group, levels = c("0\u00b0", "13\u00b0", "23\u00b0")))

curve_colours <- c("13\u00b0" = "#000000", "0\u00b0" = "#FF4086", "23\u00b0" = "#35A1FF")
curve_lty     <- c("13\u00b0" = "solid",   "0\u00b0" = "dashed",  "23\u00b0" = "dashed")

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
    legend.key.width  = unit(0.8, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y  = unit(0.05, "cm"),
    panel.grid        = element_blank(),
    plot.margin       = margin(4, 6, 4, 4),
    aspect.ratio = 1
  )

#### ── ASSEMBLE & SAVE ───────────────────────────────────────────────────── ####
fig <- p_a / (p_binned | p_bl_nature) +
  plot_layout(heights = c(1.1, 1), widths = c(1))

fig


fig <- ggdraw() +
  draw_plot(as_grob(p_a),    x = 0.1,  y = 0.45, width = 0.8,  height = 0.54) +  # narrower A
  draw_plot(p_binned,        x = 0.09, y = 0,    width = 0.45, height = 0.46) +
  draw_plot(p_bl_nature,     x = 0.43, y = 0,    width = 0.45, height = 0.46)
fig




ggsave("figure1_diff_layout.png", fig,
       width  = 220,   # narrower overall
       height = 170,   # taller relative to width
       units  = "mm",
       dpi    = 200)




