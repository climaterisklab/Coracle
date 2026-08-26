#### author: Puja Pande
#### description: builds Supplementary Figure S1 — model response curves
####              (ZOIB/ordered beta and FE-binned DHW) with fixed-effects binned
####              estimates overlaid

library(dplyr)
library(fixest)
library(betareg)
library(ggplot2)
library(tidyr)
library(glmmTMB)
library(patchwork)



#### Prep data ####
All_Bleaching_Events_Data_AllDHW <- read.csv("path/to/Coracle/Datasets/Final_Panel_Data_allDHW_01July.csv")
All_Bleaching_Events_Data_AllDHW$Proportion_Bleached <- All_Bleaching_Events_Data_AllDHW$Percent_Bleached / 100
All_Bleaching_Events_Data_AllDHW$abs_lat             <- abs(All_Bleaching_Events_Data_AllDHW$Latitude_Degrees)
All_Bleaching_Events_Data_AllDHW$mass_bleaching      <- ifelse(All_Bleaching_Events_Data_AllDHW$Percent_Bleached >= 30, 1, 0)
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Ecoregion_Month = interaction(Ecoregion_Name, Date_Month, drop = TRUE)
  )
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Proportion_Bleached     = Percent_Bleached / 100,
    Proportion_Bleached_adj = (Proportion_Bleached * (n() - 1) + 0.5) / n(),
    any_bleaching           = ifelse(Percent_Bleached > 0, 1, 0)
  )


All_Bleaching_Events_Data_AllDHW_Max <- read.csv("path/to/data/Final_Combined_Data_all_levels_01July2026_extravars.csv")
All_Bleaching_Events_Data_AllDHW_Max <- left_join(All_Bleaching_Events_Data_AllDHW, All_Bleaching_Events_Data_AllDHW_Max[, c(1:8, 17:34)], by = c("Site_ID", "Latitude_Degrees", "Longitude_Degrees", "Date_Year", "Date_Month", "Date_Day", "Ecoregion_Name", "Percent_Bleached"))




#### Fit models ####
# Linear FE
model_linear <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW_Max
)


#### ZOIB (ordered beta) with random effects ####
model_zoib <- glmmTMB(
  Proportion_Bleached ~ dhw + dhw:abs_lat + (1 | Site_ID) + (1 | Date_Year) + (1 | Ecoregion_Month),
  family = ordbeta(),
  data   = All_Bleaching_Events_Data_AllDHW_Max
)

#### coefficients ####
dhw_seq <- seq(0, max(All_Bleaching_Events_Data_AllDHW_Max$dhw, na.rm = TRUE),
               length.out = 300)

# Reference latitude the curves below are evaluated at, matching figure_1_updated.R's
# panel C "13°" curve. Both models include a dhw:abs_lat interaction, so the slope
# isn't a single number — it depends on which latitude you evaluate it at (e.g. the
# linear model's slope is 3.82 at the equator vs 2.13 at 13°, since the interaction
# coefficient is negative). Previously these curves used only the dhw main-effect
# coefficient, silently plotting the abs_lat=0 (steepest) slope instead.
PLOT_ABS_LAT <- 13

# Linear
b_lin      <- coef(model_linear)["dhw"]
b_lin_lat  <- coef(model_linear)["dhw:abs_lat"]
se_lin     <- model_linear$se["dhw"]

# ZOIB
b_zoib     <- fixef(model_zoib)$cond["dhw"]
b_zoib_lat <- fixef(model_zoib)$cond["dhw:abs_lat"]
se_zoib    <- sqrt(vcov(model_zoib)$cond["dhw", "dhw"])
int_zoib   <- fixef(model_zoib)$cond["(Intercept)"]

#### Predictions ####
pred_df <- data.frame(dhw = dhw_seq) %>%
  mutate(
    # Linear, evaluated at abs_lat = PLOT_ABS_LAT
    slope_lin    = b_lin + b_lin_lat * PLOT_ABS_LAT,
    pred_linear  = slope_lin * dhw,
    ci_lo_linear = pmax((slope_lin - 1.96 * se_lin) * dhw, 0),
    ci_hi_linear = pmin((slope_lin + 1.96 * se_lin) * dhw, 100),

    # ZOIB, evaluated at abs_lat = PLOT_ABS_LAT - subtract baseline at DHW=0
    slope_zoib    = b_zoib + b_zoib_lat * PLOT_ABS_LAT,
    baseline_zoib = plogis(int_zoib) * 100,
    pred_zoib     = plogis(int_zoib + slope_zoib * dhw) * 100 - baseline_zoib,
    ci_lo_zoib    = pmax(plogis(int_zoib + (slope_zoib - 1.96 * se_zoib) * dhw) * 100 - baseline_zoib, 0),
    ci_hi_zoib    = pmin(plogis(int_zoib + (slope_zoib + 1.96 * se_zoib) * dhw) * 100 - baseline_zoib, 100)
  ) %>%
  select(-baseline_zoib, -slope_lin, -slope_zoib)

#### Reshape to long ####
pred_long <- pred_df %>%
  pivot_longer(starts_with("pred_"),
               names_to = "model", values_to = "pred") %>%
  mutate(model = recode(model,
                        pred_linear   = "Linear with Latitude Interaction",
                        pred_zoib     = "Zero-One-Inflated Beta with Latitude\nInteraction"
  ))

ci_long <- pred_df %>%
  pivot_longer(starts_with("ci_lo_"),
               names_to = "model", values_to = "ci_lo") %>%
  mutate(model = recode(model,
                        ci_lo_linear   = "Linear with Latitude Interaction",
                        ci_lo_zoib     = "Zero-One-Inflated Beta with Latitude\nInteraction"
  )) %>%
  left_join(
    pred_df %>%
      pivot_longer(starts_with("ci_hi_"),
                   names_to = "model", values_to = "ci_hi") %>%
      mutate(model = recode(model,
                            ci_hi_linear   = "Linear with Latitude Interaction",
                            ci_hi_zoib     = "Zero-One-Inflated Beta with Latitude\nInteraction"
      )),
    by = c("dhw", "model")
  )

plot_df <- left_join(pred_long,
                     ci_long %>% select(dhw, model, ci_lo, ci_hi),
                     by = c("dhw", "model"))

#### Colours ####
model_colours <- c(
  "Linear with Latitude Interaction"                   = "#0072B2",
  "Zero-One-Inflated Beta with Latitude\nInteraction"         = "#E69F00"
)

#### Plot ####
p_models <- ggplot() +
  geom_ribbon(
    data  = plot_df,
    aes(x = dhw, ymin = ci_lo, ymax = ci_hi, fill = model),
    alpha = 0.15
  ) +
  geom_line(
    data      = plot_df,
    aes(x = dhw, y = pred, colour = model),
    linewidth = 0.8
  ) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_fill_manual(values   = model_colours, guide = "none") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 9), breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) +
  labs(
    x = "Degree heating weeks (DHW)",
    y = "Coral bleaching (%)"
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line        = element_blank(),
    panel.border      = element_rect(fill = NA, colour = "black", linewidth = 0.5),
    axis.ticks        = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks.length = unit(3, "pt"),
    axis.title.y = element_text(size = 7, colour = "black", margin = margin(r = 10)),
    axis.text        = element_text(size = 6, colour = "black"),
    axis.title       = element_text(size = 7, colour = "black"),
    legend.text      = element_text(size = 6, colour = "black"),
    legend.position  = "top",
    legend.key.width = unit(0.5, "cm"),
    panel.grid       = element_blank(),
    plot.margin      = margin(4, 6, 4, 4),
    aspect.ratio = 0.7
  )

p_models

#### Settings ####
base_path    <- "path/to/Coracle"

#### Load model outputs ####
modlin      <- read.csv(file.path(base_path, "linear_lat_conley200km_coefs.csv"))
modbin      <- read.csv(file.path(base_path, "FE_binned_DHW_conley200km_coefs.csv"))
modlin_vcov <- readRDS(file.path(base_path, "linear_lat_conley200km_vcov.rds"))
modbin_vcov <- readRDS(file.path(base_path, "FE_binned_DHW_conley200km_vcov.rds"))


b_dhw     <- modlin$x[1]
b_dhw_lat <- modlin$x[2]

#### Parse bin terms ####
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

get_left_edge <- function(b) {
  as.numeric(sub("\\[([0-9.]+),.*", "\\1", b))
}

#### Binned estimates function ####
binned_estimates <- function(abs_lat) {
  bins <- unique(modbin$bin)
  
  result <- lapply(bins, function(bin_label) {
    main_name <- paste0('DHW_bin::', bin_label, '")')
    lat_name  <- paste0('DHW_bin::', bin_label, ':abs_lat")')
    
    pred <- main_terms[[bin_label]] + lat_terms[[bin_label]] * abs_lat
    
    var_main <- modbin_vcov[main_name, main_name]
    var_lat  <- modbin_vcov[lat_name,  lat_name]
    cov_ml   <- modbin_vcov[main_name, lat_name]
    var_pred <- var_main + abs_lat^2 * var_lat + 2 * abs_lat * cov_ml
    
    data.frame(
      dhw     = get_left_edge(bin_label),
      pred    = pred,
      se      = sqrt(var_pred),
      abs_lat = abs_lat
    )
  })
  
  bind_rows(
    data.frame(dhw = 0, pred = 0, se = 0, abs_lat = abs_lat),
    bind_rows(result)
  )
}

#### Generate FE-binned estimates ####
bin_df <- binned_estimates(PLOT_ABS_LAT) %>%
  mutate(
    ci_lo = pred - 1.96 * se,
    ci_hi = pred + 1.96 * se
  )

print(bin_df)

#### Model curves + FE-binned estimates overlaid ####
p_models_binned <- p_models +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black",
             linewidth  = 0.4) +
  geom_errorbar(
    data  = bin_df,
    aes(x = dhw, ymin = pmax(ci_lo, 0), ymax = ci_hi),
    width = 0.2, linewidth = 0.4, colour = "black"
  ) +
  geom_point(
    data   = bin_df,
    aes(x  = dhw, y = pred),
    colour = "black", size = 1.5, shape = 16
  )

p_models_binned

ggsave("supp_model_comparison.png", p_models_binned,
       width  = 89,
       height = 89,
       units  = "mm",
       dpi    = 600)


