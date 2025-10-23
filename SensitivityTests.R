### author: Puja Pande
### date: 06 October 2025
### description: the script runs a variety of sensitivity test to see how the results change with different parameters and assumptions

#### load libraries ####
library(dplyr)
library(fixest)
library(ggplot2)
library(tidyr)
library(patchwork)
library(plotly)

#### load datasets ####
#all_bleaching_events <- read.csv("All_Bleaching_Events_Data.csv")
all_bleaching_events_depth_filtered <- read.csv("All_Bleaching_Events_Data_Depth_Filtered.csv")
All_Bleaching_Events_Data_AllDHW <- read.csv("All_Bleaching_Events_Data_AllDHW.csv")
all_bleaching_events_depth_filtered <- All_Bleaching_Events_Data_AllDHW %>% filter(Depth_m <= 10)



#### different model types ####
model_linear <- feols(Percent_Bleached ~ dhw | Site_ID + Date_Year, cluster = ~Ecoregion_Name, data = All_Bleaching_Events_Data_AllDHW)
model_log <- feols(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, data = All_Bleaching_Events_Data_AllDHW)
model_poisson <- fepois(Percent_Bleached ~ dhw | Site_ID + Date_Year, cluster = ~Ecoregion_Name, data = All_Bleaching_Events_Data_AllDHW)
model_poisson_log <- fepois(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                            data = All_Bleaching_Events_Data_AllDHW)

summary(model_linear)
summary(model_log)
summary(model_poisson)
summary(model_poisson_log)

## initially used felm model and the coefficients matched with feols
# summary(lfe::felm(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year | 0 | Ecoregion_Name, data = all_bleaching_events))

# ## converting and comparing coefficients to each other
# (marginal_effect_logX <- coef(model_log)["log1p(dhw)"] / (1 + mean(all_bleaching_events$dhw, na.rm = TRUE)))
# (percent_change_poisson <- 100 * (exp(coef(model_poisson)) - 1))
# (percent_change_poisson_log <- 100 * coef(model_poisson_log)) ## best - continue with this model
# 
# -------
# # Choose a meaningful change in dhw
# delta_X <- 1  # or whatever is meaningful in your context
# mean_X <- mean(all_bleaching_events$dhw, na.rm = TRUE)
# 
# # Model 1: Linear
# effect_linear <- coef(model_linear)["dhw"] * delta_X
# 
# # Model 2: Poisson (% change in E[Y])
# effect_poisson <- 100 * (exp(coef(model_poisson)["dhw"] * delta_X) - 1)
# 
# # Model 3: Linear with log(X) - evaluate at mean
# effect_log <- coef(model_log)["log1p(dhw)"] * delta_X / (1 + mean_X)
# 
# # Model 4: Poisson with log(X) - for 1% increase
# effect_poisson_log <- coef(model_poisson_log)["log1p(dhw)"]  # This IS the elasticity
# # OR for a specific % increase:
# pct_increase <- 10  # 10% increase in X
# effect_poisson_log_10pct <- pct_increase * coef(model_poisson_log)["log1p(dhw)"]
# #A 10% increase in dhw is associated with a 5.3% increase in expected bleaching counts
# 


#### different interactions ####
all_bleaching_events_edited <- All_Bleaching_Events_Data_AllDHW %>% group_by(Site_ID) %>% 
  mutate(site_temp_ave = mean(ClimSST, na.rm = TRUE)) %>% ungroup() %>% mutate(site_temp_ave_int = site_temp_ave*dhw) %>% 
  ungroup()
all_bleaching_events_edited <- all_bleaching_events_edited %>% group_by(Site_ID) %>% 
  mutate(depth_int = Depth_m*dhw) %>% ungroup() 
all_bleaching_events_edited <- all_bleaching_events_edited %>% group_by(Site_ID) %>% 
  mutate(lat_int = Latitude_Degrees*dhw) %>% ungroup() 
all_bleaching_events_edited <- all_bleaching_events_edited %>% group_by(Site_ID) %>% 
  mutate(turbidity_int = Turbidity*dhw) %>% ungroup() 

model_site_temp_ave_int <- fepois(Percent_Bleached ~ dhw + site_temp_ave_int | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                  data = all_bleaching_events_edited)
model_depth_int <- fepois(Percent_Bleached ~ dhw + depth_int | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                  data = all_bleaching_events_edited)
model_latitude_int <- fepois(Percent_Bleached ~ dhw + lat_int | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                  data = all_bleaching_events_edited)
model_turbidity_int <- fepois(Percent_Bleached ~ dhw + turbidity_int | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                             data = all_bleaching_events_edited)

summary(model_site_temp_ave_int)     ### not significant
summary(model_depth_int)     ### not significant
summary(model_latitude_int)     ### not significant
summary(model_turbidity_int)  ### significant

#### cluster sensitivity test ####
model_cluster_none <- fepois(Percent_Bleached ~ dhw | Site_ID + Date_Year, data = All_Bleaching_Events_Data_AllDHW)
model_cluster_ecoregion <- fepois(Percent_Bleached ~ dhw | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                  data = All_Bleaching_Events_Data_AllDHW)
model_cluster_site <- fepois(Percent_Bleached ~ dhw | Site_ID + Date_Year, cluster = ~Site_ID,
                                  data = All_Bleaching_Events_Data_AllDHW)
summary(model_cluster_none)        
summary(model_cluster_ecoregion)   
summary(model_cluster_site)        


#### depth sensitivity test ####
model_all_depths <- fepois(Percent_Bleached ~ dhw | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                          data = All_Bleaching_Events_Data_AllDHW)
model_depth_10andless <- fepois(Percent_Bleached ~ dhw | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                data = all_bleaching_events_depth_filtered)
summary(model_all_depths)
summary(model_depth_10andless)           # 0.002 difference between the two models, both significant

## sites with multiple measurements at the same depth
depth_same <- All_Bleaching_Events_Data_AllDHW %>%
  group_by(Site_ID, Depth_m) %>%
  filter(n() > 1) %>%
  ungroup()

model_same_depth <- fepois(Percent_Bleached ~ dhw | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                               data = depth_same)
summary(model_same_depth)  



#### lag sensitivity test ####
model_lagminus1 <- fepois(Percent_Bleached ~ dhw_lag.1 | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                 data = All_Bleaching_Events_Data_AllDHW)
model_lag0 <- fepois(Percent_Bleached ~ dhw_lag0 | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                     data = All_Bleaching_Events_Data_AllDHW)
model_lag1 <- fepois(Percent_Bleached ~ dhw_lag1 | Site_ID + Date_Year, cluster = ~Ecoregion_Name,
                     data = All_Bleaching_Events_Data_AllDHW)
model_lag2 <- fepois(Percent_Bleached ~ dhw_lag2 | Site_ID + Date_Year, cluster = ~Ecoregion_Name,
                     data = All_Bleaching_Events_Data_AllDHW)
model_lag3 <- fepois(Percent_Bleached ~ dhw_lag3 | Site_ID + Date_Year, cluster = ~Ecoregion_Name,
                     data = All_Bleaching_Events_Data_AllDHW)
model_lag4 <- fepois(Percent_Bleached ~ dhw_lag4 | Site_ID + Date_Year, cluster = ~Ecoregion_Name,
                     data = All_Bleaching_Events_Data_AllDHW)
model_lag5 <- fepois(Percent_Bleached ~ dhw_lag5 | Site_ID + Date_Year, cluster = ~Ecoregion_Name,
                     data = All_Bleaching_Events_Data_AllDHW)
summary(model_lagminus1)   # significant
summary(model_lag0)        # significant           ## best lag - 0 = lag 1 as by Chris, normal lag as per NOAA
summary(model_lag1)        # significant
summary(model_lag2)        # significant
summary(model_lag3)        # not significant
summary(model_lag4)        # not significant
summary(model_lag5)        # not significant

## max lag
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>% rowwise() %>%
  mutate(dhw_val_max = max(c_across(starts_with("dhw_lag")), na.rm = TRUE)) %>%
  ungroup()
model_lag_max <- fixest::fepois(Percent_Bleached ~ dhw_val_max | Site_ID + Date_Year, 
                                cluster = ~Ecoregion_Name,
                                data = All_Bleaching_Events_Data_AllDHW)
summary(model_lag_max) ## significant but 0 lag better


#### generating some more sensitivity tests ####
## using negative binomial instead of poisson
model_negative_binomial <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                    data = All_Bleaching_Events_Data_AllDHW)
summary(model_negative_binomial)         ## ver similar to poisson



#### rounding off percent bleaching to make it a proper count variable ####
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Percent_Bleached_Counts = case_when(
      Percent_Bleached > 0 & Percent_Bleached < 1 ~ 1,
      TRUE ~ round(Percent_Bleached)
    )
  )

final_model_counts <- fepois(Percent_Bleached_Counts ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                      data = All_Bleaching_Events_Data_AllDHW)

## final model 
final_model <- fepois(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                      data = All_Bleaching_Events_Data_AllDHW)

#### model summary ####
library(modelsummary)

# --- Comprehensive Model Summary ---
modelsummary(
  list(
    # ---- Base Models ----
    "Linear"                     = model_linear,
    "Log(DHW)"                   = model_log,
    "Poisson"                    = model_poisson,
    "Poisson Log(DHW)"           = model_poisson_log,
    "Negative Binomial (Log DHW)"= model_negative_binomial,
    
    # ---- Interaction Models ----
    "Temp × DHW"                 = model_site_temp_ave_int,
    "Depth × DHW"                = model_depth_int,
    "Latitude × DHW"             = model_latitude_int,
    "Turbidity × DHW"            = model_turbidity_int,
    
    # ---- Depth Sensitivity ----
    "All Depths"                 = model_all_depths,
    "Depth ≤ 10m"                = model_depth_10andless,
    "Same Depth Sites"           = model_same_depth,
    
    # ---- Lag Sensitivity ----
    "Lag -1 (DHW_lag.1)"         = model_lagminus1,
    "Lag 0 (DHW_lag0)"           = model_lag0,
    "Lag +1 (DHW_lag1)"          = model_lag1,
    "Lag +2 (DHW_lag2)"          = model_lag2,
    "Lag +3 (DHW_lag3)"          = model_lag3,
    "Lag +4 (DHW_lag4)"          = model_lag4,
    "Lag +5 (DHW_lag5)"          = model_lag5,
    "Max Lag"                    = model_lag_max,
    
    # ---- Clustering Sensitivity ----
    "Cluster: None"              = model_cluster_none,
    "Cluster: Ecoregion"         = model_cluster_ecoregion,
    "Cluster: Site"              = model_cluster_site,
    
    # ---- Final Models ----
    "Final Model (Percent Bleached)" = final_model,
    "Final Model (Counts)"           = final_model_counts
  ),
  
  # --- Display Options ---
  stars = c('*' = 0.1, '**' = 0.05, '***' = 0.01),
  coef_rename = c(
    "dhw"              = "DHW",
    "log1p(dhw)"       = "log(DHW + 1)",
    "site_temp_ave_int"= "Site Temp × DHW",
    "depth_int"        = "Depth × DHW",
    "lat_int"          = "Latitude × DHW",
    "turbidity_int"    = "Turbidity × DHW",
    "dhw_lag.1"        = "DHW (Lag -1)",
    "dhw_lag0"         = "DHW (Lag 0)",
    "dhw_lag1"         = "DHW (Lag +1)",
    "dhw_lag2"         = "DHW (Lag +2)",
    "dhw_lag3"         = "DHW (Lag +3)",
    "dhw_lag4"         = "DHW (Lag +4)",
    "dhw_lag5"         = "DHW (Lag +5)",
    "dhw_val_max"      = "DHW (Max)"
  ),
  statistic = NULL,  # omit SEs, z-stats, and p-values for compactness
  gof_omit = 'AIC|Log.Lik|F|RMSE|R2|R2 Adj.|R2 Within|R2 Within Adj.|BIC|Std.Errors',
  output = "model_summary_all.html"
)


#### plotting the main model ####
# First, identify which rows were used (no NAs in key variables)
model_data <- All_Bleaching_Events_Data_AllDHW %>%
  filter(!is.na(Percent_Bleached) & !is.na(dhw)) %>%
  # Add predictions - fixest will only predict for rows it can use
  mutate(predicted_bleaching = predict(final_model, newdata = ., type = "response"))

# Remove any remaining NAs from predictions (from dropped fixed effects)
model_data <- model_data %>%
  filter(!is.na(predicted_bleaching))

# Now plot
ggplot(model_data, aes(x = dhw)) +
  geom_point(aes(y = Percent_Bleached), alpha = 0.3, color = "gray50", size = 1) +
  geom_smooth(aes(y = predicted_bleaching), 
              method = "loess", color = "steelblue", size = 1.2, se = TRUE) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Percent Bleaching vs DHW",
    subtitle = "Points show actual observations, blue line shows model predictions"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12)
  )


#### Comparative model visualization ####
# Base data for prediction (no missing DHW or bleaching)
model_data <- All_Bleaching_Events_Data_AllDHW %>%
  filter(!is.na(Percent_Bleached), !is.na(dhw))

## did not add negatuve binomial no log because percent bleaching goes over 100

# Add predictions from all relevant models
model_data <- model_data %>%
  mutate(
    pred_linear                = predict(model_linear, newdata = ., type = "response"),
    pred_log                   = predict(model_log, newdata = ., type = "response"),
    pred_poisson               = predict(model_poisson, newdata = ., type = "response"),
    pred_poisson_log           = predict(model_poisson_log, newdata = ., type = "response"),
    pred_negative_binomial     = predict(model_negative_binomial, newdata = ., type = "response")
  )

# Filter out rows with missing predictions (common for FE models)
model_data <- model_data %>%
  filter(
    !is.na(pred_linear),
    !is.na(pred_log),
    !is.na(pred_poisson),
    !is.na(pred_poisson_log),
    !is.na(pred_negative_binomial)
  )

# Reshape to long format for ggplot
plot_data <- model_data %>%
  dplyr::select(
    dhw, Percent_Bleached,
    pred_linear, pred_log,
    pred_poisson, pred_poisson_log,
    pred_negative_binomial
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("pred_"),
    names_to = "model",
    values_to = "predicted_bleaching"
  ) %>%
  dplyr::mutate(model = dplyr::recode(
    model,
    pred_linear = "Linear",
    pred_log = "Log(DHW)",
    pred_poisson = "Poisson",
    pred_poisson_log = "Poisson Log(DHW)",
    pred_negative_binomial = "Negative Binomial (Log DHW)"
    ))


# Define custom color palette for publication-style clarity
model_colors <- c(
  "Linear" = "#1b9e77",
  "Log(DHW)" = "blue",
  "Poisson" = "purple",
  "Poisson Log(DHW)" = "red",
  "Negative Binomial (Log DHW)" = "orange"
)

# Generate plot
ggplot(plot_data, aes(x = dhw)) +
  geom_point(aes(y = Percent_Bleached),
             alpha = 0.25, color = "gray50", size = 0.8) +
  geom_smooth(aes(y = predicted_bleaching, color = model),
              method = "loess", se = TRUE, size = 1.1, span = 0.7) +
  scale_color_manual(values = model_colors) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Percent Bleaching vs DHW",
    subtitle = "Observed bleaching (gray points) with fitted curves from multiple model specifications",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )



#### ecoregion curves ####
#### STEP 0: Define regional groups ####

regions_caribbean <- c(
  "Hispaniola, Puerto Rico and Lesser Antilles",
  "Bahamas and Florida Keys",
  "Belize and west Caribbean",
  "Netherlands Antilles and south Caribbean",
  "Cuba and Cayman Islands",
  "Jamaica"
)

regions_gbr <- c(
  "Central and northern Great Barrier Reef",
  "Southern Great Barrier Reef",
  "Coral Sea"
)

regions_coral_triangle <- c(
  "Sunda Shelf, south-east Asia", "Gulf of Thailand", "Lesser Sunda Islands and Savu Sea",
  "Java Sea", "Banda Sea and Molucca Islands", "Gulf of Tomini, Indonesia",
  "Makassar Strait, Indonesia", "West Sumatra", "South Java",
  "Cenderawasih Bay, Papua", "Birds Head Peninsula, Papua", "Celebes Sea",
  "Sulu Sea", "South-east Philippines", "North Philippines",
  "Solomon Islands and Bougainville", "Milne Bay, Papua New Guinea", "Bismarck Sea, New Guinea"
)

regions_east_africa <- c(
  "Kenya and Tanzania coast", "North Mozambique coast",
  "North Madagascar", "South Madagascar", "Mayotte and Comoros"
)

regions_maldives <- c("Maldive Islands")

regions_focus <- c(
  regions_caribbean,
  regions_gbr,
  regions_coral_triangle,
  regions_east_africa,
  regions_maldives
)

#### STEP 1: Prepare data ####

model_data <- All_Bleaching_Events_Data_AllDHW %>%
  filter(
    !is.na(Percent_Bleached),
    !is.na(dhw),
    Ecoregion_Name %in% regions_focus
  ) %>%
  mutate(pred_poisson_log = predict(model_poisson_log, newdata = ., type = "response")) %>%
  filter(!is.na(pred_poisson_log)) 

#### STEP 2: Ecoregion summary ####

ecoregion_stats <- model_data %>%
  group_by(Ecoregion_Name) %>%
  summarise(
    n_obs = n(),
    n_sites = n_distinct(Site_ID),
    mean_dhw = mean(dhw, na.rm = TRUE),
    mean_bleaching = mean(Percent_Bleached, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_obs))

# ecoregion_stats_filtered <- ecoregion_stats %>%
#   filter(n_obs > 100)
# 
# 
# selected_ecoregions <- ecoregion_stats_filtered$Ecoregion_Name
# 
# model_data <- model_data %>%
#   filter(Ecoregion_Name %in% selected_ecoregions)


#### STEP 3: Faceted plot for focus regions ####

p_faceted <- ggplot(model_data, aes(x = dhw)) +
  geom_point(aes(y = Percent_Bleached),
             alpha = 0.4, color = "gray50", size = 0.8) +
  geom_smooth(aes(y = pred_poisson_log),
              method = "loess", color = "red", size = 1, se = TRUE, span = 0.8) +
  facet_wrap(~ Ecoregion_Name, scales = "free_y", ncol = 4) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Coral Bleaching by Ecoregion (Poisson Log Model)",
    subtitle = "Gray points = observed bleaching | Red line = model predictions with 95% CI"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  )

print(p_faceted)


#### STEP 5: Marginal effect comparison ####

ecoregion_slopes <- model_data %>%
  group_by(Ecoregion_Name) %>%
  summarise(
    n_obs = n(),
    mean_dhw = mean(dhw, na.rm = TRUE),
    mean_pred = mean(pred_poisson_log, na.rm = TRUE),
    beta = coef(model_poisson_log)["log1p(dhw)"],
    marginal_effect = beta * mean_pred / (1 + mean_dhw),
    .groups = "drop"
  ) %>%
  arrange(desc(marginal_effect))

p_slopes <- ggplot(ecoregion_slopes,
                   aes(x = reorder(Ecoregion_Name, marginal_effect),
                       y = marginal_effect, fill = marginal_effect)) +
  geom_col(alpha = 0.8) +
  coord_flip() +
  scale_fill_gradient2(
    low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c",
    midpoint = median(ecoregion_slopes$marginal_effect),
    name = "Marginal\nEffect"
  ) +
  labs(
    title = "Bleaching Sensitivity by Ecoregion (Focus Regions)",
    subtitle = "Marginal effect: % bleaching increase per 1 DHW increase",
    x = NULL,
    y = "Marginal Effect (% bleaching per DHW)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14))

print(p_slopes)



#### demeaned data and plots ####
# ===== 1. Create Demeaned Dataset =====
All_Bleaching_Events_Data_Demeaned <- All_Bleaching_Events_Data_AllDHW %>%
  group_by(Site_ID, Date_Year) %>%
  mutate(
    mean_log_dhw = mean(log1p(dhw), na.rm = TRUE),
    log_dhw_demeaned = log1p(dhw) - mean_log_dhw
  ) %>%
  ungroup()

# ===== 2. Fit Both Models =====
original_model <- fepois(
  Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)

demeaned_model <- fepois(
  Percent_Bleached ~ log_dhw_demeaned,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_Demeaned
)

# ===== 3. Generate Predictions =====
All_Bleaching_Events_Data_AllDHW$pred_original <- predict(
  original_model, 
  All_Bleaching_Events_Data_AllDHW
)

All_Bleaching_Events_Data_Demeaned$pred_demeaned <- predict(
  demeaned_model,
  All_Bleaching_Events_Data_Demeaned
)

# Prepare data for plotting
plot_data_fitted <- bind_rows(
  All_Bleaching_Events_Data_AllDHW %>%
    select(dhw, Percent_Bleached, predicted_bleaching = pred_original) %>%
    mutate(model = "Original (with FE)"),
  All_Bleaching_Events_Data_Demeaned %>%
    select(dhw, Percent_Bleached, predicted_bleaching = pred_demeaned) %>%
    mutate(model = "Demeaned")
)

# ===== 4. Create Plot =====
model_colors <- c(
  "Original (with FE)" = "#1b9e77",
  "Demeaned" = "#d95f02"
)

# Separate the data by model
plot_data_original <- plot_data_fitted %>% 
  filter(model == "Original (with FE)")

plot_data_demeaned <- plot_data_fitted %>% 
  filter(model == "Demeaned")

ggplot(plot_data_fitted, aes(x = dhw)) +
  geom_point(
    aes(y = Percent_Bleached),
    alpha = 0.25, 
    color = "gray50", 
    size = 0.8
  ) +
  # Original model - NO confidence interval
  geom_smooth(
    data = plot_data_original,
    aes(y = predicted_bleaching, color = model),
    method = "loess", 
    se = TRUE,  # CI for original model
    linewidth = 1.1, 
    span = 0.7, 
    alpha = 0.2
  ) +
  # Demeaned model - WITH confidence interval
  geom_smooth(
    data = plot_data_demeaned,
    aes(y = predicted_bleaching, color = model),
    method = "loess", 
    se = TRUE,  # CI for demeaned model
    linewidth = 1.1, 
    span = 0.7,
    alpha = 0.2
  ) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Percent Bleaching vs DHW: Original vs Demeaned Model",
    subtitle = "Observed bleaching (gray points) with fitted curves from both model specifications",
    color = "Model Type",
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

