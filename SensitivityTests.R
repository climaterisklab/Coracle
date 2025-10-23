### author: Puja Pande
### date: 06 October 2025
### description: the script runs a variety of sensitivity test to see how the results change with different parameters and assumptions

#### load libraries ####
library(dplyr)
library(fixest)
library(ggplot2)

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
## final model 
final_model <- fepois(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                      data = All_Bleaching_Events_Data_AllDHW)

# Get the coefficient and SE for log1p(dhw)
coef_dhw <- coef(final_model)["log1p(dhw)"]
se_dhw <- se(final_model)["log1p(dhw)"]

# Generate DHW range
dhw_range <- seq(0, max(All_Bleaching_Events_Data_AllDHW$dhw, na.rm = TRUE), length.out = 100)

# Calculate linear predictor and CI on log scale
log_dhw <- log1p(dhw_range)
linear_pred <- coef_dhw * log_dhw
linear_lower <- (coef_dhw - 1.96 * se_dhw) * log_dhw
linear_upper <- (coef_dhw + 1.96 * se_dhw) * log_dhw

# Transform to response scale (exp for Poisson)
# Note: This is relative to baseline (fixed effects at 0)
pred_df <- data.frame(
  dhw = dhw_range,
  relative_risk = exp(linear_pred),
  lower = exp(linear_lower),
  upper = exp(linear_upper)
)

# Plot relative risk
ggplot(pred_df, aes(x = dhw)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.3, fill = "steelblue") +
  geom_line(aes(y = relative_risk),
            color = "steelblue", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Relative Risk of Bleaching",
    title = "Effect of DHW on Coral Bleaching",
    subtitle = "Relative to baseline (DHW = 0) with 95% CI"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12)
  )


#### some more plots - predicted percent bleaching vs DHW #### 
# Get predictions and merge back with the original data
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































library(fixest)
library(dplyr)
library(ggplot2)

# Define a sequence of predictor values
SSTA_seq <- seq(min(all_bleaching_events$dhw, na.rm = TRUE),
                max(all_bleaching_events$dhw, na.rm = TRUE),
                length.out = 200)

# Example: use one existing Site_ID and Date_Year from your data
example_site <- all_bleaching_events$Site_ID[1]
example_year <- all_bleaching_events$Date_Year[1]

newdata <- data.frame(
  dhw = SSTA_seq,
  Site_ID = example_site,
  Date_Year = example_year
)

# Then predict normally
newdata$linear <- predict(model_linear, newdata, type = "response")
newdata$log <- predict(model_log, newdata, type = "response")
newdata$poisson <- predict(model_poisson, newdata, type = "response")
newdata$poisson_log <- predict(model_poisson_log, newdata, type = "response")
newdata$neg_binom <- predict(model_negative_binomial, newdata, type = "response")

# --- Convert to long format for plotting ---
plot_df <- newdata %>%
  tidyr::pivot_longer(cols = c(linear, log, poisson, poisson_log, neg_binom),
                      names_to = "Model",
                      values_to = "Predicted_Bleaching")

# --- Plot ---
ggplot(plot_df, aes(x = dhw, y = Predicted_Bleaching, color = Model)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("linear" = "#1b9e77",
                                "log" = "#d95f02",
                                "poisson" = "#7570b3",
                                "poisson_log" = "#e7298a",
                                "neg_binom" = "#66a61e")) +
  labs(title = "Comparison of Model Fits: Linear, Log, Poisson, and Negative Binomial",
       x = "Sea Surface Temperature Anomaly (dhw)",
       y = "Predicted Percent Bleached",
       color = "Model Type") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")









library(fixest)
library(dplyr)
library(ggplot2)
library(tidyr)

# Define a sequence of predictor values
SSTA_seq <- seq(min(all_bleaching_events$dhw, na.rm = TRUE),
                max(all_bleaching_events$dhw, na.rm = TRUE),
                length.out = 200)

# Use a representative site and year from your data
example_site <- all_bleaching_events$Site_ID[1]
example_year <- all_bleaching_events$Date_Year[1]

# Create newdata with all required variables
newdata <- data.frame(
  dhw = SSTA_seq,
  Site_ID = example_site,
  Date_Year = example_year
)

# Get predictions from each model
# For fixest models, we need to handle fixed effects properly
newdata$log <- predict(model_log, newdata = newdata, type = "response")
newdata$poisson <- predict(model_poisson, newdata = newdata, type = "response")
newdata$poisson_log <- predict(model_poisson_log, newdata = newdata, type = "response")
newdata$neg_binom <- predict(model_negative_binomial, newdata = newdata, type = "response")

# Convert to long format for plotting - ALL FOUR MODELS
plot_df <- newdata %>%
  pivot_longer(cols = c(log, poisson, poisson_log, neg_binom),
               names_to = "Model",
               values_to = "Predicted_Bleaching")

# Create the plot with all four models
ggplot(plot_df, aes(x = dhw, y = Predicted_Bleaching, color = Model)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("log" = "#d95f02",
               "poisson" = "#7570b3",
               "poisson_log" = "#e7298a",
               "neg_binom" = "#66a61e"),
    labels = c("log" = "Log-Linear (feols, log1p(X))",
               "poisson" = "Poisson (fepois, X)",
               "poisson_log" = "Poisson (fepois, log1p(X))",
               "neg_binom" = "Negative Binomial (fenegbin, log1p(X))")
  ) +
  labs(
    title = "Comparison of Model Predictions: Four Model Specifications",
    x = "Sea Surface Temperature Anomaly (SSTA DHW, °C-weeks)",
    y = "Predicted Percent Bleached (%)",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10)
  )
