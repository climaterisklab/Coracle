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
library(betareg)
library(glmmTMB)
library(sandwich)
library(lmtest)
library(fixest)

#### load datasets ####
#all_bleaching_events <- read.csv("All_Bleaching_Events_Data.csv")
#all_bleaching_events_depth_filtered <- read.csv("All_Bleaching_Events_Data_Depth_Filtered.csv")
All_Bleaching_Events_Data_AllDHW <- read.csv("Final_Relevant_Scripts/All_Bleaching_Events_Data_AllDHW.csv")
All_Bleaching_Events_Data_AllDHW <- read.csv("Coracle_backup/Final_Relevant_Scripts/Merged_Mermaid_Panel_Bleaching_Data.csv")

All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Realm_Name = case_when(
      Realm_Name == "Eastern Indo-Pacific" ~ "Tropical Eastern Pacific",
      Realm_Name == "Temperate Australasia" ~ "Central Indo-Pacific",
      TRUE ~ Realm_Name  # Keep all other values as is
    )
  )


#Mermaid <- read.csv("Final_Relevant_Scripts/MermaidBleachingData_WithDHW.csv")
#all_bleaching_events_depth_filtered <- All_Bleaching_Events_Data_AllDHW %>% filter(Depth_m <= 10)



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
summary(model_negative_binomial)         ## very similar to poisson

## beta models
# Fit beta model with site and year fixed effects
model_beta_fe <- betareg(
  Proportion_Bleached_Trans ~ log1p(dhw) + Site_ID_factor + Date_Year_factor,
  data = model_data
)
# Fit ZOIB model (handles true 0s and 1s, no transformation needed)
model_zoib <- glmmTMB(Proportion_Bleached ~ log1p(dhw) + (1 | Site_ID) + (1 | Date_Year),
                      family = ordbeta(),
                      data = model_data)

## logit model
# Fit logit model with intercept (no fixed effects for comparison)
model_logit_intercept <- feglm(
  Proportion_Bleached ~ 1 + log1p(dhw),
  family = binomial(link = "logit"),
  cluster = ~Ecoregion_Name,
  data = model_data
)


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

## final model - old 
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


#----------------------------------------------------------------------------------
#sensitivity checks from 20/02/2026

#hausmann test on linear model
# Your current fixed effects model
fe_model <- feols(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, 
                  cluster = ~Ecoregion_Name,
                  data = All_Bleaching_Events_Data_AllDHW)

# Random effects approximation (no fixed effects)
re_model <- feols(Percent_Bleached ~ log1p(dhw), 
                  cluster = ~Ecoregion_Name,
                  data = All_Bleaching_Events_Data_AllDHW)

# Compare models
etable(fe_model, re_model)

# Or use AIC/BIC
AIC(fe_model, re_model)
BIC(fe_model, re_model)


library(plm)
library(lmtest)

# Aggregate to one observation per site-year (annual mean)
annual_data <- All_Bleaching_Events_Data_AllDHW %>%
  group_by(Site_ID, Date_Year, Ecoregion_Name) %>%
  summarise(
    Percent_Bleached = mean(Percent_Bleached, na.rm = TRUE),
    dhw = mean(dhw, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(Percent_Bleached), !is.na(dhw))

# Fixed effects model
fixed <- plm(Percent_Bleached ~ log1p(dhw), 
             data = annual_data, 
             index = c("Site_ID", "Date_Year"), 
             model = "within")

# Random effects model
random <- plm(Percent_Bleached ~ log1p(dhw), 
              data = annual_data, 
              index = c("Site_ID", "Date_Year"), 
              model = "random")

# Hausman test
phtest(fixed, random)


# bleaching dependence test -> lag the maximum bleaching observed at a site in the previous year and see if it predicts bleaching in the current year
#### Create lagged maximum bleaching variable ####
# Step 1: Calculate maximum bleaching per site-year
annual_max_bleaching <- All_Bleaching_Events_Data_AllDHW %>%
  group_by(Site_ID, Date_Year) %>%
  summarise(
    Max_Bleaching = max(Percent_Bleached, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Create lagged variable (previous year's max)
annual_max_bleaching <- annual_max_bleaching %>%
  group_by(Site_ID) %>%
  arrange(Site_ID, Date_Year) %>%
  mutate(
    Lag_Max_Bleaching = lag(Max_Bleaching, n = 1)
  ) %>%
  ungroup()

annual_max_bleaching <- annual_max_bleaching %>% filter(!is.infinite(Lag_Max_Bleaching))  # Keep only rows with valid lag


# Step 3: Merge back to original data
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  left_join(
    annual_max_bleaching %>% select(Site_ID, Date_Year, Lag_Max_Bleaching),
    by = c("Site_ID", "Date_Year")
  )


# Model 1: Without lag (your baseline)
model_no_lag <- fepois(
  Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)

# Model 2: With lag (tests bleaching dependence)
model_with_lag <- fepois(
  Percent_Bleached ~ log1p(dhw) + Lag_Max_Bleaching | Site_ID + Date_Year,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)

# Test if lag is significant
summary(model_with_lag)
summary(model_no_lag)


#### Interpretation ####

# Extract lag coefficient
lag_coef <- coef(model_with_lag)["Lag_Max_Bleaching"]
lag_se <- se(model_with_lag)["Lag_Max_Bleaching"]
lag_pval <- pvalue(model_with_lag)["Lag_Max_Bleaching"]

cat("\n=== Bleaching Dependence Test ===\n")
cat("Lag coefficient:", round(lag_coef, 4), "\n")
cat("Standard error:", round(lag_se, 4), "\n")
cat("P-value:", round(lag_pval, 4), "\n\n")

if(lag_pval < 0.05) {
  if(lag_coef > 0) {
    cat("Result: POSITIVE dependence - sites with high bleaching last year\n")
    cat("        bleach MORE this year (sensitization/damage)\n")
  } else {
    cat("Result: NEGATIVE dependence - sites with high bleaching last year\n")
    cat("        bleach LESS this year (adaptation/mortality selection)\n")
  }
} else {
  cat("Result: NO significant dependence - previous bleaching doesn't\n")
  cat("        predict current bleaching (independent events)\n")
}


#### Alternative: Test multiple lags ####

# Try lags up to 3 years
annual_max_bleaching <- annual_max_bleaching %>%
  group_by(Site_ID) %>%
  arrange(Site_ID, Date_Year) %>%
  mutate(
    Lag1_Max_Bleaching = lag(Max_Bleaching, n = 1),
    Lag2_Max_Bleaching = lag(Max_Bleaching, n = 2),
    Lag3_Max_Bleaching = lag(Max_Bleaching, n = 3)
  ) %>%
  ungroup()

All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  select(-starts_with("Lag")) %>%
  left_join(
    annual_max_bleaching %>% select(Site_ID, Date_Year, starts_with("Lag")),
    by = c("Site_ID", "Date_Year")
  )

# Model with multiple lags
model_multi_lag <- fepois(
  Percent_Bleached ~ log1p(dhw) + Lag1_Max_Bleaching + Lag2_Max_Bleaching + Lag3_Max_Bleaching | Site_ID + Date_Year,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)

etable(model_no_lag, model_with_lag, model_multi_lag)


#### Check data ####
cat("\n=== Data Summary ===\n")
cat("Total observations:", nrow(All_Bleaching_Events_Data_AllDHW), "\n")
cat("With lag data:", sum(!is.na(All_Bleaching_Events_Data_AllDHW$Lag_Max_Bleaching)), "\n")
cat("Missing lag:", sum(is.na(All_Bleaching_Events_Data_AllDHW$Lag_Max_Bleaching)), "\n\n")

# Show example of lag construction
cat("Example for one site:\n")
All_Bleaching_Events_Data_AllDHW %>%
  filter(Site_ID == first(Site_ID)) %>%
  select(Site_ID, Date_Year, Date_Month, Percent_Bleached, Lag_Max_Bleaching) %>%
  arrange(Date_Year, Date_Month) %>%
  head(20)




#----------------------------------------------------------------------------------

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

n_obs <- nrow(All_Bleaching_Events_Data_AllDHW %>% filter(!is.na(Percent_Bleached), !is.na(dhw)))
model_data <- All_Bleaching_Events_Data_AllDHW %>%
  filter(!is.na(Percent_Bleached), !is.na(dhw)) %>%
  mutate(
    Proportion_Bleached = Percent_Bleached / 100,
    Proportion_Bleached_Trans = (Proportion_Bleached * (n_obs - 1) + 0.5) / n_obs,
    Site_ID_factor = factor(Site_ID),
    Date_Year_factor = factor(Date_Year)
  )

# Fit beta model with site and year fixed effects
model_beta_fe <- betareg(
  Proportion_Bleached_Trans ~ log1p(dhw) + Site_ID_factor + Date_Year_factor,
  data = model_data
)

# Fit ZOIB model (handles true 0s and 1s, no transformation needed)
model_zoib <- glmmTMB(Proportion_Bleached ~ log1p(dhw) + (1 | Site_ID) + (1 | Date_Year),
                      family = ordbeta(),
                      data = model_data)

# Fit logit model with intercept (no fixed effects for comparison)
model_logit_intercept <- feglm(
  Proportion_Bleached ~ 1 + log1p(dhw),
  family = binomial(link = "logit"),
  cluster = ~Ecoregion_Name,
  data = model_data
)
cat("Logit model (with intercept) fitted successfully!\n")
cat("Intercept:", coef(model_logit_intercept)["(Intercept)"], "\n")
cat("Coefficient (log1p(dhw)):", coef(model_logit_intercept)["log1p(dhw)"], "\n\n")


#### Comparative model visualization ####
model_data <- model_data %>%
  mutate(
    pred_linear            = predict(model_linear, newdata = ., type = "response"),
    pred_log               = predict(model_log, newdata = ., type = "response"),
    pred_poisson           = predict(model_poisson, newdata = ., type = "response"),
    pred_poisson_log       = predict(model_poisson_log, newdata = ., type = "response"),
    pred_negative_binomial = predict(model_negative_binomial, newdata = ., type = "response"),
    pred_beta_fe           = predict(model_beta_fe, newdata = ., type = "response") * 100,
    pred_zoib              = predict(model_zoib, newdata = ., type = "response") * 100,
    pred_logit             = predict(model_logit_intercept, newdata = ., type = "response") * 100
  )

# Filter out rows with missing predictions
model_data <- model_data %>%
  filter(
    !is.na(pred_linear),
    !is.na(pred_log),
    !is.na(pred_poisson),
    !is.na(pred_poisson_log),
    !is.na(pred_negative_binomial),
    !is.na(pred_beta_fe),
    !is.na(pred_zoib),
    !is.na(pred_logit)
  )

# Reshape to long format
plot_data <- model_data %>%
  dplyr::select(
    dhw, Percent_Bleached,
    pred_linear, pred_log,
    pred_poisson, pred_poisson_log,
    pred_negative_binomial,
    pred_beta_fe,
    pred_zoib,
    pred_logit
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("pred_"),
    names_to = "model",
    values_to = "predicted_bleaching"
  ) %>%
  dplyr::mutate(model = dplyr::recode(
    model,
    pred_linear                = "Linear",
    pred_log                   = "Log(DHW)",
    pred_poisson               = "Poisson",
    pred_poisson_log           = "Poisson Log(DHW)",
    pred_negative_binomial     = "Negative Binomial (Log DHW)",
    pred_beta_fe               = "Beta Regression (FE)",
    pred_zoib                  = "ZOIB Ordbeta (RE)",
    pred_logit                 = "Logistic"
  ))

# Color palette
model_colors <- c(
  "Linear"                       = "#1b9e77",
  "Log(DHW)"                     = "blue",
  "Poisson"                      = "purple",
  "Poisson Log(DHW)"             = "red",
  "Negative Binomial (Log DHW)"  = "orange",
  "Beta Regression (FE)"         = "darkred",
  "ZOIB Ordbeta (RE)"            = "darkblue",
  "Logistic"         = "cyan4"
)

# Generate plot
ggplot(plot_data, aes(x = dhw)) +
  geom_point(aes(y = Percent_Bleached),
             alpha = 0.25, color = "gray50", size = 0.8) +
  geom_smooth(aes(y = predicted_bleaching, color = model),
              method = "loess", se = TRUE, linewidth = 1.1, span = 0.7) +
  scale_color_manual(values = model_colors) +
  scale_y_continuous(limits = c(0, 100)) +
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



## did not add negative binomial no log because percent bleaching goes over 100

# Add predictions from all relevant models

# model_data <- model_data %>%
#   mutate(
#     pred_linear                = predict(model_linear, newdata = ., type = "response"),
#     pred_log                   = predict(model_log, newdata = ., type = "response"),
#     pred_poisson               = predict(model_poisson, newdata = ., type = "response"),
#     pred_poisson_log           = predict(model_poisson_log, newdata = ., type = "response"),
#     pred_negative_binomial     = predict(model_negative_binomial, newdata = ., type = "response")
#   )
# 
# # Filter out rows with missing predictions (common for FE models)
# model_data <- model_data %>%
#   filter(
#     !is.na(pred_linear),
#     !is.na(pred_log),
#     !is.na(pred_poisson),
#     !is.na(pred_poisson_log),
#     !is.na(pred_negative_binomial)
#   )
# 
# # Reshape to long format for ggplot
# plot_data <- model_data %>%
#   dplyr::select(
#     dhw, Percent_Bleached,
#     pred_linear, pred_log,
#     pred_poisson, pred_poisson_log,
#     pred_negative_binomial
#   ) %>%
#   tidyr::pivot_longer(
#     cols = dplyr::starts_with("pred_"),
#     names_to = "model",
#     values_to = "predicted_bleaching"
#   ) %>%
#   dplyr::mutate(model = dplyr::recode(
#     model,
#     pred_linear = "Linear",
#     pred_log = "Log(DHW)",
#     pred_poisson = "Poisson",
#     pred_poisson_log = "Poisson Log(DHW)",
#     pred_negative_binomial = "Negative Binomial (Log DHW)"
#     ))
# 
# 
# # Define custom color palette for publication-style clarity
# model_colors <- c(
#   "Linear" = "#1b9e77",
#   "Log(DHW)" = "blue",
#   "Poisson" = "purple",
#   "Poisson Log(DHW)" = "red",
#   "Negative Binomial (Log DHW)" = "orange"
# )
# 
# # Generate plot
# ggplot(plot_data, aes(x = dhw)) +
#   geom_point(aes(y = Percent_Bleached),
#              alpha = 0.25, color = "gray50", size = 0.8) +
#   geom_smooth(aes(y = predicted_bleaching, color = model),
#               method = "loess", se = TRUE, size = 1.1, span = 0.7) +
#   scale_color_manual(values = model_colors) +
#   labs(
#     x = "Degree Heating Weeks (DHW)",
#     y = "Percent Bleached",
#     title = "Predicted Percent Bleaching vs DHW",
#     subtitle = "Observed bleaching (gray points) with fitted curves from multiple model specifications",
#     color = "Model Type"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(face = "bold", size = 15),
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold"),
#     panel.grid.minor = element_blank()
#   )


#### poisson or binomial ####
dispersion <- model$deviance / model$df.residual
model_poisson_log$deviance / model_poisson_log$df.residual

model_negative_binomial$theta

mean_bleach <- mean(All_Bleaching_Events_Data_AllDHW$Percent_Bleached, na.rm = TRUE)
var_bleach <- var(All_Bleaching_Events_Data_AllDHW$Percent_Bleached, na.rm = TRUE)
dispersion_raw <- var_bleach / mean_bleach
print(paste("Raw data dispersion ratio:", round(dispersion_raw, 2)))


# Compare AIC/BIC between Poisson and Negative Binomial
AIC(model_poisson_log)
AIC(model_negative_binomial)  # Your negbin model with log(dhw)

BIC(model_poisson_log)
BIC(model_negative_binomial)
# Much lower AIC/BIC for Negbin = strong evidence of overdispersion


etable(model_linear, model_log, model_poisson, model_poisson_log, model_negative_binomial, 
       dict = c("dhw" = "DHW", "log_dhw" = "log(DHW)"), 
       fitstat = c("aic", "bic", "r2"))


#### ecoregion curves ####
#### STEP 0: Define regional groups ####

data_per_ecoregion <- All_Bleaching_Events_Data_AllDHW %>%
  group_by(Ecoregion_Name) %>%
  summarise(n_obs = n(), n_sites = n_distinct(Site_ID), .groups = "drop") %>%
  arrange(desc(n_obs))

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



## with neg bin


p_faceted <- ggplot(model_data, aes(x = dhw)) +
  geom_point(aes(y = Percent_Bleached, color = "Observed"),
             alpha = 0.4, size = 0.8) +
  # Poisson log model with loess smoothing (constrained to >= 0)
  geom_smooth(aes(y = pmax(pred_poisson_log, 0), color = "Poisson Log"),
              method = "loess", size = 1, se = TRUE, span = 0.8) +
  # Negative Binomial with loess smoothing (constrained to >= 0)
  geom_smooth(aes(y = pmax(pred_negbin_log, 0), color = "Negative Binomial"),
              method = "loess", size = 1.2, se = TRUE, span = 0.8) +
  scale_color_manual(
    name = "Model Type",
    values = c("Observed" = "gray50", 
               "Poisson Log" = "red", 
               "Negative Binomial" = "orange")
  ) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +  # Constrain both axes to >= 0
  facet_wrap(~ Ecoregion_Name, scales = "free_y", ncol = 4) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Coral Bleaching by Ecoregion",
    subtitle = "Both models shown with loess smoothing and 95% CI"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
print(p_faceted)


# Identify top 8 ecoregions by number of observations
top_ecoregions <- model_data %>%
  count(Ecoregion_Name, sort = TRUE) %>%
  slice_head(n = 8) %>%
  pull(Ecoregion_Name)

# Filter data to top 8 ecoregions
model_data_top8 <- model_data %>%
  filter(Ecoregion_Name %in% top_ecoregions)

# Plot with both x and y constrained to >= 0
p_faceted_top8 <- ggplot(model_data_top8, aes(x = dhw)) +
  geom_point(aes(y = Percent_Bleached, color = "Observed"),
             alpha = 0.4, size = 0.8) +
  # Poisson log model with loess smoothing (constrained to >= 0)
  geom_smooth(aes(y = pmax(pred_poisson_log, 0), color = "Poisson Log"),
              method = "loess", size = 1, se = TRUE, span = 0.8) +
  # Negative Binomial with loess smoothing (constrained to >= 0)
  geom_smooth(aes(y = pmax(pred_negbin_log, 0), color = "Negative Binomial"),
              method = "loess", size = 1.2, se = TRUE, span = 0.8) +
  scale_color_manual(
    name = "Model Type",
    values = c("Observed" = "gray50", 
               "Poisson Log" = "red", 
               "Negative Binomial" = "orange")
  ) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +  # Constrain both axes to >= 0
  facet_wrap(~ Ecoregion_Name, scales = "free_y", ncol = 4) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Coral Bleaching by Ecoregion (Top 8 Regions)",
    subtitle = "Both models shown with loess smoothing and 95% CI"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p_faceted_top8)

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

# 
# 
# #### demeaned data and plots ####
# # ===== 1. Create Demeaned Dataset =====
# All_Bleaching_Events_Data_Demeaned <- All_Bleaching_Events_Data_AllDHW %>%
#   group_by(Site_ID, Date_Year) %>%
#   mutate(
#     mean_log_dhw = mean(log1p(dhw), na.rm = TRUE),
#     log_dhw_demeaned = log1p(dhw) - mean_log_dhw
#   ) %>%
#   ungroup()
# 
# # ===== 2. Fit Both Models =====
# original_model <- fepois(
#   Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year,
#   cluster = ~Ecoregion_Name,
#   data = All_Bleaching_Events_Data_AllDHW
# )
# 
# demeaned_model <- fepois(
#   Percent_Bleached ~ log_dhw_demeaned,
#   cluster = ~Ecoregion_Name,
#   data = All_Bleaching_Events_Data_Demeaned
# )
# 
# # ===== 3. Generate Predictions =====
# All_Bleaching_Events_Data_AllDHW$pred_original <- predict(
#   original_model, 
#   All_Bleaching_Events_Data_AllDHW
# )
# 
# All_Bleaching_Events_Data_Demeaned$pred_demeaned <- predict(
#   demeaned_model,
#   All_Bleaching_Events_Data_Demeaned
# )
# 
# # Prepare data for plotting
# plot_data_fitted <- bind_rows(
#   All_Bleaching_Events_Data_AllDHW %>%
#     select(dhw, Percent_Bleached, predicted_bleaching = pred_original) %>%
#     mutate(model = "Original (with FE)"),
#   All_Bleaching_Events_Data_Demeaned %>%
#     select(dhw, Percent_Bleached, predicted_bleaching = pred_demeaned) %>%
#     mutate(model = "Demeaned")
# )
# 
# # ===== 4. Create Plot =====
# model_colors <- c(
#   "Original (with FE)" = "#1b9e77",
#   "Demeaned" = "#d95f02"
# )
# 
# # Separate the data by model
# plot_data_original <- plot_data_fitted %>% 
#   filter(model == "Original (with FE)")
# 
# plot_data_demeaned <- plot_data_fitted %>% 
#   filter(model == "Demeaned")
# 
# ggplot(plot_data_fitted, aes(x = dhw)) +
#   geom_point(
#     aes(y = Percent_Bleached),
#     alpha = 0.25, 
#     color = "gray50", 
#     size = 0.8
#   ) +
#   # Original model - NO confidence interval
#   geom_smooth(
#     data = plot_data_original,
#     aes(y = predicted_bleaching, color = model),
#     method = "loess", 
#     se = TRUE,  # CI for original model
#     linewidth = 1.1, 
#     span = 0.7, 
#     alpha = 0.2
#   ) +
#   # Demeaned model - WITH confidence interval
#   geom_smooth(
#     data = plot_data_demeaned,
#     aes(y = predicted_bleaching, color = model),
#     method = "loess", 
#     se = TRUE,  # CI for demeaned model
#     linewidth = 1.1, 
#     span = 0.7,
#     alpha = 0.2
#   ) +
#   scale_color_manual(values = model_colors) +
#   scale_fill_manual(values = model_colors) +
#   labs(
#     x = "Degree Heating Weeks (DHW)",
#     y = "Percent Bleached",
#     title = "Predicted Percent Bleaching vs DHW: Original vs Demeaned Model",
#     subtitle = "Observed bleaching (gray points) with fitted curves from both model specifications",
#     color = "Model Type",
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(face = "bold", size = 15),
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold"),
#     panel.grid.minor = element_blank()
#   )
# 


#### realm plots ####
model_data_realm <- All_Bleaching_Events_Data_AllDHW %>%
  filter(
    !is.na(Percent_Bleached),
    !is.na(dhw)) %>%
  mutate(
    pred_poisson_log = predict(model_poisson_log, newdata = ., type = "response"),
    pred_negbin = predict(model_negative_binomial, newdata = ., type = "response"),
    pred_logit = predict(model_logit, newdata = ., type = "response") * 100  # MULTIPLY BY 100
  ) %>%
  filter(!is.na(pred_poisson_log))


ggplot(model_data_realm, aes(x = dhw)) +
  geom_point(aes(y = Percent_Bleached),
             alpha = 0.4, color = "gray50", size = 0.8) +
  geom_smooth(aes(y = pred_poisson_log),
              method = "loess", color = "red", size = 1, se = TRUE, span = 0.8) +
  geom_smooth(aes(y = pred_negbin),
              method = "loess", color = "blue", size = 1, se = TRUE, span = 0.8) +
  geom_smooth(aes(y = pred_logit),
              method = "loess", color = "green", size = 1, se = TRUE, span = 0.8) +
  facet_wrap(~ Realm_Name, scales = "free_y", ncol = 2) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Coral Bleaching by Realm (Poisson Log Model)",
    subtitle = "Gray points = observed bleaching | Red line = model predictions with 95% CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )


## binned
library(ggplot2)
library(dplyr)

# Define bins
custom_bins <- c(0, 4, 8, 12, 16, 20, 30)

# Create binned data with actual bin boundaries
# model_data_binned <- model_data_realm %>%
#   mutate(dhw_bin = cut(dhw, breaks = custom_bins, include.lowest = TRUE)) %>%
#   group_by(Realm_Name, dhw_bin) %>%
#   summarise(
#     mean_bleached = mean(Percent_Bleached, na.rm = TRUE),
#     n = n(),
#     .groups = "drop"
#   ) %>%
#   filter(!is.na(dhw_bin)) %>%
#   # Extract the actual bin boundaries from the factor levels
#   mutate(
#     dhw_lower = as.numeric(sub("\\((.+),.*", "\\1", dhw_bin)),
#     dhw_lower = ifelse(is.na(dhw_lower), as.numeric(sub("\\[(.+),.*", "\\1", dhw_bin)), dhw_lower),
#     dhw_upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", dhw_bin))
#   )


binned_data <- model_data_realm %>%
  mutate(dhw_bin = cut(dhw, breaks = custom_bins, include.lowest = TRUE)) %>%
  group_by(Realm_Name, dhw_bin) %>%
  summarise(
    mean_bleached = mean(Percent_Bleached, na.rm = TRUE),
    mean_pred_poisson = mean(pred_poisson_log, na.rm = TRUE),
    mean_pred_negbin = mean(pred_negbin, na.rm = TRUE),
    mean_pred_logit = mean(pred_logit, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(dhw_bin)) %>%
  mutate(
    bin_clean = gsub("\\[|\\]|\\(|\\)", "", dhw_bin),
    dhw_lower = as.numeric(sub(",.*", "", bin_clean)),
    dhw_upper = as.numeric(sub(".*,", "", bin_clean)),
    dhw_mid = (dhw_lower + dhw_upper) / 2
  )

# Create the plot
ggplot(model_data_realm, aes(x = dhw)) +
  # Raw data points
  geom_point(aes(y = Percent_Bleached),
             alpha = 0.3, color = "gray50", size = 0.8) +
  # Binned bars showing mean bleaching
  geom_rect(
    data = model_data_binned,
    aes(
      xmin = dhw_lower, 
      xmax = dhw_upper,
      ymin = 0,
      ymax = mean_bleached
    ),
    fill = "tan",
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  # Poisson predictions
  geom_smooth(aes(y = pred_poisson_log),
              method = "loess", color = "red", size = 1, se = TRUE, span = 0.8) +
  # Negative Binomial predictions
  geom_smooth(aes(y = pred_negbin),
              method = "loess", color = "blue", size = 1, se = TRUE, span = 0.8) +
  geom_smooth(aes(y = pred_logit),
              method = "loess", color = "green", size = 1, se = TRUE, span = 0.8) +
  facet_wrap(~ Realm_Name, scales = "free_y", ncol = 2) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Coral Bleaching by Realm",
    subtitle = "Gray points = raw data | Tan bars = binned mean | Red = Poisson | Blue = Negative Binomial | Green = Logit"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )




## confine y axis
# Create the plot
ggplot(model_data_realm, aes(x = dhw)) +
  # Raw data points
  geom_point(aes(y = Percent_Bleached),
             alpha = 0.3, color = "gray50", size = 0.8) +
  # Binned bars showing mean bleaching
  geom_rect(
    data = model_data_binned,
    aes(
      xmin = dhw_lower, 
      xmax = dhw_upper,
      ymin = 0,
      ymax = mean_bleached
    ),
    fill = "tan",
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  # Poisson predictions
  geom_smooth(aes(y = pred_poisson_log),
              method = "loess", color = "red", size = 1, se = TRUE, span = 0.8) +
  # Negative Binomial predictions
  geom_smooth(aes(y = pred_negbin),
              method = "loess", color = "blue", size = 1, se = TRUE, span = 0.8) +
  geom_smooth(aes(y = pred_logit),
              method = "loess", color = "green", size = 1, se = TRUE, span = 0.8) +
  facet_wrap(~ Realm_Name, scales = "free_y", ncol = 2) +
  coord_cartesian(ylim = c(0, 100)) +  # Add this line
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Coral Bleaching by Realm",
    subtitle = "Gray points = raw data | Tan bars = binned mean | Red = Poisson | Blue = Negative Binomial | Green = Logit"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )



#### map of realms ####

# Load required libraries
library(ggplot2)
library(maps)
library(dplyr)

# Get unique site locations with realm information
site_locations <- All_Bleaching_Events_Data_AllDHW %>%
  select(Site_ID, Latitude_Degrees, Longitude_Degrees, Realm_Name, 
         Country_Name, Ecoregion_Name) %>%
  distinct(Site_ID, .keep_all = TRUE)

# Get world map data
world_map <- map_data("world")

# Define color palette for realms
realm_colors <- c(
  "Central Indo-Pacific" = "#FF6B6B",
  "Tropical Atlantic" = "#FFEAA7",
  "Temperate Australasia" = "#96CEB4",
  "Western Indo-Pacific" = "#4ECDC4",
  "Eastern Indo-Pacific" = "#45B7D1",
  "Tropical Eastern Pacific" = "#DDA15E",
  "Temperate Northern Pacific" = "#A29BFE"
)

# Create the map
map_plot <- ggplot() +
  # Add world map
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70", linewidth = 0.3) +
  # Add site points
  geom_point(data = site_locations,
             aes(x = Longitude_Degrees, y = Latitude_Degrees, 
                 color = Realm_Name, fill = Realm_Name),
             size = 1, alpha = 0.8, shape = 21, stroke = 0.5) +
  # Customize colors
  scale_color_manual(values = realm_colors, name = "Marine Realm") +
  scale_fill_manual(values = realm_colors, name = "Marine Realm") +
  # Set coordinate system
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 60)) +
  # Customize theme
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "aliceblue", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Global Coral Bleaching Monitoring Sites",
    subtitle = paste0("Distribution of ", nrow(site_locations), 
                      " monitoring sites across marine realms"),
    x = "Longitude",
    y = "Latitude"
  ) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 4)))

# Display the map
print(map_plot)


# Print summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total unique sites:", nrow(site_locations), "\n\n")
cat("Sites by Marine Realm:\n")
realm_summary <- site_locations %>%
  group_by(Realm_Name) %>%
  summarise(n_sites = n()) %>%
  arrange(desc(n_sites))
print(realm_summary)

cat("\n\nSites by Country (Top 10):\n")
country_summary <- site_locations %>%
  group_by(Country_Name) %>%
  summarise(n_sites = n()) %>%
  arrange(desc(n_sites)) %>%
  head(10)
print(country_summary)


# Print summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total unique sites:", nrow(site_locations), "\n")
cat("Total observations:", nrow(All_Bleaching_Events_Data_AllDHW), "\n\n")

cat("Sites and Observations by Marine Realm:\n")
realm_summary <- All_Bleaching_Events_Data_AllDHW %>%
  group_by(Realm_Name) %>%
  summarise(
    n_observations = n(),
    n_sites = n_distinct(Site_ID)
  ) %>%
  arrange(desc(n_observations))
print(realm_summary)

cat("\n\nSites by Country (Top 10):\n")
country_summary <- site_locations %>%
  group_by(Country_Name) %>%
  summarise(n_sites = n()) %>%
  arrange(desc(n_sites)) %>%
  head(10)
print(country_summary)

# Create a bar plot showing observations by realm
obs_plot <- ggplot(realm_summary, aes(x = reorder(Realm_Name, n_observations), 
                                      y = n_observations, fill = Realm_Name)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = realm_colors) +
  coord_flip() +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "Number of Bleaching Observations by Marine Realm",
    x = "Marine Realm",
    y = "Number of Observations"
  ) +
  geom_text(aes(label = n_observations), hjust = -0.2, size = 3.5)

print(obs_plot)




#### density plot ####
# Load required libraries
library(ggplot2)
library(maps)
library(dplyr)
library(viridis)

# Get world map data
world_map <- map_data("world")

# Prepare the data - all observations with coordinates
obs_data <- All_Bleaching_Events_Data_AllDHW %>%
  select(Latitude_Degrees, Longitude_Degrees, Realm_Name) %>%
  filter(!is.na(Latitude_Degrees) & !is.na(Longitude_Degrees))

# Function to convert lat/lon to approximate 5km grid cells
# At the equator, 1 degree ≈ 111 km, so 5km ≈ 0.045 degrees
# We'll use 0.045 degrees as our grid size (approximately 5km)
grid_size <- 1

# Create 5km pixel grid
obs_data_grid <- obs_data %>%
  mutate(
    # Round to nearest grid cell
    lon_grid = round(Longitude_Degrees / grid_size) * grid_size,
    lat_grid = round(Latitude_Degrees / grid_size) * grid_size
  ) %>%
  group_by(lon_grid, lat_grid) %>%
  summarise(
    n_observations = n(),
    realms = paste(unique(Realm_Name), collapse = ", "),
    .groups = "drop"
  )

# Create 5km pixel density map
pixel_density_map <- ggplot() +
  # Add world map
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray70", linewidth = 0.2) +
  # Add 5km pixels as tiles
  geom_tile(data = obs_data_grid,
            aes(x = lon_grid, y = lat_grid, fill = n_observations),
            width = grid_size, height = grid_size, alpha = 0.85) +
  # Use viridis color scale with log transformation
  scale_fill_viridis(option = "inferno", 
                     name = "Observations\nper 10km pixel",
                     trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000),
                     labels = c("1", "10", "100", "1,000", "10,000")) +
  # Set coordinate system
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 60)) +
  # Customize theme
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "aliceblue", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Coral Bleaching Observation Density (10km Grid)",
    subtitle = paste0("Each pixel represents a ~10km × 10km area | Total: ", 
                      nrow(obs_data), " observations in ", 
                      nrow(obs_data_grid), " grid cells"),
    x = "Longitude",
    y = "Latitude"
  )

print(pixel_density_map)






#### creating own regions and condensing realms ####
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Condensed_Realm = case_when(
      Realm_Name %in% c("Central Indo-Pacific", "Temperate Australasia") ~ "Central Indo-Pacific",
      TRUE ~ Realm_Name
    )
  )

All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Condensed_Realm2 = case_when(
      Realm_Name %in% c("Central Indo-Pacific", "Temperate Australasia", "Eastern Indo-Pacific") ~ "Central Indo-Pacific",
      TRUE ~ Realm_Name
    )
  )

model_negative_binomial <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                    data = All_Bleaching_Events_Data_AllDHW)
model_region_condensed_CIP <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                       data = All_Bleaching_Events_Data_AllDHW)
model_region_condensed_CIP2 <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                       data = All_Bleaching_Events_Data_AllDHW)
summary(model_negative_binomial)
summary(model_region_condensed_CIP)
summary(model_region_condensed_CIP2)


#### Condensed Realm Plots (Version 1: CIP + Temperate Australasia) ####

# Prepare data with predictions for Condensed_Realm
model_data_realm_condensed <- All_Bleaching_Events_Data_AllDHW %>%
  filter(
    !is.na(Percent_Bleached),
    !is.na(dhw)) %>%
  mutate(pred_poisson_log = predict(model_poisson_log, newdata = ., type = "response")) %>%
  mutate(pred_negbin = predict(model_region_condensed_CIP, newdata = ., type = "response")) %>%
  filter(!is.na(pred_poisson_log))

# Define bins
custom_bins <- c(0, 4, 8, 12, 16, 20, 30)

# Create binned data
model_data_binned_condensed <- model_data_realm_condensed %>%
  mutate(dhw_bin = cut(dhw, breaks = custom_bins, include.lowest = TRUE)) %>%
  group_by(Condensed_Realm, dhw_bin) %>%
  summarise(
    mean_bleached = mean(Percent_Bleached, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(dhw_bin)) %>%
  mutate(
    dhw_lower = as.numeric(sub("\\((.+),.*", "\\1", dhw_bin)),
    dhw_lower = ifelse(is.na(dhw_lower), as.numeric(sub("\\[(.+),.*", "\\1", dhw_bin)), dhw_lower),
    dhw_upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", dhw_bin))
  )

# Create the plot with y-axis confined
ggplot(model_data_realm_condensed, aes(x = dhw)) +
  geom_point(aes(y = Percent_Bleached),
             alpha = 0.3, color = "gray50", size = 0.8) +
  geom_rect(
    data = model_data_binned_condensed,
    aes(
      xmin = dhw_lower, 
      xmax = dhw_upper,
      ymin = 0,
      ymax = mean_bleached
    ),
    fill = "tan",
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  geom_smooth(aes(y = pred_poisson_log),
              method = "loess", color = "red", size = 1, se = TRUE, span = 0.8) +
  geom_smooth(aes(y = pred_negbin),
              method = "loess", color = "blue", size = 1, se = TRUE, span = 0.8) +
  facet_wrap(~ Condensed_Realm, scales = "free_y", ncol = 2) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Coral Bleaching by Condensed Realm (v1)",
    subtitle = "Gray points = raw data | Tan bars = binned mean | Red = Poisson | Blue = Negative Binomial"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )


#### Condensed Realm Plots (Version 2: CIP + Temperate Australasia + Eastern Indo-Pacific) ####

# Prepare data with predictions for Condensed_Realm2
model_data_realm_condensed2 <- All_Bleaching_Events_Data_AllDHW %>%
  filter(
    !is.na(Percent_Bleached),
    !is.na(dhw)) %>%
  mutate(pred_poisson_log = predict(model_poisson_log, newdata = ., type = "response")) %>%
  mutate(pred_negbin = predict(model_region_condensed_CIP2, newdata = ., type = "response")) %>%
  filter(!is.na(pred_poisson_log))

# Create binned data
model_data_binned_condensed2 <- model_data_realm_condensed2 %>%
  mutate(dhw_bin = cut(dhw, breaks = custom_bins, include.lowest = TRUE)) %>%
  group_by(Condensed_Realm2, dhw_bin) %>%
  summarise(
    mean_bleached = mean(Percent_Bleached, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(dhw_bin)) %>%
  mutate(
    dhw_lower = as.numeric(sub("\\((.+),.*", "\\1", dhw_bin)),
    dhw_lower = ifelse(is.na(dhw_lower), as.numeric(sub("\\[(.+),.*", "\\1", dhw_bin)), dhw_lower),
    dhw_upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", dhw_bin))
  )

# Create the plot with y-axis confined
ggplot(model_data_realm_condensed2, aes(x = dhw)) +
  geom_point(aes(y = Percent_Bleached),
             alpha = 0.3, color = "gray50", size = 0.8) +
  geom_rect(
    data = model_data_binned_condensed2,
    aes(
      xmin = dhw_lower, 
      xmax = dhw_upper,
      ymin = 0,
      ymax = mean_bleached
    ),
    fill = "tan",
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  geom_smooth(aes(y = pred_poisson_log),
              method = "loess", color = "red", size = 1, se = TRUE, span = 0.8) +
  geom_smooth(aes(y = pred_negbin),
              method = "loess", color = "blue", size = 1, se = TRUE, span = 0.8) +
  facet_wrap(~ Condensed_Realm2, scales = "free_y", ncol = 2) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Coral Bleaching by Condensed Realm (v2)",
    subtitle = "Gray points = raw data | Tan bars = binned mean | Red = Poisson | Blue = Negative Binomial"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

