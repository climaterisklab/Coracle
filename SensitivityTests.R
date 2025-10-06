### author: Puja Pande
### date: 06 October 2025
### description: the script runs a variety of sensitivity test to see how the results change with different parameters and assumptions

#### load libraries ####
library(dplyr)
library(fixest)


#### load datasets ####
all_bleaching_events <- read.csv("All_Bleaching_Events_Data.csv")
all_bleaching_events_depth_filtered <- read.csv("All_Bleaching_Events_Data_Depth_Filtered.csv")
All_Bleaching_Events_Data_AllDHW <- read.csv("All_Bleaching_Events_Data_AllDHW.csv")



#### different model types ####
model_linear <- feols(Percent_Bleached ~ SSTA_DHW | Site_ID + Date_Year, cluster = ~Ecoregion_Name, data = all_bleaching_events)
model_log <- feols(Percent_Bleached ~ log1p(SSTA_DHW) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, data = all_bleaching_events)
model_poisson <- fepois(Percent_Bleached ~ SSTA_DHW | Site_ID + Date_Year, cluster = ~Ecoregion_Name, data = all_bleaching_events)
model_poisson_log <- fepois(Percent_Bleached ~ log1p(SSTA_DHW) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                            data = all_bleaching_events)

summary(model_linear)
summary(model_log)
summary(model_poisson)
summary(model_poisson_log)

## initially used felm model and the coefficients matched with feols
# summary(lfe::felm(Percent_Bleached ~ log1p(SSTA_DHW) | Site_ID + Date_Year | 0 | Ecoregion_Name, data = all_bleaching_events))

# ## converting and comparing coefficients to each other
# (marginal_effect_logX <- coef(model_log)["log1p(SSTA_DHW)"] / (1 + mean(all_bleaching_events$SSTA_DHW, na.rm = TRUE)))
# (percent_change_poisson <- 100 * (exp(coef(model_poisson)) - 1))
# (percent_change_poisson_log <- 100 * coef(model_poisson_log)) ## best - continue with this model
# 
# -------
# # Choose a meaningful change in SSTA_DHW
# delta_X <- 1  # or whatever is meaningful in your context
# mean_X <- mean(all_bleaching_events$SSTA_DHW, na.rm = TRUE)
# 
# # Model 1: Linear
# effect_linear <- coef(model_linear)["SSTA_DHW"] * delta_X
# 
# # Model 2: Poisson (% change in E[Y])
# effect_poisson <- 100 * (exp(coef(model_poisson)["SSTA_DHW"] * delta_X) - 1)
# 
# # Model 3: Linear with log(X) - evaluate at mean
# effect_log <- coef(model_log)["log1p(SSTA_DHW)"] * delta_X / (1 + mean_X)
# 
# # Model 4: Poisson with log(X) - for 1% increase
# effect_poisson_log <- coef(model_poisson_log)["log1p(SSTA_DHW)"]  # This IS the elasticity
# # OR for a specific % increase:
# pct_increase <- 10  # 10% increase in X
# effect_poisson_log_10pct <- pct_increase * coef(model_poisson_log)["log1p(SSTA_DHW)"]
# #A 10% increase in SSTA_DHW is associated with a 5.3% increase in expected bleaching counts
# 


#### different interactions ####
all_bleaching_events_edited <- all_bleaching_events %>% group_by(Site_ID) %>% 
  mutate(site_temp_ave = mean(ClimSST, na.rm = TRUE)) %>% ungroup() %>% mutate(site_temp_ave_int = site_temp_ave*SSTA_DHW) %>% 
  ungroup()
all_bleaching_events_edited <- all_bleaching_events_edited %>% group_by(Site_ID) %>% 
  mutate(depth_int = Depth_m*SSTA_DHW) %>% ungroup() 
all_bleaching_events_edited <- all_bleaching_events_edited %>% group_by(Site_ID) %>% 
  mutate(lat_int = Latitude_Degrees*SSTA_DHW) %>% ungroup() 

model_site_temp_ave_int <- fepois(Percent_Bleached ~ SSTA_DHW + site_temp_ave_int | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                  data = all_bleaching_events_edited)
model_depth_int <- fepois(Percent_Bleached ~ SSTA_DHW + depth_int | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                  data = all_bleaching_events_edited)
model_latitude_int <- fepois(Percent_Bleached ~ SSTA_DHW + lat_int | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                  data = all_bleaching_events_edited)

summary(model_site_temp_ave_int)     ### not significant
summary(model_depth_int)     ### not significant
summary(model_latitude_int)     ### not significant


#### depth sensitivity test ####
model_all_depths <- fepois(Percent_Bleached ~ SSTA_DHW | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                          data = all_bleaching_events)
model_depth_10andless <- fepois(Percent_Bleached ~ SSTA_DHW | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                data = all_bleaching_events_depth_filtered)
summary(model_all_depths)
summary(model_depth_10andless)           # 0.002 difference between the two models, both significant


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




