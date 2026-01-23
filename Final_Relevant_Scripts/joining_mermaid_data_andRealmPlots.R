## author: Puja Pande
## date: 10 December 2025
## description: doing all the tasks from Slack 

#### load libraries ####
library(tidyverse)


#### files ####
final_completed_dataset <- read_csv("Final_Relevant_Scripts/Panel_Data_AllDHW.csv")
final_completed_dataset <- final_completed_dataset %>% filter(!is.na(Percent_Bleached))
final_completed_dataset <- final_completed_dataset %>% filter(Bleaching_Level == "Population")
final_completed_dataset <- final_completed_dataset %>% select(Site_ID, Latitude_Degrees, Longitude_Degrees, Ocean_Name, Realm_Name, 
                                                              Ecoregion_Name, Country_Name, Date_Day, Date_Month, Date_Year, 
                                                              Depth_m, Percent_Bleached, dhw)
## joining mermaid data

#### data ####
mermaid_data <- readRDS("Final_Relevant_Scripts/MermaidBleachingData_WithDHW.RDS")
mermaid_data <- mermaid_data %>%
  mutate(Site_ID = as.numeric(as.factor(site)) + 19999)

#### merging data ####
final_completed_dataset <- final_completed_dataset %>% select(c(Site_ID, Date_Year, Date_Month, Latitude_Degrees, Longitude_Degrees, Ecoregion_Name, Country_Name, 
                                                                                  Realm_Name, Ocean_Name, Country_Name, Date_Day, Depth_m, Percent_Bleached, dhw))
colnames(mermaid_data)
mermaid_data <- mermaid_data %>% select(c(country, Site_ID, month, latitude, longitude, year, colonies_bleached_percent_bleached_avg, dhw, day))
mermaid_data <- mermaid_data %>% rename(Country_Name = country,
                                        Date_Year = year,
                                        Date_Month = month,
                                        Latitude_Degrees = latitude,
                                        Longitude_Degrees = longitude,
                                        Percent_Bleached = colonies_bleached_percent_bleached_avg, 
                                        Date_Day = day)

# 
# # Create a common key for matching
# test1 <- final_completed_dataset %>%
#   mutate(lat_long_key = paste(round(Latitude_Degrees, 6), 
#                               round(Longitude_Degrees, 6), 
#                               sep = "_"))
# 
# test2 <- mermaid_data %>%
#   mutate(lat_long_key = paste(round(Latitude_Degrees, 6), 
#                               round(Longitude_Degrees, 6), 
#                               sep = "_"))
# 
# # Find matches
# common_locations <- inner_join(
#   test1 %>% select(lat_long_key) %>% distinct(),
#   test2 %>% select(lat_long_key) %>% distinct(),
#   by = "lat_long_key"
# )


joined_merblea <- full_join(final_completed_dataset, mermaid_data) 


joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Fiji" & is.na(Ecoregion_Name), "Fiji", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Fiji" & is.na(Realm_Name), "Central Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Fiji" & is.na(Ocean_Name), "Pacific", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Indonesia" & is.na(Ecoregion_Name), "Java Sea", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Indonesia" & is.na(Realm_Name), "Central Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Indonesia" & is.na(Ocean_Name), "Pacific", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Kenya" & is.na(Ecoregion_Name), "Kenya and Tanzania coast", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Kenya" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Kenya" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Philippines" & is.na(Ecoregion_Name), "South-east Philippines", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Philippines" & is.na(Realm_Name), "Central Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Philippines" & is.na(Ocean_Name), "Pacific", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "United States" & is.na(Ecoregion_Name), "Eastern Hawaii", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "United States" & is.na(Realm_Name), "Eastern Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "United States" & is.na(Ocean_Name), "Pacific", Ocean_Name)) 


joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Madagascar" & is.na(Ecoregion_Name) & Latitude_Degrees > -15, "North Madagascar", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Madagascar" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Madagascar" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Madagascar" & is.na(Ecoregion_Name) & Latitude_Degrees < -15, "South Madagascar", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Madagascar" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Madagascar" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Mauritius" & is.na(Ecoregion_Name), "Mascarene Islands", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Mauritius" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Mauritius" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Mayotte" & is.na(Ecoregion_Name), "Mayotte and Comoros", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Mayotte" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Mayotte" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Mozambique" & is.na(Ecoregion_Name) & Latitude_Degrees > -15, "North Mozambique", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Mozambique" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Mozambique" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Mozambique" & is.na(Ecoregion_Name) & Latitude_Degrees < -15, "South Mozambique", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Mozambique" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Mozambique" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Tanzania" & is.na(Ecoregion_Name), "Kenya and Tanzania coast", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Tanzania" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Tanzania" & is.na(Ocean_Name), "Indian", Ocean_Name)) 


joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Australia" & is.na(Ecoregion_Name), "Moreton Bay, eastern Australia", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Australia" & is.na(Realm_Name), "Temperate Australasia", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Australia" & is.na(Ocean_Name), "Pacific", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Belize" & is.na(Ecoregion_Name), "Belize and west Caribbean", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Belize" & is.na(Realm_Name), "Tropical Atlantic", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Belize" & is.na(Ocean_Name), "Atlantic", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "India" & is.na(Ecoregion_Name) & Longitude_Degrees < 85, "Lakshadweep Islands", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "India" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "India" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "India" & is.na(Ecoregion_Name) & Longitude_Degrees > 85, "Andaman Islands", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "India" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "India" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Malaysia" & is.na(Ecoregion_Name) & Longitude_Degrees > 85, "Sulu Sea", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Malaysia" & is.na(Realm_Name), "Central Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Malaysia" & is.na(Ocean_Name), "Pacific", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Maldives" & is.na(Ecoregion_Name) & Longitude_Degrees < 85, "Maldive Islands", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Maldives" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Maldives" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Solomon Islands" & is.na(Ecoregion_Name), "Solomon Islands and Bougainville", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Solomon Islands" & is.na(Realm_Name), "Central Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Solomon Islands" & is.na(Ocean_Name), "Pacific", Ocean_Name)) 

joined_merblea <- joined_merblea %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Comoros" & is.na(Ecoregion_Name), "Mayotte and Comoros", Ecoregion_Name)) %>%
  mutate(Realm_Name = if_else(Country_Name == "Comoros" & is.na(Realm_Name), "Western Indo-Pacific", Realm_Name)) %>%
  mutate(Ocean_Name = if_else(Country_Name == "Comoros" & is.na(Ocean_Name), "Indian", Ocean_Name)) 

joined_merblea$Ecoregion_Name <- ifelse(joined_merblea$Ecoregion_Name == "North Mozambique coast", "North Mozambique",
                         joined_merblea$Ecoregion_Name)

sort(unique(joined_merblea$Ecoregion_Name))
unique(joined_merblea$Realm_Name)
unique(joined_merblea$Ocean_Name)

#### saving merged data ####
saveRDS(joined_merblea, "Merged_Mermaid_Bleaching_Data.RDS")
write.csv(joined_merblea, "Merged_Mermaid_Bleaching_Data.csv", row.names = FALSE)


## filtering for panel
joined_merblea <- joined_merblea %>% 
  group_by(Site_ID) %>% 
  filter(n() > 1) %>% 
  ungroup()

saveRDS(joined_merblea, "Merged_Mermaid_Panel_Bleaching_Data.RDS")
write.csv(joined_merblea, "Merged_Mermaid_Panel_Bleaching_Data.csv", row.names = FALSE)



#### running through impact model ####
library(fixest)
model_negative_binomial <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                    data = final_completed_dataset)
summary(model_negative_binomial)

model_negative_binomial_merblea <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                            data = joined_merblea)
summary(model_negative_binomial_merblea)





#### realm plots ####
model_data_realm <- joined_merblea %>%
  filter(
    !is.na(Percent_Bleached),
    !is.na(dhw)) %>%
  mutate(pred_poisson_log = predict(model_poisson_log, newdata = ., type = "response")) %>%
  mutate(pred_negbin = predict(model_negative_binomial, newdata = ., type = "response")) %>%
  filter(!is.na(pred_poisson_log))


ggplot(model_data_realm, aes(x = dhw)) +
  geom_point(aes(y = Percent_Bleached),
             alpha = 0.4, color = "gray50", size = 0.8) +
  geom_smooth(aes(y = pred_poisson_log),
              method = "loess", color = "red", size = 1, se = TRUE, span = 0.8) +
  geom_smooth(aes(y = pred_negbin),
              method = "loess", color = "blue", size = 1, se = TRUE, span = 0.8) +
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
model_data_binned <- model_data_realm %>%
  mutate(dhw_bin = cut(dhw, breaks = custom_bins, include.lowest = TRUE)) %>%
  group_by(Realm_Name, dhw_bin) %>%
  summarise(
    mean_bleached = mean(Percent_Bleached, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(dhw_bin)) %>%
  # Extract the actual bin boundaries from the factor levels
  mutate(
    dhw_lower = as.numeric(sub("\\((.+),.*", "\\1", dhw_bin)),
    dhw_lower = ifelse(is.na(dhw_lower), as.numeric(sub("\\[(.+),.*", "\\1", dhw_bin)), dhw_lower),
    dhw_upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", dhw_bin))
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
  facet_wrap(~ Realm_Name, scales = "free_y", ncol = 2) +
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Coral Bleaching by Realm",
    subtitle = "Gray points = raw data | Tan bars = binned mean | Red = Poisson | Blue = Negative Binomial"
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
  facet_wrap(~ Realm_Name, scales = "free_y", ncol = 2) +
  coord_cartesian(ylim = c(0, 100)) +  # Add this line
  labs(
    x = "Degree Heating Weeks (DHW)",
    y = "Percent Bleached",
    title = "Predicted Coral Bleaching by Realm",
    subtitle = "Gray points = raw data | Tan bars = binned mean | Red = Poisson | Blue = Negative Binomial"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )


##. condensed realms



#### creating own regions and condensing realms ####
All_Bleaching_Events_Data_AllDHW <- joined_merblea %>%
  mutate(
    Condensed_Realm = case_when(
      Realm_Name %in% c("Central Indo-Pacific", "Temperate Australasia") ~ "Central Indo-Pacific",
      TRUE ~ Realm_Name
    )
  )

All_Bleaching_Events_Data_AllDHW <- joined_merblea %>%
  mutate(
    Condensed_Realm = case_when(
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
