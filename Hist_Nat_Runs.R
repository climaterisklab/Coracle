#### author: Puja Pande
#### date: 27 October 2025
#### description: historical and natural runs

#### libraries ####
library(fixest)
library(ncdf4)
library(dplyr)
library(ggplot2)
library(reshape2)


#### load data and model ####
All_Bleaching_Events_Data_AllDHW <- read.csv("All_Bleaching_Events_Data_AllDHW.csv")
final_model <- fepois(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                      data = All_Bleaching_Events_Data_AllDHW)

nc_counts_counterfactuals <- nc_open("~/Library/CloudStorage/Dropbox/Coracle/Counterfactuals_Updated/counts_counterfactuals.nc")
nc_counts_historical <- nc_open("~/Library/CloudStorage/Dropbox/Coracle/Counterfactuals_Updated/counts_historical.nc")
#nc_lambda_counterfactuals <- nc_open("~/Library/CloudStorage/Dropbox/Coracle/Counterfactuals_Updated/lambda_counterfactuals.nc")
#nc_lambda_historical <- nc_open("~/Library/CloudStorage/Dropbox/Coracle/Counterfactuals_Updated/lambda_historical.nc")

#### eda ####
## historical
# check dimension values
months <- ncvar_get(nc_counts_historical, "month")
years <- ncvar_get(nc_counts_historical, "time")
locs <- ncvar_get(nc_counts_historical, "location_id")

dhw_hist <- ncvar_get(nc_counts_historical, "DHW")
dim(dhw_hist)

months <- 1:12
years <- seq(1981, 2020)  # since time size = 40
locs <- 1:13776

# Convert to dataframe
hist_df <- melt(dhw_hist, varnames = c("month", "location_id", "time"), value.name = "dhw")

# Replace indices with actual numeric coordinates
hist_df$month <- months[hist_df$month]
hist_df$time <- years[hist_df$time]
hist_df$location_id <- locs[hist_df$location_id]

# Add metadata
hist_df$Scenario <- "Historical"



## counterfactual
dhw_cf <- ncvar_get(nc_counts_counterfactuals, "DHW")
dim(dhw_cf)

# reshape
cf_df <- melt(dhw_cf, varnames = c("month", "location_id", "time", "firm"), value.name = "dhw")

# replace indices with coordinates
cf_df$month <- months[cf_df$month]
cf_df$time <- years[cf_df$time]
cf_df$location_id <- locs[cf_df$location_id]

# label scenarios
cf_df$Scenario <- paste0("CF_", cf_df$firm)

if (!"firm" %in% names(hist_df)) {
  hist_df$firm <- NA
}

cf_df <- cf_df[, names(hist_df)]

all_dhw <- rbind(hist_df, cf_df)



#### making predictions ####
# Merge with bleaching model predictors
all_dhw <- merge(all_dhw, All_Bleaching_Events_Data_AllDHW[, c("Site_ID", "Date_Year", "Percent_Bleached")],
                 by.x = c("location_id", "time"), by.y = c("Site_ID", "Date_Year"), all.x = TRUE)
all_dhw$dhw <- log1p(all_dhw$dhw)
all_dhw$Site_ID <- all_dhw$location_id
all_dhw$Date_Year <- all_dhw$time

all_dhw$pred_bleaching <- predict(final_model, newdata = all_dhw, type = "response")

all_dhw_no_nas <- all_dhw %>% filter(!is.na(pred_bleaching))


#### plotting ####
# Aggregate by year and scenario
annual_bleaching <- all_dhw %>%
  group_by(time, Scenario) %>%
  summarise(
    mean_pred = mean(pred_bleaching, na.rm = TRUE),
    median_pred = median(pred_bleaching, na.rm = TRUE)
  )

annual_bleaching <- annual_bleaching %>%
  filter(time >= 2000)

ggplot(annual_bleaching, aes(x = time, y = mean_pred, color = Scenario)) +
  geom_line(size = 1) +
  labs(title = "Predicted Bleaching: Historical vs Counterfactuals",
       x = "Year", y = "Mean Predicted Bleaching Probability") +
  theme_minimal() +
  scale_color_brewer(palette = "Set2")


# Filter for one site
site_id <- 12  # replace with desired location_id
site_data <- all_dhw %>% filter(location_id == site_id) %>% filter(time >= 2000)

# Plot predicted bleaching over time by scenario
ggplot(site_data, aes(x = time, y = pred_bleaching, color = Scenario)) +
  geom_line(size = 1.2) +
  labs(
    title = paste("Predicted Coral Bleaching Over Time - Site", site_id),
    x = "Year",
    y = "Predicted Bleaching (%)"
  ) +
  theme_minimal(base_size = 14)

