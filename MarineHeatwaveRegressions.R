### author: Puja Pande
### date: 06 October 2025
### GMST - marine heatwave regressions

#### load libraries ####
library(dplyr)
library(lubridate)
library(purrr)

#### load gmst and DHW data ####
gmst_data <- read.csv("~/Library/CloudStorage/Dropbox/Coracle/ERA5_GMT.csv") %>% select(-X)
dhw_data <- read.csv("Panel_Data_AllDHW.csv") %>% select(-X)

#### quick plot ####
plot(gmst_data$Year, gmst_data$GMT, type = "l", main = "Observed GMST")

#### running poisson regression model for every pixel ####
## creating monthly dataframes
monthly_data <- map(1:12, ~ dhw_data %>% filter(Date_Month == .x))
names(monthly_data) <- month.name

plot(monthly_data$January$dhw)

#### Compute yearly mean DHW per pixel for each month ####
monthly_data <- map(1:12, ~ {
  dhw_data %>%
    filter(Date_Month == .x) %>%
    group_by(Date_Year, Latitude_Degrees, Longitude_Degrees) %>%
    summarise(mean_dhw = mean(dhw, na.rm = TRUE), .groups = "drop")
})
names(monthly_data) <- month.name

#### Join GMST data to each monthly dataset ####
monthly_data <- map(monthly_data, ~ 
                      inner_join(.x, gmst_data, by = c("Date_Year" = "Year"))
)

# 
# ## computing the yearly mean DHW for each month
# monthly_data <- map(monthly_data, ~ {
#   # Compute the yearly mean DHW for that month
#   yearly_means <- .x %>%
#     group_by(Date_Year, Latitude_Degrees, Longitude_Degrees) %>%
#     summarise(mean_dhw = mean(dhw, na.rm = TRUE, .groups = "drop"))
#   
#   # Join the yearly mean back to the original monthly data
#   .x %>%
#     left_join(yearly_means, by = c("Date_Year", "Latitude_Degrees", "Longitude_Degrees"))
# })
# 
# ## creating a new dataset with this information 
# monthly_mean_dhw <- map(monthly_data, ~ {
#   .x %>%
#     group_by(Date_Year) %>%
#     summarise(mean_dhw = mean(dhw, na.rm = TRUE))
# })
# 
# ## joining the GMST data to each monthly dataframe by year 
# monthly_mean_dhw <- map(monthly_mean_dhw, ~ {
#   .x %>%
#     inner_join(gmst_data, by = c("Date_Year" = "Year"))
# })


#### running poisson regression models for each month #### 
monthly_poisson_models <- map(monthly_data, ~ {
  glm(mean_dhw ~ GMT, family = "poisson", data = .x)
})
model_summaries <- map(monthly_poisson_models, summary)
model_summaries









preds <- predict(lm_fit, gmst[47:83, ], type = "response")
