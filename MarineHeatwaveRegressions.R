### author: Puja Pande
### date: 06 October 2025
### GMST - marine heatwave regressions

#### load libraries ####
library(dplyr)
library(lubridate)
library(purrr)

#### load gmst and DHW data ####
gmst_data <- read.csv("~/Library/CloudStorage/Dropbox/Coracle/ERA5_GMT.csv")
gmst_data$X <- NULL

dhw_data <- read.csv("Panel_Data_AllDHW.csv")
dhw_data$X <- NULL

#### quick plot ####
plot(gmst_data$Year, gmst_data$GMT, type = "l", main = "Observed GMST")

#### running poisson regression model for every pixel ####
## creating monthly dataframes
monthly_data <- map(1:12, ~ dhw_data %>% filter(Date_Month == .x))
names(monthly_data) <- month.name

plot(monthly_data$January$dhw)

## computing the yearly mean DHW for each month
monthly_data <- map(monthly_data, ~ {
  # Compute the yearly mean DHW for that month
  yearly_means <- .x %>%
    group_by(Date_Year) %>%
    summarise(mean_dhw = mean(dhw, na.rm = TRUE))
  
  # Join the yearly mean back to the original monthly data
  .x %>%
    left_join(yearly_means, by = "Date_Year")
})

## creating a new dataset with this information 
monthly_mean_dhw <- map(monthly_data, ~ {
  .x %>%
    group_by(Date_Year) %>%
    summarise(mean_dhw = mean(dhw, na.rm = TRUE))
})

## joining the GMST data to each monthly dataframe by year 
monthly_mean_dhw <- map(monthly_mean_dhw, ~ {
  .x %>%
    inner_join(gmst_data, by = c("Date_Year" = "Year"))
})




lm_fit <- glm(mean_dhw ~ GMT, family = "poisson", data = merged_df)
summary(lm_fit)
preds <- predict(lm_fit, gmst[47:83, ], type = "response")
