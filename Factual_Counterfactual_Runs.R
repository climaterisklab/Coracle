#### author: Puja Pande
#### date: 27 October 2025
#### description: factual and counterfactual runs 

#### libraries ####
library(fixest)
library(ncdf4)
library(dplyr)
library(tidyr)

#### load data and model ####
final_model <- fepois(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                      data = All_Bleaching_Events_Data_AllDHW)

#### test factual ####
factual_preds <- predict(final_model, 
                          newdata = All_Bleaching_Events_Data_AllDHW %>%
                            mutate(dhw = dhw),
                          type = "response")

plot(factual_preds)
#### test counterfactual ####