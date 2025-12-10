#### author: Puja Pande
#### date: 08/12/2025
#### description: merging of panel and mermaid data + eda

#### libraries ####
library(dplyr)




#### data ####
mermaid_data <- readRDS("MermaidBleachingData_WithDHW.RDS")
mermaid_data <- mermaid_data %>% 
  group_by(site) %>% 
  filter(n() > 1) %>% 
  ungroup()

mermaid_data <- mermaid_data %>%
  mutate(Site_ID = as.numeric(as.factor(site)) + 19999)

## panel data 
All_Bleaching_Events_Data_AllDHW <- read.csv("All_Bleaching_Events_Data_AllDHW.csv")
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>% 
  group_by(Site_ID) %>% 
  filter(n() > 1) %>% 
  ungroup()

#### merging data ####
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>% select(c(Site_ID, Date_Year, Date_Month, Latitude_Degrees, Longitude_Degrees, Ecoregion_Name, Country_Name, 
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

joined_merblea <- full_join(All_Bleaching_Events_Data_AllDHW, mermaid_data) 


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


unique(joined_merblea$Ecoregion_Name)
unique(joined_merblea$Realm_Name)
unique(joined_merblea$Ocean_Name)

#### saving merged data ####
saveRDS(joined_merblea, "Merged_Mermaid_Panel_Bleaching_Data.RDS")
write.csv(joined_merblea, "Merged_Mermaid_Panel_Bleaching_Data.csv", row.names = FALSE)



#### running through impact model ####
library(fixest)
model_negative_binomial <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                    data = All_Bleaching_Events_Data_AllDHW)
summary(model_negative_binomial)

model_negative_binomial_merblea <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
                                    data = joined_merblea)
summary(model_negative_binomial_merblea)
