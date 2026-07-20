### author: Puja Pande
### date: 21 January 2026
### join all usable data from Mermaid and GCBD
### covers: (1) joining/cleaning the GCBD SQLite tables, (2) downloading and joining
###         MERMAID data, (3) merging MERMAID + GCBD into a single site-year panel,
###         and (4) extracting NOAA DHW (degree heating week) values, incl. lagged
###         DHW terms, for each panel observation

#### load all necessary libraries #####
library(dplyr)
library(DBI)
library(RSQLite)
library(tidyverse)

#### read in the data #####
# Connect to your SQLite database - this should work for everyone using Dropbox
con <- dbConnect(SQLite(), "path/to/Coracle/GCBD_Data/Global_Coral_Bleaching_Database_SQLite_11_24_21.db")

# Get all table names
table_names <- dbListTables(con)

# Read each table into a named list of data frames
tables_list <- lapply(table_names, function(tbl) dbReadTable(con, tbl))
names(tables_list) <- table_names

# Done: disconnect
dbDisconnect(con)

# reading each table into its own data frame
for (tbl in table_names) {
  assign(tbl, tables_list[[tbl]])
}


#### joining the data to create one appropriate dataframe ####
## join necessary tables to Site_Info_tbl
Site_Info_tbl_updated <- left_join(Site_Info_tbl, Ocean_Name_LUT, by = join_by(Ocean_Name == Ocean_ID)) %>% 
  select(-Ocean_Name) %>% rename(Ocean_Name = Ocean_Name.y)
Site_Info_tbl_updated <- left_join(Site_Info_tbl_updated, Realm_Name_LUT, by = join_by(Realm_Name == Realm_ID)) %>% 
  select(-Realm_Name) %>% rename(Realm_Name = Realm_Name.y)
Site_Info_tbl_updated <- left_join(Site_Info_tbl_updated, Ecoregion_Name_LUT, by = join_by(Ecoregion_Name == Ecoregion_ID)) %>% 
  select(-Ecoregion_Name) %>% rename(Ecoregion_Name = Ecoregion_Name.y)
Site_Info_tbl_updated <- left_join(Site_Info_tbl_updated, Country_Name_LUT, by = join_by(Country_Name == Country_ID)) %>% 
  select(-Country_Name) %>% rename(Country_Name = Country_Name.y)
Site_Info_tbl_updated <- left_join(Site_Info_tbl_updated, State_Island_Province_Name_LUT, 
                                   by = join_by(State_Island_Province_Name == State_Island_Province_ID)) %>% 
  select(-State_Island_Province_Name) %>% 
  rename(State_Island_Province_Name = State_Island_Province_Name.y)
Site_Info_tbl_updated <- left_join(Site_Info_tbl_updated, City_Town_Name_LUT, by = join_by(City_Town_Name == City_Town_ID)) %>% 
  select(-City_Town_Name) %>% 
  rename(City_Town_Name = City_Town_Name.y)
Site_Info_tbl_updated <- left_join(Site_Info_tbl_updated, City_Town_Name_LUT, by = join_by(City_Town_Name_2 == City_Town_ID)) %>% 
  select(-City_Town_Name_2) %>% 
  rename(City_Town_Name_2 = City_Town_Name.y)
Site_Info_tbl_updated <- left_join(Site_Info_tbl_updated, City_Town_Name_LUT, by = join_by(City_Town_Name_3 == City_Town_ID)) %>% 
  select(-City_Town_Name_3) %>% 
  rename(City_Town_Name_3 = City_Town_Name)
Site_Info_tbl_updated <- left_join(Site_Info_tbl_updated, City_Town_Name_LUT, by = join_by(City_Town_Name_4 == City_Town_ID)) %>% 
  select(-City_Town_Name_4) %>% 
  rename(City_Town_Name_4 = City_Town_Name)
Site_Info_tbl_updated <- left_join(Site_Info_tbl_updated, Exposure_LUT, by = join_by(Exposure == Exposure_ID)) %>% select(-Exposure) %>% 
  rename(Exposure = Exposure.y)
# removing all columns starting with 'TRIAL'
Site_Info_tbl_updated <- Site_Info_tbl_updated %>% select(-starts_with("TRIAL"))

## joining to Cover_tbl
Cover_tbl_updated <- left_join(Cover_tbl, Substrate_Type_LUT, by = join_by(Substrate_Type == Substrate_ID)) %>% 
  select(-Substrate_Type) %>% rename(Substrate_Type = Substrate_Name)
Cover_tbl_updated <- Cover_tbl_updated %>% select(-starts_with("TRIAL"))
Cover_tbl_updated <- Cover_tbl_updated %>%
  rename(
    S1_Cover = S1,
    S2_Cover = S2,
    S3_Cover = S3,
    S4_Cover = S4
  )

## joining to Bleaching_tbl
Bleaching_tbl_updated <- left_join(Bleaching_tbl, Bleaching_Prevalence_Score_LUT, 
                                   by = join_by(Bleaching_Prevalence_Score == Bleaching_Prevalence_ID)) %>% 
  select(-Bleaching_Prevalence_Score) %>% 
  rename(Bleaching_Prevalence_Score = Bleaching_Prevalence_Score.y)
Bleaching_tbl_updated <- left_join(Bleaching_tbl_updated, Bleaching_Level_LUT, by = join_by(Bleaching_Level == Bleaching_Level_ID)) %>% 
  select(-Bleaching_Level) %>% 
  rename(Bleaching_Level = Bleaching_Level.y)
Bleaching_tbl_updated <- left_join(Bleaching_tbl_updated, Severity_Code_LUT, by = join_by(Severity_Code == Severity_ID)) %>% 
  select(-Severity_Code) %>% 
  rename(Severity_Code = Severity_Code.y)
Bleaching_tbl_updated <- Bleaching_tbl_updated %>% select(-starts_with("TRIAL"))
Bleaching_tbl_updated <- Bleaching_tbl_updated %>%
  rename(
    S1_Bleaching = S1,
    S2_Bleaching = S2,
    S3_Bleaching = S3,
    S4_Bleaching = S4
  )

## joining to Sample_Event_tbl
Sample_Event_tbl_updated <- left_join(Sample_Event_tbl, Environmental_tbl, by = "Sample_ID") 

# some Sample_IDs have multiple Environmental measurements
nrow(Sample_Event_tbl_updated) 
n_distinct(Sample_Event_tbl_updated$Sample_ID) 

# checking which ones
env_counts <- Sample_Event_tbl_updated %>%
  group_by(Sample_ID) %>%
  summarise(count = n()) %>%
  filter(count > 1)
# extracting rows of repeated sample ids
env_repeats <- Sample_Event_tbl_updated %>%
  filter(Sample_ID %in% env_counts$Sample_ID) %>%        
  arrange(Sample_ID) %>%                                
  filter(!Sample_ID %in% c(21121, 10323635, 10330547))
# removing the repeats from the main dataframe
Sample_Event_tbl_updated <- Sample_Event_tbl_updated %>%
  arrange(Sample_ID) %>%
  group_by(Sample_ID) %>%
  filter(
    Sample_ID %in% c(10323635, 10330547) | row_number() == 1
  ) %>%
  ungroup()
Sample_Event_tbl_updated <- Sample_Event_tbl_updated %>% select(-starts_with("TRIAL"))

## joining Sample_Event_tbl to Site_Info_tbl
Site_Info_tbl_updated <- full_join(Site_Info_tbl_updated, Sample_Event_tbl_updated, 
                                   by = c("Site_ID"))

## joining bleaching and cover table to site info table updated
Site_Info_tbl_updated <- full_join(Site_Info_tbl_updated, Bleaching_tbl_updated, by = "Sample_ID") 
Site_Info_tbl_final <- full_join(Site_Info_tbl_updated, Cover_tbl_updated, by = "Sample_ID")

## use site_info_tbl_final for raw, unfiltered data
############## joining complete

#### cleaning the data ####
## removing irrelevant columns
Site_Info_tbl_final <- Site_Info_tbl_final %>% select(-c(Data_Source, Site_Name, Comments.x, City_Town_Name_2, 
                                                         City_Town_Name_3, City_Town_Name_4, Comments.y, 
                                                         Comments, Percent_Bleaching_Old_Method, Bleaching_ID, 
                                                         Environmental_ID, Reef_ID, City_Town_Name.x, Cover_ID, Sample_ID, 
                                                         Quadrat_No))
## redefining some columns
Site_Info_tbl_final <- Site_Info_tbl_final %>%
  rowwise() %>%
  mutate(mean_bleaching = mean(c_across(S1_Bleaching:S4_Bleaching), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Severity_Code = case_when(
      mean_bleaching > 0  & mean_bleaching <= 10 ~ "Mild (1-10% Bleached)",
      mean_bleaching > 10 & mean_bleaching <= 50 ~ "Moderate (11-50% Bleached)",
      mean_bleaching > 50                        ~ "Severe (>50% Bleached)",
      mean_bleaching <= 0                        ~ "No Bleaching",
      TRUE ~ Severity_Code   # keep existing if none of the above
    )
  )

Site_Info_tbl_final <- Site_Info_tbl_final %>%
  mutate(Severity_Code = case_when(
    Percent_Bleached > 0 & Percent_Bleached <= 10 ~ "Mild (1-10% Bleached)",
    Percent_Bleached > 10 & Percent_Bleached <= 50 ~ "Moderate (11-50% Bleached)",
    Percent_Bleached > 50 ~ "Severe (>50% Bleached)",
    Percent_Bleached <= 0 ~ "No Bleaching",
    TRUE ~ Severity_Code
  ))


Site_Info_tbl_final <- Site_Info_tbl_final %>%
  mutate(Severity_Code = case_when(
    Bleaching_Prevalence_Score == ">50% Reef Area Bleached" ~ "Severe (>50% Bleached)",
    Bleaching_Prevalence_Score == "<= 10% Reef Area Bleached" ~ "Mild (1-10% Bleached)",
    Bleaching_Prevalence_Score == "25-50% Reef Area Bleached" ~ "Moderate (11-50% Bleached)",
    Bleaching_Prevalence_Score == "10-25% Reef Area Bleached" ~ "Moderate (11-50% Bleached)",
    Bleaching_Prevalence_Score == "No Bleaching" ~ "No Bleaching",
    TRUE ~ Severity_Code
  ))
Site_Info_tbl_final$Percent_Bleached <- ifelse(is.na(Site_Info_tbl_final$Percent_Bleached) == TRUE, 
                                               Site_Info_tbl_final$mean_bleaching, Site_Info_tbl_final$Percent_Bleached)
Site_Info_tbl_final <- Site_Info_tbl_final %>% select(-c(Bleaching_Prevalence_Score, mean_bleaching))
table(Site_Info_tbl_final$Percent_Bleached, useNA = "ifany")
table(Site_Info_tbl_final$Severity_Code, useNA = "ifany")
Site_Info_tbl_final <- Site_Info_tbl_final %>% filter(Percent_Bleached != 'NaN')

Site_Info_tbl_final <- Site_Info_tbl_final %>% filter(Bleaching_Level == "Population")

summary(Site_Info_tbl_final)

write.csv(Site_Info_tbl_final, "path/to/Coracle/Datasets/GCBD_all_levels_01July.csv", row.names = FALSE)
#write.csv(Site_Info_tbl_final_pop, "path/to/Coracle/Datasets/GCBD_pop_only_01July.csv", row.names = FALSE)



### mermaid
####  Load packages and libraries ####
## If this is the first time using mermaidr, install the package through "remotes"
remotes::install_github("data-mermaid/mermaidr")

library(mermaidr) #package to download data from datamermaid.org
library(tidyverse) #package that makes it easier to work with data
library(plotly) #for interactive plotting
library(htmlwidgets) #for saving plots at html files
library(dplyr)
library(lubridate)

#### Get data from MERMAID for creating aggregate visualizations ####
projects <- mermaid_get_my_projects()
allMermaidSampEventsTBL <- mermaidr::mermaid_get_summary_sampleevents()

allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% filter(!is.na(colonies_bleached_percent_bleached_avg))

allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>%
  select(!starts_with(c("beltfish", "benthicpqt_", "management", "data_policy")))

allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% select(!ends_with(c("Sand", "Rubble", "Seagrass", "Macroalgae", "Turf algae", 
                                                                           "Cyanobacteria", "Bare substrate", 
                                                                           "Other invertebrates", "Crustose coralline algae")))

allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% select(!starts_with(c("benthiclit_percent_cover_life_", 
                                                                             "benthicpit_percent_cover_life_histories",
                                                                             "habitatcomplexity", 
                                                                             "colonies_bleached_percent_cover_life_histories")))

allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% select(-c(project_id, project, tags, project_notes, site_notes, 
                                                                 observers, contact_link, benthicpit_sample_unit_count, 
                                                                 `benthicpit_percent_cover_benthic_category_sd_Hard coral`, 
                                                                 `benthicpit_percent_cover_benthic_category_sd_Soft coral`, 
                                                                 quadrat_benthic_percent_percent_algae_avg_avg), 
                                                              benthiclit_sample_unit_count, 
                                                              `benthiclit_percent_cover_benthic_category_sd_Hard coral`, 
                                                              `benthiclit_percent_cover_benthic_category_sd_Soft coral`)
allMermaidSampEventsTBL$colonies_bleached_sample_unit_count <- NULL


allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>%
  mutate(
    date = ymd(sample_date),
    year  = year(sample_date),
    month = month(sample_date),
    day   = day(sample_date)
  )


allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% select(c(country, site_id, site, latitude, longitude, reef_type, 
                                                                reef_zone, reef_exposure, sample_date, 
                                                                colonies_bleached_percent_bleached_avg, date, year, month, day))
allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% select(-c(site_id, reef_type, reef_zone, reef_exposure, sample_date))


## mermaid ecoregions
# Load ecoregions shapefile
library(sf)

# Load ecoregions (in Mercator)
ecoregions <- st_read("path/to/Coracle/Ecoregion_shapefiles/ecoregion_dataPolygon.shp")

# Transform to WGS84 (EPSG:4326) to match Mermaid coordinates
ecoregions <- st_transform(ecoregions, crs = 4326)

# Verify transformation worked
cat("Ecoregions CRS after transform:", st_crs(ecoregions)$input, "\n")
cat("Bounding box:", paste(st_bbox(ecoregions), collapse = ", "), "\n")

# Make geometries valid (important!)
ecoregions <- st_make_valid(ecoregions)

# Convert Mermaid to sf
allMermaidSampEventsTBL_sf <- allMermaidSampEventsTBL %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Now do the spatial join
allMermaidSampEventsTBL_joined <- allMermaidSampEventsTBL_sf %>%
  st_join(ecoregions, join = st_within)

# Drop geometry and clean up
allMermaidSampEventsTBL <- allMermaidSampEventsTBL_joined %>%
  st_drop_geometry() %>%
  select(-c(Confirmed, Predicted, Doubtful, Absent, nid, id))  # Remove ecoregion metadata columns

# Verify
allMermaidSampEventsTBL %>%
  count(Ecoregion) %>%
  arrange(desc(n)) %>%
  head(10)


allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>%
  filter(!is.na(colonies_bleached_percent_bleached_avg))

allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>%
  mutate(Site_ID = as.numeric(as.factor(site)) + 19999)

#### joining mermaid data with GCBD data ####
colnames(allMermaidSampEventsTBL)
colnames(Site_Info_tbl_final)

allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% 
  mutate(Ecoregion = if_else(country == "Kenya" & is.na(Ecoregion), "Kenya and Tanzania coast", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Madagascar" & is.na(Ecoregion), "North and north-east Madagascar", 
                             Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Madagascar" & Ecoregion == "North Madagascar", 
                             "North and north-east Madagascar", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Tanzania" & is.na(Ecoregion), "Kenya and Tanzania coast", 
                             Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "United States" & is.na(Ecoregion), "Eastern Hawaii", 
                             Ecoregion)) 

allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% mutate(Ecoregion = if_else(country == "Mozambique" & is.na(Ecoregion), 
                                                                                  "North Mozambique coast", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Philippines" & is.na(Ecoregion), "South-east Philippines", 
                             Ecoregion)) 


allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% 
  mutate(Ecoregion = if_else(country == "Tanzania" & is.na(Ecoregion), "Kenya and Tanzania coast", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "United States" & is.na(Ecoregion), 
                             "Eastern Hawaii", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Solomon Islands" & is.na(Ecoregion),
                             "Solomon Islands and Bougainville", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Philippines" & is.na(Ecoregion) & latitude > 10,
                             "Sulu Sea, Philippines", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Philippines" & is.na(Ecoregion) & latitude < 10,
                             "South-east Philippines", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Mozambique" & is.na(Ecoregion) & latitude > -15,
                             "North Mozambique coast", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Mozambique" & is.na(Ecoregion) & latitude < -15,
                             "South Mozambique coast", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Mayotte" & is.na(Ecoregion),
                             "Mayotte and Comoro Islands", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Mauritius" & is.na(Ecoregion),
                             "Mascarene Islands", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Malaysia" & is.na(Ecoregion),
                             "Sulu Sea, Philippines", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Madagascar" & is.na(Ecoregion),
                             "North and north-east Madagascar", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Kenya" & is.na(Ecoregion), "Kenya and Tanzania coast", Ecoregion)) 


allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% 
  mutate(Ecoregion = if_else(country == "Belize" & is.na(Ecoregion), 
                             "Belize and west Caribbean", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Comoros" & is.na(Ecoregion),
                             "Mayotte and Comoro Islands", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "India" & is.na(Ecoregion) & longitude > 80,
                             "Andaman Islands", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "India" & is.na(Ecoregion) & longitude < 80,
                             "Lakshadweep islands", Ecoregion))
allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% 
  mutate(Ecoregion = if_else(country == "Indonesia" & is.na(Ecoregion) & longitude < 120,
                             "Lesser Sunda Islands and Savu Sea", Ecoregion)) %>%
  mutate(Ecoregion = if_else(country == "Indonesia" & is.na(Ecoregion) & longitude > 120,
                             "Celebes Sea, Indonesia", Ecoregion))  

allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% 
  mutate(Ecoregion = if_else(country == "Fiji" & is.na(Ecoregion),
                             "Fiji", Ecoregion)) 
allMermaidSampEventsTBL <- allMermaidSampEventsTBL %>% filter(year <= 2024)


write.csv(allMermaidSampEventsTBL, "path/to/Coracle/Datasets/Mermaid_01July.csv", row.names = FALSE)







#### Join Mermaid with GCBD data ####

# First, standardize column names for joining
allMermaidSampEventsTBL_clean <- allMermaidSampEventsTBL %>%
  rename(
    Latitude_Degrees = latitude,
    Longitude_Degrees = longitude,
    Date_Year = year,
    Date_Month = month,
    Date_Day = day,
    Percent_Bleached = colonies_bleached_percent_bleached_avg,
    Ecoregion_Name = Ecoregion
  ) %>%
  mutate(
    Site_ID = Site_ID,  # Use Mermaid's id as Site_ID
    Source = "Mermaid"
  )

Site_Info_tbl_final <- Site_Info_tbl_final %>%
  mutate(Source = "GCBD")

# Combine the two datasets
Combined_Bleaching_Data <- bind_rows(
  allMermaidSampEventsTBL_clean %>%
    select(Site_ID, Latitude_Degrees, Longitude_Degrees, Date_Year, Date_Month, Date_Day,
           Percent_Bleached, Ecoregion_Name, Country_Name = country, Source),
  
  Site_Info_tbl_final %>%
    select(Site_ID, Latitude_Degrees, Longitude_Degrees, Date_Year, Date_Month, Date_Day,
           Percent_Bleached, Ecoregion_Name, Country_Name, Source, Distance_to_Shore, Turbidity, Ocean_Name, Realm_Name, 
           Country_Name, Depth_m, Bleaching_Level)
)

# Check the combined data
cat("From Mermaid:", sum(Combined_Bleaching_Data$Source == "Mermaid"), "\n")
cat("From GCBD:", sum(Combined_Bleaching_Data$Source == "GCBD"), "\n")

# Remove duplicates if any (same site, date, and bleaching value)
Combined_Bleaching_Data <- Combined_Bleaching_Data %>%
  distinct(Latitude_Degrees, Longitude_Degrees, Date_Year, Date_Month, Date_Day, 
           Percent_Bleached, .keep_all = TRUE)

Combined_Data <- Combined_Bleaching_Data
write.csv(Combined_Data, "path/to/Coracle/GCBD_Data/Combined_GCBD_Mermaid_Data_pop_only_01July.csv", row.names = FALSE)





Combined_Data <- Combined_Data %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Kenya" & is.na(Ecoregion_Name), "Kenya and Tanzania coast", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Madagascar" & is.na(Ecoregion_Name), "North and north-east Madagascar", 
                                  Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Madagascar" & Ecoregion_Name == "North Madagascar", 
                                  "North and north-east Madagascar", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Tanzania" & is.na(Ecoregion_Name), "Kenya and Tanzania coast", 
                                  Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "United States" & is.na(Ecoregion_Name), "Eastern Hawaii", 
                                  Ecoregion_Name)) 

Combined_Data <- Combined_Data %>% mutate(Ecoregion_Name = if_else(Country_Name == "Mozambique" & is.na(Ecoregion_Name), 
                                                                   "North Mozambique coast", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Philippines" & is.na(Ecoregion_Name), "South-east Philippines", 
                                  Ecoregion_Name)) 


Combined_Data <- Combined_Data %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Tanzania" & is.na(Ecoregion_Name), "Kenya and Tanzania coast", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "United States" & is.na(Ecoregion_Name), 
                                  "Eastern Hawaii", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Solomon Islands" & is.na(Ecoregion_Name),
                                  "Solomon Islands and Bougainville", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Philippines" & is.na(Ecoregion_Name) & Latitude_Degrees > 10,
                                  "Sulu Sea, Philippines", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Philippines" & is.na(Ecoregion_Name) & Latitude_Degrees < 10,
                                  "South-east Philippines", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Mozambique" & is.na(Ecoregion_Name) & Latitude_Degrees > -15,
                                  "North Mozambique coast", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Mozambique" & is.na(Ecoregion_Name) & Latitude_Degrees < -15,
                                  "South Mozambique coast", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Mayotte" & is.na(Ecoregion_Name),
                                  "Mayotte and Comoro Islands", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Mauritius" & is.na(Ecoregion_Name),
                                  "Mascarene Islands", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Malaysia" & is.na(Ecoregion_Name),
                                  "Sulu Sea, Philippines", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Madagascar" & is.na(Ecoregion_Name),
                                  "North and north-east Madagascar", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Kenya" & is.na(Ecoregion_Name), "Kenya and Tanzania coast", Ecoregion_Name)) 


Combined_Data <- Combined_Data %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Belize" & is.na(Ecoregion_Name), 
                                  "Belize and west Caribbean", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Comoros" & is.na(Ecoregion_Name),
                                  "Mayotte and Comoro Islands", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "India" & is.na(Ecoregion_Name) & Longitude_Degrees > 80,
                                  "Andaman Islands", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "India" & is.na(Ecoregion_Name) & Longitude_Degrees < 80,
                                  "Lakshadweep islands", Ecoregion_Name))
Combined_Data <- Combined_Data %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Indonesia" & is.na(Ecoregion_Name) & Longitude_Degrees < 120,
                                  "Lesser Sunda Islands and Savu Sea", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Indonesia" & is.na(Ecoregion_Name) & Longitude_Degrees > 120,
                                  "Celebes Sea, Indonesia", Ecoregion_Name))  

Combined_Data <- Combined_Data %>% 
  mutate(Ecoregion_Name = if_else(Country_Name == "Fiji" & is.na(Ecoregion_Name),
                                  "Fiji", Ecoregion_Name)) 

Combined_Data <- Combined_Data %>% filter(Date_Year <= 2024)

#write.csv(Combined_Data, "Final_Combined_Data_14May2026.csv", row.names = FALSE)
#write.csv(Combined_Data, "Final_Combined_Data_all_levels_01July2026.csv", row.names = FALSE)


##### panel data - 02/03/2026
panel_data <- Combined_Bleaching_Data %>%
  filter(!is.na(Percent_Bleached))

# removing site_id with only one observation
panel_data <- panel_data %>%
  group_by(Site_ID) %>%
  filter(n() > 1) %>%
  ungroup()

## fixing ecoregions
#change madagascar to north-east madagascar
panel_data <- panel_data %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Kenya" & is.na(Ecoregion_Name), "Kenya and Tanzania coast", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Madagascar" & is.na(Ecoregion_Name), "North and north-east Madagascar",
                                  Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Madagascar" & Ecoregion_Name == "North Madagascar",
                                  "North and north-east Madagascar", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Tanzania" & is.na(Ecoregion_Name), "Kenya and Tanzania coast",
                                  Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "United States" & is.na(Ecoregion_Name), "Eastern Hawaii",
                                  Ecoregion_Name))

panel_data <- panel_data %>% mutate(Ecoregion_Name = if_else(Country_Name == "Mozambique" & is.na(Ecoregion_Name),
                                                             "North Mozambique coast", Ecoregion_Name)) %>%
  mutate(Ecoregion_Name = if_else(Country_Name == "Philippines" & is.na(Ecoregion_Name), "South-east Philippines",
                                  Ecoregion_Name))
panel_data <- panel_data %>%
  group_by(Site_ID, Date_Year, Date_Month) %>%
  slice_max(Percent_Bleached, n = 1, with_ties = FALSE) %>%
  ungroup()

write.csv(panel_data, "path/to/Coracle/Datasets/panel_data_final_01July.csv", row.names = FALSE)


#### extracting DHW values ####
#### load libraries ####
library(terra)
library(lubridate)
library(dplyr)
library(pbapply)


#### creating the rasterbrick ####
root_dir <- "path/to/Coracle/NOAA_DHW_Data/"
rasters <- list.files(
  root_dir,
  pattern = "ct5km_dhw-max_v3.1_.*\\.nc$",
  full.names = TRUE,
  recursive = TRUE
)
rasters <- sort(rasters)
rasterbrick_dhw <- rast(rasters, subds = "degree_heating_week")

# check
rasterbrick_dhw
nlyr(rasterbrick_dhw)
# Define start and end dates
start_date <- as.Date("1985-04-01")
end_date   <- as.Date("2025-12-01")
# Generate a monthly sequence
date_seq <- seq(start_date, end_date, by = "month")
# Create names like DHW4-1985, DHW5-1985, ... , DHW10-2023
names(rasterbrick_dhw) <- paste0("DHW", format(date_seq, "%m-%Y"))

#writeRaster(rasterbrick_dhw, filename = "rasterbrick_dhw.tif")

#### reading in the panel data to extract DHW values or closest non-NA DHW values ####
## doing this for the entire GBCD dataset and the all_bleaching_events dataset
# full_panel_data <- read.csv("Final_Relevant_Scripts/GBCD_full_data_cleaned.csv")
# full_panel_data <- read.csv("Final_Relevant_Scripts/All_Bleaching_Events_Data.csv")
# full_panel_data <- read.csv("path/to/Coracle/Datasets/panel_data_final.csv")

full_panel_data <- panel_data

## keeps sites that appear more than once
full_panel_data <- full_panel_data %>%
  group_by(Site_ID) %>%
  filter(n() > 1) %>%
  ungroup()


#### extracting DHW values for each row in the panel data ####
# Add year-month column
full_panel_data$ym <- format(
  ymd(paste(full_panel_data$Date_Year, full_panel_data$Date_Month, "01", sep = "-")),
  "%Y-%m"
)
library(terra)
library(lubridate)

# Check what's in the raster
cat("Number of layers:", nlyr(rasterbrick_dhw), "\n")

# Extract dates from layer names
layer_names <- names(rasterbrick_dhw)

# Parse the date format from names (assuming "DHW04-1985" = April 1985)
raster_dhw_ym <- sapply(layer_names, function(name) {
  # Extract month and year from "DHW04-1985" format
  parts <- strsplit(name, "-")[[1]]
  month_part <- gsub("DHW", "", parts[1])  # Get "04" from "DHW04"
  year_part <- parts[2]  # Get "1985"
  
  paste(year_part, month_part, sep = "-")  # Return "1985-04"
})

# Now redo the matching
idx <- match(full_panel_data$ym, raster_dhw_ym)

# Now extract DHW
full_panel_data$dhw <- pbsapply(1:nrow(full_panel_data), function(i) {
  layer_idx <- idx[i]
  if (is.na(layer_idx)) return(NA)
  
  point <- terra::vect(
    cbind(full_panel_data$Longitude_Degrees[i], full_panel_data$Latitude_Degrees[i]), 
    crs = "EPSG:4326"
  )
  
  val <- terra::extract(rasterbrick_dhw[[layer_idx]], point)[1, 2]
  
  if (is.na(val)) {
    buff <- terra::buffer(point, width = 10000)
    vals <- terra::extract(rasterbrick_dhw[[layer_idx]], buff)[, 2]
    non_na_vals <- vals[!is.na(vals)]
    if (length(non_na_vals) > 0) {
      val <- non_na_vals[1]
    }
  }
  
  return(val)
})

#### creating DHW lags from -3 - 5 months for full_panel_data ####

# Pre-extract full time series for all sites (faster than point extraction in loop)
coords <- cbind(full_panel_data$Longitude_Degrees, full_panel_data$Latitude_Degrees)
sites_vect <- terra::vect(coords, crs = "EPSG:4326")
all_vals <- terra::extract(rasterbrick_dhw, sites_vect)[, -1]  # drop ID column
# all_vals: nrow = nrow(full_panel_data), ncol = nlyr(rasterbrick_dhw)

lags <- -3:5

lagged_df <- pblapply(lags, function(L) {
  shifted_idx <- idx - L
  valid <- !is.na(shifted_idx) & shifted_idx > 0 & shifted_idx <= ncol(all_vals)
  vals <- rep(NA_real_, length(idx))
  vals[valid] <- all_vals[cbind(seq_along(idx)[valid], shifted_idx[valid])]
  vals
})

# Attach lagged columns and fill NAs with buffer
for (j in seq_along(lags)) {
  col_name <- paste0("dhw_lag", lags[j])
  full_panel_data[[col_name]] <- lagged_df[[j]]
  
  na_indices <- which(is.na(full_panel_data[[col_name]]))
  if (length(na_indices) > 0) {
    full_panel_data[[col_name]][na_indices] <- pbsapply(na_indices, function(i) {
      layer_idx <- idx[i] - lags[j]
      if (is.na(layer_idx) || layer_idx < 1 || layer_idx > nlyr(rasterbrick_dhw)) return(NA)
      
      point <- terra::vect(
        cbind(full_panel_data$Longitude_Degrees[i], full_panel_data$Latitude_Degrees[i]),
        crs = "EPSG:4326"
      )
      
      buff <- terra::buffer(point, width = 10000)
      vals <- terra::extract(rasterbrick_dhw[[layer_idx]], buff)[, 2]
      non_na_vals <- vals[!is.na(vals)]
      if (length(non_na_vals) > 0) return(non_na_vals[1])
      return(NA)
    })
  }
}

# Check results
for (j in seq_along(lags)) {
  col_name <- paste0("dhw_lag", lags[j])
  cat(col_name, "- Non-NA:", sum(!is.na(full_panel_data[[col_name]])), "\n")
}


## taking the maximum bleaching value for each site, year and month and keep all variables
final_panel_data <- full_panel_data %>%
  group_by(Site_ID, Date_Year, Date_Month) %>%
  slice_max(Percent_Bleached, n = 1, with_ties = FALSE) %>%
  ungroup()

final_panel_data <- final_panel_data %>%
  group_by(Site_ID) %>%
  filter(n() > 1) %>%
  ungroup()

write.csv(final_panel_data, "path/to/Coracle/Datasets/Final_Panel_Data_allDHW_01July.csv", row.names = FALSE)
