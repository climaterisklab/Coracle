### author: Puja Pande
### date: 30 September 2025
### description: script to read in GBCD data, manipulate and filter it and perform some light EDA

#### load all necessary libraries #####
library(DBI)
library(RSQLite)
library(tidyverse)

#### read in the data #####
# Connect to your SQLite database - this should work for everyone using Dropbox
con <- dbConnect(SQLite(), "~/Library/CloudStorage/Dropbox/Coracle/GCBD Data/Global_Coral_Bleaching_Database_SQLite_11_24_21.db")

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
nrow(Sample_Event_tbl_updated) ##27025
n_distinct(Sample_Event_tbl_updated$Sample_ID) #27008 
# checking which ones
env_counts <- Sample_Event_tbl_updated %>%
  group_by(Sample_ID) %>%
  summarise(count = n()) %>%
  filter(count > 1)
# extracting rows of repeated sample ids
env_repeats <- Sample_Event_tbl_updated %>%
  filter(Sample_ID %in% env_counts$Sample_ID) %>%        #### 10323635, 10330547 do not have the same values for bleaching
  arrange(Sample_ID) %>%                                #### 21121 does not have any values
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
#write_csv(Site_Info_tbl_final, "GBCD_full_data_cleaned.csv")
#### site_info_tbl_final now contains all the cleaned data

#### making a dataset that only contains data for bleaching events ####
## changing substrate and cover from wide to long format - easier for analysis and to read
all_bleach_info <- Site_Info_tbl_final 
all_bleach_info <- all_bleach_info %>%
  rowwise() %>%
  mutate(total_cover = sum(c_across(S1_Cover:S4_Cover), na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = Substrate_Type,
    values_from = total_cover,
    values_fill = NA,             # fill missing with 0
    values_fn = sum              # <-- combine duplicates by summing
  )

all_bleach_info$`Hard Coral` <- ifelse(is.na(all_bleach_info$S1_Cover) == TRUE & is.na(all_bleach_info$S2_Cover) == TRUE & 
                                         is.na(all_bleach_info$S3_Cover) == TRUE & 
                                         is.na(all_bleach_info$S4_Cover) == TRUE, NA, all_bleach_info$`Hard Coral`)

all_bleach_info$`Nutrient Indicator Algae` <- ifelse(is.na(all_bleach_info$S1_Cover) == TRUE & is.na(all_bleach_info$S2_Cover) == TRUE & 
                                                       is.na(all_bleach_info$S3_Cover) == TRUE & 
                                                       is.na(all_bleach_info$S4_Cover) == TRUE, NA, all_bleach_info$`Nutrient Indicator Algae`)

all_bleach_info$`Fleshy Seaweed` <- ifelse(is.na(all_bleach_info$S1_Cover) == TRUE & is.na(all_bleach_info$S2_Cover) == TRUE & 
                                             is.na(all_bleach_info$S3_Cover) == TRUE & 
                                             is.na(all_bleach_info$S4_Cover) == TRUE, NA, all_bleach_info$`Fleshy Seaweed`)

all_bleach_info <- all_bleach_info %>% select(-c(`NA`, S1_Cover, S2_Cover, S3_Cover, S4_Cover))


## i want to condense hard coral, nutrient indicator algae and fleshy seaweed into one row if site, year, month and bleaching level are the same
substrate_cols <- c("Hard Coral", "Nutrient Indicator Algae", "Fleshy Seaweed")

all_bleach_info <- all_bleach_info %>%
  group_by(Site_ID, Date_Year, Date_Month, Bleaching_Level) %>% 
  summarise(
    across(all_of(substrate_cols), ~ {
      if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE)
    }),
    across(-all_of(substrate_cols), ~ first(.x)),   # keep all other cols
    .groups = "drop"
  )
## filtering to only population level bleaching data
all_bleach_info <- all_bleach_info %>%
  filter(Bleaching_Level == 'Population')

## keeping sites with more than one bleaching event
# Count how often bleaching occurs per site (e.g., Percent_Bleached > 0)
bleach_counts <- Site_Info_tbl_final %>%
  group_by(Site_ID, Longitude_Degrees, Latitude_Degrees) %>%
  summarise(
    bleach_events = sum(!is.na(Severity_Code), na.rm = TRUE),
    .groups = "drop"
  )
all_bleach_info <- all_bleach_info %>%
  filter(Site_ID %in% bleach_counts$Site_ID[bleach_counts$bleach_events > 1])

## there aren't any rows that are NaN but this is just to be sure
all_bleach_info <- all_bleach_info %>%
  filter(Percent_Bleached != 'NaN') 

## creating a percent_bleached_max column for sensitivity test 
all_bleach_info <- all_bleach_info %>%
  rowwise() %>%
  mutate(Percent_Bleached_Max = max(c_across(S1_Bleaching:S4_Bleaching), na.rm = TRUE)) %>%
  ungroup()

## removing unnecessary columns
all_bleach_info <- all_bleach_info %>% select(-c(S1_Bleaching, S2_Bleaching, S3_Bleaching, S4_Bleaching, Number__Bleached_Colonies, 
                                                 bleach_intensity, Percent_Hard_Coral, Percent_Macroalgae))

## depth filtered <10
all_bleach_info_depth_filtered <- all_bleach_info %>%
  filter(Depth_m <= 10)

#### datasets ready to be used in a panel regression model ####
write_csv(Site_Info_tbl_final, "GBCD_full_data_cleaned.csv")
write_csv(all_bleach_info, "All_Bleaching_Events_Data.csv")
write_csv(all_bleach_info_depth_filtered, "All_Bleaching_Events_Data_Depth_Filtered.csv")









