### author: Puja Pande
### date: 30 September 2024
### description: script to read in GBCD data, manipulate and filter it and perform some light EDA

#### load all necessary libraries #####
library(DBI)
library(RSQLite)
library(tidyverse)
##library(RColorBrewer)  # for better palettes - might not need this
library(sf)
library(mapview)

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







