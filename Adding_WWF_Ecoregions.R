#### author: Puja Pande
#### date: 04 November 2025
#### description: Adding WWF Ecoregions to bleaching events data

#### libraries ####
library(sf)
library(dplyr)

#### loading bleaching data ####
all_bleaching_events_depth_filtered <- read.csv("All_Bleaching_Events_Data_Depth_Filtered.csv")
All_Bleaching_Events_Data_AllDHW <- read.csv("All_Bleaching_Events_Data_AllDHW.csv")

#### loading the WWF ecoregions shapefile ####
wwf_ecoregions <- read_sf("~/Library/CloudStorage/Dropbox/Coracle/Global_200_Marine/Global_200_Marine.shp")

#### combining ####
bleaching_sf <- st_as_sf(
  all_bleaching_events_depth_filtered,
  coords = c("Longitude_Degrees", "Latitude_Degrees"),
  crs = 4326  # assuming coordinates are in WGS84
)

wwf_ecoregions <- st_transform(wwf_ecoregions, crs = st_crs(bleaching_sf))
bleaching_with_ecoregion <- st_join(bleaching_sf, wwf_ecoregions, join = st_intersects)
bleaching_with_ecoregion_df <- st_drop_geometry(bleaching_with_ecoregion)
bleaching_with_ecoregion_df <- bleaching_with_ecoregion_df %>% dplyr::select(-c(AREA, PERIMETER, G200_MAR_, G200_MAR_I, ECOREGION, 
                                                                         ECO_CODE, FOCAL, G200_NUMBE, ERBC, G200_NUM, MHT_NUM, 
                                                                         OBJECTID, Shape_Leng, Shape_Area))

#write.csv(bleaching_with_ecoregion_df, "All_Bleaching_Events_With_Ecoregion.csv", row.names = FALSE)
