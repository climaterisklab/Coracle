### author: Puja Pande
### date: 12 May 2026 - updated and new
### description: the script runs a variety of sensitivity test to see how the results change with different parameters and assumptions
### covers: model type comparisons, MPA interactions, cluster/lag/polynomial (quadratic,
###         cubic) sensitivity, Conley (200km spherical) SE checks, seasonality FE,
###         and DHW specification (DHW vs DHW_adj) sensitivity, with summary tables

#### load libraries ####
library(dplyr)
library(fixest)
library(ggplot2)
library(tidyr)
library(patchwork)
library(plotly)
library(sf)
library(glmmTMB)
library(modelsummary)
library(flextable)
library(officer)


#### load datasets ####
All_Bleaching_Events_Data_AllDHW <- read.csv("path/to/Coracle/Datasets/Final_Panel_Data_allDHW_01July.csv")
All_Bleaching_Events_Data_AllDHW$Proportion_Bleached <- All_Bleaching_Events_Data_AllDHW$Percent_Bleached / 100
All_Bleaching_Events_Data_AllDHW$abs_lat             <- abs(All_Bleaching_Events_Data_AllDHW$Latitude_Degrees)
All_Bleaching_Events_Data_AllDHW$mass_bleaching      <- ifelse(All_Bleaching_Events_Data_AllDHW$Percent_Bleached >= 30, 1, 0)
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Ecoregion_Month = interaction(Ecoregion_Name, Date_Month, drop = TRUE)
  )

All_Bleaching_Events_Data_AllDHW_Max <- read.csv("path/to/data/Final_Combined_Data_all_levels_01July2026_extravars.csv")
All_Bleaching_Events_Data_AllDHW_Max <- left_join(All_Bleaching_Events_Data_AllDHW, All_Bleaching_Events_Data_AllDHW_Max[, c(1:8, 17:34)], by = c("Site_ID", "Latitude_Degrees", "Longitude_Degrees", "Date_Year", "Date_Month", "Date_Day", "Ecoregion_Name", "Percent_Bleached"))


#### different model types ####
model_linear <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  cluster = ~Ecoregion_Name, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_turbidity <- feols(
  Percent_Bleached ~ dhw + dhw:Turbidity | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)

model_linear_depth <- feols(
  Percent_Bleached ~ dhw + dhw:Depth_m | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)

model_linear_dhw_std <- feols(
  Percent_Bleached ~ dhw + dhw:DHW_1985.2005_std | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW_Max
)

model_linear_hotspot_std <- feols(
  Percent_Bleached ~ dhw + dhw:hotspot_1985.2005_std | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW_Max
)

model_linear_hotspot_warming <- feols(
  Percent_Bleached ~ dhw + dhw:hotspot_warming | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW_Max
)

model_linear_lat_depth <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat + dhw:Depth_m | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW_Max
)

model_linear_lat_hotspot <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat + dhw:hotspot_1985.2005_std | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW_Max
)

#### Adding MPA interactions to bleaching events data ####
# Read all files
poly_files <- c(
  "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_0/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp",
  "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_1/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp",
  "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_2/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp"
)

point_files <- c(
  "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_0/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp",
  "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_1/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp",
  "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_2/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp"
)

# Merge polygons
poly_list <- lapply(poly_files, st_read)
merged_poly <- do.call(rbind, poly_list)

# Merge points
point_list <- lapply(point_files, st_read)
merged_points <- do.call(rbind, point_list)

# Get all unique column names
all_cols <- union(names(merged_poly), names(merged_points))

# Add missing columns to each dataset
for (col in all_cols) {
  if (!col %in% names(merged_poly)) {
    merged_poly[[col]] <- NA
  }
  if (!col %in% names(merged_points)) {
    merged_points[[col]] <- NA
  }
}

# Reorder columns to match
merged_poly <- merged_poly[, all_cols]
merged_points <- merged_points[, all_cols]

# Now combine
all_marine <- rbind(merged_poly, merged_points)

# Convert your data to sf object (points)
merged_sf <- st_as_sf(All_Bleaching_Events_Data_AllDHW, 
                      coords = c("Longitude_Degrees", "Latitude_Degrees"),  # Adjust column names
                      crs = 4326)  # WGS84

# Fix invalid geometries in MPA data
all_marine <- st_make_valid(all_marine)

# Spatial join
sf_use_s2(FALSE)
bleaching_with_mpa <- st_join(merged_sf, all_marine["SITE_ID"], left = TRUE)
sf_use_s2(TRUE)

# Add binary MPA indicator
bleaching_with_mpa$in_mpa <- !is.na(bleaching_with_mpa$SITE_ID)

# Drop geometry and deduplicate
# st_join can create duplicates if a point falls in multiple MPAs
All_Bleaching_Events_Data_AllDHW_MPA <- bleaching_with_mpa %>%
  st_drop_geometry() %>%
  group_by(across(-c(SITE_ID, in_mpa))) %>%
  summarise(in_mpa = any(in_mpa), .groups = "drop")

#### Run MPA model ####
model_linear_mpa <- feols(
  Percent_Bleached ~ dhw + dhw:in_mpa | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW_MPA
)


#### cluster sensitivity test ####
model_linear_nocluster <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_ecoregion <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  cluster = ~Ecoregion_Name, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_site <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  cluster = ~Site_ID, 
  data = All_Bleaching_Events_Data_AllDHW)


#### lag sensitivity test ####
model_lagminus3 <- feols(Percent_Bleached ~ dhw_lag.3 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                         data = All_Bleaching_Events_Data_AllDHW)
model_lagminus2 <- feols(Percent_Bleached ~ dhw_lag.2 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                         data = All_Bleaching_Events_Data_AllDHW)
model_lagminus1 <- feols(Percent_Bleached ~ dhw_lag.1 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                         data = All_Bleaching_Events_Data_AllDHW)
model_lag0 <- feols(Percent_Bleached ~ dhw_lag0 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)
model_lag1 <- feols(Percent_Bleached ~ dhw_lag1 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)
model_lag2 <- feols(Percent_Bleached ~ dhw_lag2 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)
model_lag3 <- feols(Percent_Bleached ~ dhw_lag3 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)
model_lag4 <- feols(Percent_Bleached ~ dhw_lag4 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)
model_lag5 <- feols(Percent_Bleached ~ dhw_lag5 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)



#### Quadratic DHW models ####
model_quadratic <- feols(
  Percent_Bleached ~ dhw + I(dhw^2) | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

model_quadratic_lat <- feols(
  Percent_Bleached ~ dhw + I(dhw^2) + dhw:abs_lat + I(dhw^2):abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

#### Cubic DHW models ####
model_cubic <- feols(
  Percent_Bleached ~ dhw + I(dhw^2) + I(dhw^3) | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

model_cubic_lat <- feols(
  Percent_Bleached ~ dhw + I(dhw^2) + I(dhw^3) +
    dhw:abs_lat + I(dhw^2):abs_lat + I(dhw^3):abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)






################## conley errors ################
#### Compute Conley (200km, spherical) vcov for each model ####
conley_vcov <- function(model) {
  vcov_conley(model,
              lat      = "Latitude_Degrees",
              lon      = "Longitude_Degrees",
              cutoff   = 200,
              distance = "spherical")
}

#### Different model types ####
model_linear_vcov <- conley_vcov(model_linear)
model_linear_lat_vcov <- conley_vcov(model_linear_lat)
model_linear_turbidity_vcov <- conley_vcov(model_linear_turbidity)
model_linear_depth_vcov <- conley_vcov(model_linear_depth)
model_linear_dhw_std_vcov <- conley_vcov(model_linear_dhw_std)
model_linear_hotspot_std_vcov <- conley_vcov(model_linear_hotspot_std)
model_linear_hotspot_warming_vcov <- conley_vcov(model_linear_hotspot_warming)
model_linear_lat_depth_vcov <- conley_vcov(model_linear_lat_depth)
model_linear_lat_hotspot_vcov <- conley_vcov(model_linear_lat_hotspot)

All_Bleaching_Events_Data_AllDHW_MPA <- All_Bleaching_Events_Data_AllDHW_MPA %>%
  left_join(
    All_Bleaching_Events_Data_AllDHW %>% 
      select(Site_ID, Date_Year, Latitude_Degrees, Longitude_Degrees),
    by = c("Site_ID", "Date_Year")
  )
model_linear_mpa_vcov <- conley_vcov(model_linear_mpa)


#### Lag models ####
model_lagminus3_vcov <- conley_vcov(model_lagminus3)
model_lagminus2_vcov <- conley_vcov(model_lagminus2)
model_lagminus1_vcov <- conley_vcov(model_lagminus1)
model_lag0_vcov <- conley_vcov(model_lag0)
model_lag1_vcov <- conley_vcov(model_lag1)
model_lag2_vcov <- conley_vcov(model_lag2)
model_lag3_vcov <- conley_vcov(model_lag3)
model_lag4_vcov <- conley_vcov(model_lag4)
model_lag5_vcov <- conley_vcov(model_lag5)


#### Quadratic models ####
model_quadratic_vcov <- conley_vcov(model_quadratic)
model_quadratic_lat_vcov <- conley_vcov(model_quadratic_lat)

#### Cubic models ####
model_cubic_vcov <- conley_vcov(model_cubic)
model_cubic_lat_vcov <- conley_vcov(model_cubic_lat)



#### Custom styling function ####
nature_table_style <- function(ft) {
  ft %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 7, part = "all") %>%
    align(align = "center", part = "header") %>%
    align(j = 1, align = "left", part = "body") %>%
    align(j = 2:ncol_keys(ft), align = "center", part = "body") %>%
    bold(part = "header") %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 1.5), part = "all") %>%
    hline_bottom(border = fp_border(width = 1.5), part = "all") %>%
    hline(i = 1, border = fp_border(width = 1), part = "header") %>%
    autofit() %>%
    width(width = 1.2, unit = "in")
}

#### Interaction models ####
modelsummary(
  list(
    "Latitude"      = model_linear_lat,
    "Turbidity"     = model_linear_turbidity,
    "Depth"         = model_linear_depth,
    "DHW var."      = model_linear_dhw_std,
    "SST var."  = model_linear_hotspot_std,
    "SST warm." = model_linear_hotspot_warming,
    #"Lat+Depth"     = model_linear_lat_depth,
    "Lat+SST var."   = model_linear_lat_hotspot,
    "MPA"           = model_linear_mpa
  ),
  vcov = list(      # Pass vcov objects separately here
    model_linear_lat_vcov,
    model_linear_turbidity_vcov,
    model_linear_depth_vcov,
    model_linear_dhw_std_vcov,
    model_linear_hotspot_std_vcov,
    model_linear_hotspot_warming_vcov,
    #model_linear_lat_depth_vcov,
    model_linear_lat_hotspot_vcov,
    model_linear_mpa_vcov
  ),
  stars     = c('*' = 0.05, '**' = 0.01),
  fmt       = 2,
  statistic = "({std.error})",
  coef_omit = "Intercept",
  coef_rename = c(
    "dhw"                       = "DHW",
    "dhw:abs_lat"               = "DHW × Latitude",
    "dhw:Turbidity"             = "DHW × Turbidity",
    "dhw:Depth_m"               = "DHW × Depth",
    "dhw:DHW_1985.2005_std"     = "DHW × DHW var",
    "dhw:hotspot_1985.2005_std" = "DHW × SST var",
    "dhw:hotspot_warming"       = "DHW × SST warm",
    "dhw:in_mpaTRUE"            = "DHW × MPA"
  ),
  gof_omit  = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes     = list(
    "Standard errors in parentheses, using Conley 200km errors.",
    "All models include site and year fixed effects."
  ),
  title  = "Interaction models (DHW with site characteristics)",
  output = "flextable"
) %>%
  nature_table_style() %>%
  save_as_docx(path = "table_interactions_nature.docx")


#### Lag sensitivity ####
modelsummary(
  list(
    "−3" = model_lagminus3, "−2" = model_lagminus2, "−1" = model_lagminus1,
    "0"  = model_lag0,
    "+1" = model_lag1, "+2" = model_lag2, "+3" = model_lag3,
    "+4" = model_lag4, "+5" = model_lag5
  ),
  vcov = list(
    model_lagminus3_vcov, model_lagminus2_vcov, model_lagminus1_vcov,
    model_lag0_vcov,
    model_lag1_vcov, model_lag2_vcov, model_lag3_vcov,
    model_lag4_vcov, model_lag5_vcov
  ),
  stars     = c('*' = 0.05, '**' = 0.01),
  fmt       = 2,
  statistic = "({std.error})",
  coef_omit = "Intercept",
  coef_rename = c(
    "dhw_lag.3" = "DHW", "dhw_lag.2" = "DHW", "dhw_lag.1" = "DHW",
    "dhw_lag0"  = "DHW",
    "dhw_lag1"  = "DHW", "dhw_lag2"  = "DHW", "dhw_lag3"  = "DHW",
    "dhw_lag4"  = "DHW", "dhw_lag5"  = "DHW"
  ),
  gof_omit  = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes     = list(
    "Standard errors in parentheses, using Conley 200km errors.",
    "All models include site and year fixed effects.",
    "Column headers show lag in months relative to bleaching observation."
  ),
  title  = "Temporal lag sensitivity analysis",
  output = "flextable"
) %>%
  nature_table_style() %>%
  save_as_docx(path = "table_lags_nature.docx")

#### Polynomial sensitivity ####
modelsummary(
  list(
    "Linear"        = model_linear_lat,
    "Quadratic"     = model_quadratic,
    "Quadratic+Lat" = model_quadratic_lat,
    "Cubic"         = model_cubic,
    "Cubic+Lat"     = model_cubic_lat
  ),
  vcov = list(
    model_linear_lat_vcov,
    model_quadratic_vcov,
    model_quadratic_lat_vcov,
    model_cubic_vcov,
    model_cubic_lat_vcov
  ),
  stars     = c('*' = 0.05, '**' = 0.01),
  fmt       = 2,
  statistic = "({std.error})",
  coef_omit = "Intercept",
  coef_rename = c(
    "dhw"              = "DHW",
    "I(dhw^2)"         = "DHW²",
    "I(dhw^3)"         = "DHW³",
    "dhw:abs_lat"      = "DHW × Lat",
    "I(dhw^2):abs_lat" = "DHW² × Lat",
    "I(dhw^3):abs_lat" = "DHW³ × Lat"
  ),
  gof_omit = "R2$|R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes    = list(
    "* p < 0.05, ** p < 0.01",
    "Standard errors in parentheses, using Conley 200km errors.",
    "All models include site and year fixed effects."
  ),
  title  = "Polynomial DHW sensitivity analysis",
  output = "flextable"
) %>%
  nature_table_style() %>%
  save_as_docx(path = "table_polynomial_sensitivity.docx")






model_linear_nocluster <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_ecoregion <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  cluster = ~Ecoregion_Name, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_site <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  cluster = ~Site_ID, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_nocluster_vcov <- conley_vcov(model_linear_nocluster)
model_linear_ecoregion_vcov <- conley_vcov(model_linear_ecoregion)
model_linear_site_vcov <- conley_vcov(model_linear_site)


modelsummary(
  list(
    "Conley 200km"  = model_linear, 
    "Ecoregion" = model_linear_ecoregion,
    "Site"      = model_linear_site,
    "None"      = model_linear_nocluster
    
  ),
  vcov = list(
    model_linear_vcov,
    ~Ecoregion_Name,           # use model's own ecoregion clustering
    ~Site_ID,                  # use model's own site clustering
    "iid"                      # no clustering
  ),
  stars     = c('*' = 0.05, '**' = 0.01),
  fmt       = 2,
  statistic = "({std.error})",
  coef_omit = "Intercept",
  coef_rename = c("dhw" = "DHW"),
  gof_omit  = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes     = list(
    "Standard errors in parentheses.",
    "All models include site and year fixed effects.",
    "Clustering level varies by column as indicated."
  ),
  title  = "Standard error clustering sensitivity",
  output = "flextable"
) %>%
  nature_table_style() %>%
  save_as_docx(path = "table_clustering_nature.docx")






#### seasonality #### 
#### Model with seasonality FE ####
model_linear_lat_season <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)


#### Compare to baseline model ####
model_linear_lat <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)


model_linear_lat_vcov        <- conley_vcov(model_linear_lat)
model_linear_lat_season_vcov <- conley_vcov(model_linear_lat_season)

#### Comparison table ####
modelsummary(
  list(
    "Baseline (Site + Year FE)"            = model_linear_lat,
    "+ Seasonality (Ecoregion×Month FE)"   = model_linear_lat_season
  ),
  vcov = list(
    model_linear_lat_vcov,
    model_linear_lat_season_vcov
  ),
  stars     = c('*' = 0.05, '**' = 0.01),
  fmt       = 2,
  statistic = "({std.error})",
  coef_omit = "Intercept",
  coef_rename = c(
    "dhw"         = "DHW",
    "dhw:abs_lat" = "DHW × |Latitude|"
  ),
  gof_omit  = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes     = list(
    "Standard errors in parentheses, Conley (200km) spatial correction.",
    "Both models include Site_ID and Date_Year fixed effects.",
    "Ecoregion×Month FE absorbs seasonal differences in bleaching timing between hemispheres."
  ),
  title  = "Sensitivity to seasonality (ecoregion-month fixed effects)",
  output = "flextable"
) %>%
  nature_table_style() %>%
  save_as_docx(path = "table_seasonality_sensitivity.docx")









#### DHW specification sensitivity (DHW vs DHW_adj) ####

All_Bleaching_Events_Data_AllDHW_test <- All_Bleaching_Events_Data_AllDHW
All_Bleaching_Events_Data_AllDHW_test <- left_join(All_Bleaching_Events_Data_AllDHW_test, DF[, c(1, 2, 3, 4, 5, 6, 18)])
#### Fit base + adjusted DHW models (Conley SEs computed post-estimation) ####
model_linear_dhwspec <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month,
  data = All_Bleaching_Events_Data_AllDHW_Max)
model_linear_dhwspec_vcov <- conley_vcov(model_linear_dhwspec)

model_linear_lat_dhwspec <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  data = All_Bleaching_Events_Data_AllDHW_Max)
model_linear_lat_dhwspec_vcov <- conley_vcov(model_linear_lat_dhwspec)

model_linear_adj <- feols(
  Percent_Bleached ~ DHW_adj | Site_ID + Date_Year + Ecoregion_Month,
  data = All_Bleaching_Events_Data_AllDHW_Max)
model_linear_adj_vcov <- conley_vcov(model_linear_adj)

model_linear_adj_lat <- feols(
  Percent_Bleached ~ DHW_adj + DHW_adj:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  data = All_Bleaching_Events_Data_AllDHW_Max)
model_linear_adj_lat_vcov <- conley_vcov(model_linear_adj_lat)

#### Table: DHW specification sensitivity ####
modelsummary(
  list(
    "DHW"               = model_linear_dhwspec,
    "DHW + Lat."        = model_linear_lat_dhwspec,
    "DHW (adj.)"        = model_linear_adj,
    "DHW (adj.) + Lat." = model_linear_adj_lat
  ),
  vcov = list(
    model_linear_dhwspec_vcov,
    model_linear_lat_dhwspec_vcov,
    model_linear_adj_vcov,
    model_linear_adj_lat_vcov
  ),
  stars     = c('*' = 0.05, '**' = 0.01),
  fmt       = 2,
  statistic = "({std.error})",
  coef_omit = "Intercept",
  coef_rename = c(
    "DHW"             = "DHW",
    "DHW_adj"         = "DHW (adjusted)",
    "DHW:abs_lat"     = "DHW × |Latitude|",
    "DHW_adj:abs_lat" = "DHW (adj.) × |Latitude|"
  ),
  gof_omit  = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes     = list(
    "Standard errors in parentheses, Conley (200km) spatial correction.",
    "All models include site and year fixed effects.",
    "DHW (adjusted) uses the thermal anomaly adjustment from Ainsworth et al. (2016)."
  ),
  title  = "Sensitivity to DHW specification",
  output = "flextable"
) %>%
  nature_table_style() %>%
  save_as_docx(path = "table_sensitivity_dhw_spec.docx")
