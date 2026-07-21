## author: Puja Pande
## date: June 2026
## description: SETUP — Linear FE model with Conley standard errors
##              Parametric bootstrap using mvrnorm from Conley vcov
##              Prediction grid: all sites x 1985-2024 x 12 months

library(fixest)
library(MASS)
library(dplyr)
library(lubridate)
library(ncdf4)
library(parallel)

n_cores    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 4))
SCRATCH    <- "/gpfs/scratch/bsc32/bsc963510"
NOAA_DIR   <- file.path(SCRATCH, "NOAA_DHW")
OUT_PREFIX <- "linear_lat_conley_season"
cat("Using", n_cores, "cores | Model:", OUT_PREFIX, "\n")

#### data ####
panel_data <- read.csv(file.path(SCRATCH, "Final_Combined_Data_14May2026.csv"))
cat("Panel data rows:", nrow(panel_data), "\n")
cat("Unique sites:", length(unique(panel_data$Site_ID)), "\n")

#### extract DHW for survey rows ####
cat("Extracting DHW for survey rows\n")
rasters <- sort(list.files(NOAA_DIR, pattern = "ct5km_dhw-max_v3.1_.*\\.nc$",
                            full.names = TRUE, recursive = TRUE))
raster_dates  <- as.Date(paste0(sub(".*_(\\d{4})(\\d{2})\\.nc$", "\\1-\\2",
                                     basename(rasters)), "-01"))
raster_dhw_ym <- sprintf("%04d-%02d", year(raster_dates), month(raster_dates))
cat("NOAA rasters:", length(rasters), "\n")

nc_ref  <- nc_open(rasters[1])
nc_lons <- ncvar_get(nc_ref, "lon")
nc_lats <- ncvar_get(nc_ref, "lat")
nc_close(nc_ref)

survey_sites <- panel_data %>%
  group_by(Site_ID) %>% slice(1) %>% ungroup() %>%
  mutate(
    xi = sapply(Longitude_Degrees, function(lon) which.min(abs(nc_lons - lon))),
    yi = sapply(Latitude_Degrees,  function(lat) which.min(abs(nc_lats - lat)))
  )

panel_data <- panel_data %>%
  left_join(survey_sites %>% select(Site_ID, xi, yi), by = "Site_ID") %>%
  mutate(ym = sprintf("%04d-%02d", Date_Year, Date_Month))

unique_yms_survey <- unique(panel_data$ym)

cl <- makeCluster(n_cores)
clusterExport(cl, c("panel_data", "rasters", "raster_dhw_ym", "unique_yms_survey"))
clusterEvalQ(cl, library(ncdf4))

dhw_survey <- parLapply(cl, unique_yms_survey, function(ym) {
  r_idx <- which(raster_dhw_ym == ym)
  if (length(r_idx) == 0) return(list(idx=integer(0), dhw=numeric(0)))
  idx <- which(panel_data$ym == ym)
  nc  <- nc_open(rasters[r_idx[1]])
  ndims <- nc$var$degree_heating_week$ndims
  dhw <- sapply(idx, function(i)
    tryCatch(ncvar_get(nc, "degree_heating_week",
                       start = if(ndims==2) c(panel_data$xi[i], panel_data$yi[i]) else c(panel_data$xi[i], panel_data$yi[i], 1),
                       count = if(ndims==2) c(1,1) else c(1,1,1)),
             error=function(e) NA_real_))
  nc_close(nc)
  list(idx=idx, dhw=dhw)
})
stopCluster(cl)

panel_data$dhw <- NA_real_
for (res in dhw_survey) {
  if (length(res$idx) > 0) panel_data$dhw[res$idx] <- res$dhw
}
panel_data$ym <- NULL
panel_data$xi <- NULL
panel_data$yi <- NULL
cat("Survey DHW extracted. NAs:", sum(is.na(panel_data$dhw)), "/", nrow(panel_data), "\n")

#### model ####
model_data <- panel_data %>%
  filter(!is.na(Percent_Bleached), !is.na(dhw)) %>%
  mutate(abs_lat = abs(Latitude_Degrees),
         Ecoregion_Month = paste0(Ecoregion_Name, "_", Date_Month))
cat("Model data rows:", nrow(model_data), "\n")

model_linear_lat <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = model_data)
cat("Model fitted\n")
print(summary(model_linear_lat))

#### Conley standard errors ####
cat("Computing Conley standard errors (cutoff=200km)\n")
vcov_con <- vcov_conley(model_linear_lat,
                        lat      = "Latitude_Degrees",
                        lon      = "Longitude_Degrees",
                        cutoff   = 200,
                        distance = "spherical")
cat("Conley vcov:\n")
print(vcov_con)
print(summary(model_linear_lat, vcov_con))

# save vcov and coefficients
write.csv(vcov_con,
          file.path(SCRATCH, paste0("vcov_", OUT_PREFIX, ".csv")))
write.csv(as.data.frame(t(model_linear_lat$coefficients)),
          file.path(SCRATCH, paste0("coefs_", OUT_PREFIX, ".csv")),
          row.names = FALSE)

#### parametric bootstrap: 840 samples from Conley vcov ####
set.seed(42)
n_samples <- 840
mu    <- as.numeric(model_linear_lat$coefficients)
sigma <- matrix(as.numeric(vcov_con), nrow = 2)

cat("Drawing", n_samples, "samples from Conley vcov\n")
boot_samples <- MASS::mvrnorm(n = n_samples, mu = mu, Sigma = sigma)
boot_results_clean <- as.data.frame(boot_samples)
colnames(boot_results_clean) <- c("dhw", "dhw:abs_lat")

cat("Samples drawn:", nrow(boot_results_clean), "\n")
cat("dhw range:", range(boot_results_clean$dhw), "\n")

write.csv(boot_results_clean,
          file.path(SCRATCH, paste0("Bootstrap_Coefficients_", OUT_PREFIX, ".csv")),
          row.names = FALSE)
saveRDS(boot_results_clean,
        file.path(SCRATCH, paste0("boot_results_", OUT_PREFIX, ".rds")))
cat("Bootstrap saved\n")

#### prediction grid: all sites x 1985-2024 x 12 months ####
site_info <- panel_data %>%
  select(Site_ID, Latitude_Degrees, Longitude_Degrees, Ecoregion_Name) %>%
  distinct(Site_ID, .keep_all = TRUE)

all_years <- 1985:2024

cat("Building prediction grid\n")
prediction_grid <- expand.grid(
  Site_ID    = unique(site_info$Site_ID),
  Date_Year  = all_years,
  Date_Month = 1:12) %>%
  left_join(site_info, by = "Site_ID") %>%
  mutate(abs_lat = abs(Latitude_Degrees))
cat("Prediction grid:", nrow(prediction_grid), "rows\n")

prediction_grid <- prediction_grid %>%
  mutate(time_idx = (Date_Year - 1985) * 12 + Date_Month)
prediction_grid$ym <- sprintf("%04d-%02d", prediction_grid$Date_Year, prediction_grid$Date_Month)

unique_sites <- site_info %>%
  group_by(Site_ID) %>% slice(1) %>% ungroup() %>%
  mutate(
    xi = sapply(Longitude_Degrees, function(lon) which.min(abs(nc_lons - lon))),
    yi = sapply(Latitude_Degrees,  function(lat) which.min(abs(nc_lats - lat)))
  )
prediction_grid <- prediction_grid %>%
  left_join(unique_sites %>% select(Site_ID, xi, yi), by = "Site_ID")

unique_yms <- unique(prediction_grid$ym)
pg_light   <- prediction_grid[, c("xi", "yi")]

cl <- makeCluster(n_cores)
clusterExport(cl, c("pg_light", "prediction_grid", "rasters", "raster_dhw_ym", "unique_yms"))
clusterEvalQ(cl, library(ncdf4))

dhw_by_ym <- parLapply(cl, unique_yms, function(ym) {
  r_idx <- which(raster_dhw_ym == ym)
  if (length(r_idx) == 0) return(list(idx=integer(0), dhw=numeric(0)))
  idx <- which(prediction_grid$ym == ym)
  nc  <- nc_open(rasters[r_idx[1]])
  ndims <- nc$var$degree_heating_week$ndims
  dhw <- sapply(idx, function(i)
    tryCatch(ncvar_get(nc, "degree_heating_week",
                       start = if(ndims==2) c(pg_light$xi[i], pg_light$yi[i]) else c(pg_light$xi[i], pg_light$yi[i], 1),
                       count = if(ndims==2) c(1,1) else c(1,1,1)),
             error=function(e) NA_real_))
  nc_close(nc)
  list(idx=idx, dhw=dhw)
})
stopCluster(cl)

prediction_grid$dhw <- NA_real_
for (res in dhw_by_ym) {
  if (length(res$idx) > 0) prediction_grid$dhw[res$idx] <- res$dhw
}
prediction_grid$ym <- NULL
prediction_grid$xi <- NULL
prediction_grid$yi <- NULL

cat("Prediction grid DHW extracted. NAs:", sum(is.na(prediction_grid$dhw)),
    "/", nrow(prediction_grid), "\n")
print(prediction_grid %>% filter(Date_Year >= 2023) %>%
      group_by(Date_Year) %>% summarise(n=n(), n_nonNA=sum(!is.na(dhw))))

saveRDS(prediction_grid,
        file.path(SCRATCH, paste0("prediction_grid_", OUT_PREFIX, ".rds")))
cat("Prediction grid saved:", nrow(prediction_grid), "rows\n")
