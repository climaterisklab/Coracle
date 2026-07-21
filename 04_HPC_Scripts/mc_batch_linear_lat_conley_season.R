## author: Puja Pande
## date: June 2026
## description: MC BATCH — Linear FE model with Conley parametric bootstrap
##              Coefficients drawn from Conley vcov via mvrnorm
##              run as SLURM array (batches 1-84, one per CF file)

library(dplyr)
library(ncdf4)
library(parallel)

n_cores    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 4))
batch      <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", 1))
SCRATCH    <- "/gpfs/scratch/bsc32/bsc963510"
CF_DIR     <- "/esarchive/Coracle/climate_data/DHW/count"
OUT_PREFIX <- "linear_lat_conley_season"
OUT_DIR    <- file.path(SCRATCH, paste0(OUT_PREFIX, "_nd"))
n_samples  <- 10
n_time     <- 480

cat("Batch (CF file)", batch, "| Using", n_cores, "cores\n")

prediction_grid    <- readRDS(file.path(SCRATCH, paste0("prediction_grid_", OUT_PREFIX, ".rds")))
boot_results_clean <- readRDS(file.path(SCRATCH, paste0("boot_results_", OUT_PREFIX, ".rds")))

counterfactual_files <- sort(list.files(CF_DIR,
  pattern = "DHW_count_source_mm_[0-9]+\\.nc$", full.names = TRUE))

nc_tmp    <- nc_open(counterfactual_files[1])
scenarios <- as.character(ncvar_get(nc_tmp, "firm"))
site_ids  <- as.integer(ncvar_get(nc_tmp, "Site_ID"))
nc_close(nc_tmp)
n_sites     <- length(site_ids)
n_scenarios <- length(scenarios)
cat("Sites:", n_sites, "| Scenarios:", n_scenarios, "\n")

cf_file      <- counterfactual_files[batch]
boot_indices <- ((batch - 1) * n_samples + 1):(batch * n_samples)
boot_indices <- ((boot_indices - 1) %% nrow(boot_results_clean)) + 1
cat("CF file:", basename(cf_file), "\n")

pg_sites   <- site_ids[site_ids %in% unique(prediction_grid$Site_ID)]
pg_years   <- sort(unique(prediction_grid$Date_Year))
n_pg_sites <- length(pg_sites)
n_pg_years <- length(pg_years)
cat("Prediction grid sites:", n_pg_sites, "| Years:", n_pg_years, "\n")

nc_cf    <- nc_open(cf_file)
site_dim <- as.integer(ncvar_get(nc_cf, "Site_ID"))
firm_dim <- as.character(ncvar_get(nc_cf, "firm"))
nc_close(nc_cf)

cl <- makeCluster(n_cores)
clusterExport(cl, c("prediction_grid", "boot_results_clean", "cf_file", "boot_indices",
                     "scenarios", "pg_sites", "pg_years", "n_pg_sites", "n_pg_years",
                     "site_dim", "firm_dim", "n_samples"))
clusterEvalQ(cl, library(ncdf4))

batch_results <- parLapply(cl, 1:n_samples, function(s) {
  boot_idx <- boot_indices[s]
  beta_dhw <- boot_results_clean$dhw[boot_idx]
  beta_lat <- boot_results_clean$`dhw:abs_lat`[boot_idx]

  pred_factual_all <- ifelse(!is.na(prediction_grid$dhw),
    beta_dhw * prediction_grid$dhw +
      beta_lat * prediction_grid$dhw * prediction_grid$abs_lat,
    NA_real_)

  max_pred     <- matrix(NA_real_,    nrow=n_pg_sites, ncol=n_pg_years)
  max_time_idx <- matrix(NA_integer_, nrow=n_pg_sites, ncol=n_pg_years)

  for (si in seq_along(pg_sites)) {
    for (yi in seq_along(pg_years)) {
      rows <- which(prediction_grid$Site_ID  == pg_sites[si] &
                    prediction_grid$Date_Year == pg_years[yi] &
                    !is.na(pred_factual_all))
      if (length(rows) == 0) next
      best_row             <- rows[which.max(pred_factual_all[rows])]
      max_pred[si, yi]     <- pred_factual_all[best_row]
      max_time_idx[si, yi] <- prediction_grid$time_idx[best_row]
    }
  }

  nc_cf_open <- nc_open(cf_file)
  cf_preds <- lapply(scenarios, function(scenario) {
    firm_idx <- which(firm_dim == scenario)
    cf_mat   <- matrix(NA_real_, nrow=n_pg_sites, ncol=n_pg_years)
    for (si in seq_along(pg_sites)) {
      site_pos <- which(site_dim == pg_sites[si])
      if (length(site_pos) == 0) next
      lat <- prediction_grid$abs_lat[prediction_grid$Site_ID == pg_sites[si]][1]
      for (yi in seq_along(pg_years)) {
        ti <- max_time_idx[si, yi]
        if (is.na(ti) || length(firm_idx) == 0) next
        cf_dhw <- tryCatch(
          ncvar_get(nc_cf_open, "DHW", start=c(ti, site_pos, s, firm_idx), count=c(1,1,1,1)),
          error = function(e) NA_real_)
        if (!is.na(cf_dhw))
          cf_mat[si, yi] <- beta_dhw * cf_dhw + beta_lat * cf_dhw * lat
      }
    }
    cf_mat
  })
  names(cf_preds) <- scenarios
  nc_close(nc_cf_open)
  list(factual = max_pred, cf = cf_preds)
})
stopCluster(cl)

dir.create(OUT_DIR, showWarnings = FALSE)
out_file <- file.path(OUT_DIR, sprintf("Predictions_%s_batch%03d.nc", OUT_PREFIX, batch))

dim_site     <- ncdim_def("site",     "", 1:n_pg_sites,  create_dimvar=FALSE)
dim_year     <- ncdim_def("year",     "", 1:n_pg_years,  create_dimvar=FALSE)
dim_sample   <- ncdim_def("sample",   "", 1:n_samples,   create_dimvar=FALSE)
dim_scenario <- ncdim_def("scenario", "", 1:n_scenarios, create_dimvar=FALSE)

var_lon     <- ncvar_def("lon", "degrees_east",  list(dim_site), -9999.0, prec="float")
var_lat     <- ncvar_def("lat", "degrees_north", list(dim_site), -9999.0, prec="float")
var_factual <- ncvar_def("pred_factual", "%",
                          list(dim_site, dim_year, dim_sample), -9999.0, prec="float",
                          compression=4, longname="Factual bleaching at max month")
var_cf      <- ncvar_def("pred_cf", "%",
                          list(dim_site, dim_year, dim_sample, dim_scenario), -9999.0, prec="float",
                          compression=4, longname="Counterfactual bleaching at max month")

nc_out <- nc_create(out_file, list(var_lon, var_lat, var_factual, var_cf))
ncatt_put(nc_out, 0, "scenarios",         paste(scenarios, collapse=","))
ncatt_put(nc_out, 0, "site_ids",          paste(pg_sites,  collapse=","))
ncatt_put(nc_out, 0, "years",             paste(pg_years,  collapse=","))
ncatt_put(nc_out, 0, "model",             OUT_PREFIX)
ncatt_put(nc_out, 0, "batch",             batch)
ncatt_put(nc_out, 0, "bootstrap_indices", paste(boot_indices, collapse=","))

pg_lons <- prediction_grid$Longitude_Degrees[match(pg_sites, prediction_grid$Site_ID)]
pg_lats <- prediction_grid$Latitude_Degrees[match(pg_sites, prediction_grid$Site_ID)]
ncvar_put(nc_out, var_lon, pg_lons)
ncvar_put(nc_out, var_lat, pg_lats)

for (s in 1:n_samples) {
  ncvar_put(nc_out, var_factual, batch_results[[s]]$factual,
            start=c(1,1,s), count=c(n_pg_sites, n_pg_years, 1))
  for (sc_idx in seq_along(scenarios)) {
    ncvar_put(nc_out, var_cf, batch_results[[s]]$cf[[scenarios[sc_idx]]],
              start=c(1,1,s,sc_idx), count=c(n_pg_sites, n_pg_years, 1, 1))
  }
}
nc_close(nc_out)
cat("Saved:", out_file, "\n")
cat("Batch", batch, "complete\n")
