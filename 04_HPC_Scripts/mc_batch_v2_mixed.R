## author: Puja Pande
## date: November 2026
## description: MC BATCH — v2_mixed impact model (province-specific linear / sigmoid)
##              Parameter draws from bleaching_impact_model/parameter_draws.csv (1000 per model)
##              Factual: observed annual max DHW; counterfactual: annual max DHW by emission source
##              run as SLURM array (batches 1-84, one per CF seed file; draws (batch-1)*10 + 1:10)

library(ncdf4)

batch      <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", 1))
SCRATCH    <- "/gpfs/scratch/bsc32/bsc963510"
V2_DIR     <- "/esarchive/Coracle/climate_data/v2_sites"
IMPACT_DIR <- "/esarchive/Coracle/bleaching_impact_model"
OUT_PREFIX <- Sys.getenv("OUT_PREFIX", "v2_mixed")
OUT_DIR    <- file.path(SCRATCH, paste0(OUT_PREFIX, "_nd"))
n_samples  <- 10
years_all  <- 1985:2024
# provinces without a bleaching model (Tropical Eastern Pacific) are absent from the sites file -> NaN
# set EXCLUDE_PROVINCES="Red_Sea,Central_Pacific" to mirror the forecast paper (impact.py NO_BLEACHING)
exclude    <- Sys.getenv("EXCLUDE_PROVINCES", "")
exclude    <- if (nzchar(exclude)) strsplit(exclude, ",")[[1]] else character(0)

cat("Batch (CF seed file)", batch, "| model:", OUT_PREFIX, "| excluded provinces:",
    if (length(exclude)) paste(exclude, collapse = ",") else "none", "\n")

#### sites -> province (one model per cell, coord_id) ####
sites <- read.csv(file.path(V2_DIR, "attribution_sites_v2.csv"))
cells <- unique(sites[, c("coord_id", "province", "model_form", "coord_lat", "coord_lon")])
stopifnot("a coord_id maps to more than one province / form" = !anyDuplicated(cells$coord_id))

#### counterfactual file for this batch ####
cf_files <- sort(list.files(file.path(V2_DIR, "DHW_count"),
                            pattern = "^DHW_count_source_ann_[0-9]+\\.nc$", full.names = TRUE))
cf_files <- cf_files[order(as.integer(sub(".*_ann_([0-9]+)\\.nc$", "\\1", cf_files)))]
stopifnot(length(cf_files) == 84)
cf_file <- cf_files[batch]
cat("CF file:", basename(cf_file), "\n")

nc_cf     <- nc_open(cf_file)
cf_ids    <- as.integer(ncvar_get(nc_cf, "coord_id"))
scenarios <- as.character(ncvar_get(nc_cf, "firm"))
cf        <- ncvar_get(nc_cf, "DHW")                    # [coord_id, sample, firm, time]
nc_close(nc_cf)
stopifnot(dim(cf)[2] == n_samples, dim(cf)[4] == length(years_all))
n_scenarios <- length(scenarios)

#### observed (factual) annual max DHW ####
nc_obs  <- nc_open(file.path(V2_DIR, "DHW_obs_ann_1985-2024.nc"))
obs_ids <- as.integer(ncvar_get(nc_obs, "coord_id"))
obs     <- ncvar_get(nc_obs, "DHW")                     # [coord_id, time]
nc_close(nc_obs)
obs     <- obs[match(cf_ids, obs_ids), , drop = FALSE]  # align to the CF cells
stopifnot(!anyNA(match(cf_ids, obs_ids)))

cells <- cells[match(cf_ids, cells$coord_id), ]
stopifnot(!anyNA(cells$coord_id))
cat("Cells:", length(cf_ids), "| Years:", length(years_all), "| Scenarios:", n_scenarios, "\n")

#### parameter draws: draw (batch-1)*10 + s for realisation s ####
draws        <- read.csv(file.path(IMPACT_DIR, "parameter_draws.csv"))
draw_indices <- (batch - 1) * n_samples + 1:n_samples
stopifnot(max(draw_indices) <= max(draws$draw))
models <- setdiff(unique(cells$province), exclude)
stopifnot("province without parameter draws" = all(models %in% draws$model))

f_model <- function(D, p) {
  if (p$form == "linear") return(p$b * D)
  A  <- min(max(p$A, 0), 100); s <- max(p$s, 0.1)
  A / (1 + exp(-(D - p$x0) / s)) - A / (1 + exp(p$x0 / s))   # f(D) - f(0)
}

n_cells <- length(cf_ids)
pred_factual <- array(NA_real_, c(n_cells, length(years_all), n_samples))
pred_cf      <- array(NA_real_, c(n_cells, length(years_all), n_samples, n_scenarios))

for (m in models) {
  idx <- which(cells$province == m)
  for (s in 1:n_samples) {
    p <- draws[draws$model == m & draws$draw == draw_indices[s], ]
    stopifnot(nrow(p) == 1)
    pred_factual[idx, , s] <- f_model(obs[idx, , drop = FALSE], p)
    for (sc in seq_len(n_scenarios)) {
      pred_cf[idx, , s, sc] <- f_model(matrix(cf[idx, s, sc, ], nrow = length(idx)), p)
    }
  }
  cat("  ", m, ":", length(idx), "cells\n")
}

#### write ####
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(OUT_DIR, sprintf("Predictions_%s_batch%03d.nc", OUT_PREFIX, batch))

dim_site     <- ncdim_def("site",     "", 1:n_cells,       create_dimvar = FALSE)
dim_year     <- ncdim_def("year",     "", 1:length(years_all), create_dimvar = FALSE)
dim_sample   <- ncdim_def("sample",   "", 1:n_samples,     create_dimvar = FALSE)
dim_scenario <- ncdim_def("scenario", "", 1:n_scenarios,   create_dimvar = FALSE)

var_lon     <- ncvar_def("lon", "degrees_east",  list(dim_site), -9999.0, prec = "float")
var_lat     <- ncvar_def("lat", "degrees_north", list(dim_site), -9999.0, prec = "float")
var_factual <- ncvar_def("pred_factual", "%",
                         list(dim_site, dim_year, dim_sample), -9999.0, prec = "float",
                         compression = 4, longname = "Factual bleaching at annual max DHW")
var_cf      <- ncvar_def("pred_cf", "%",
                         list(dim_site, dim_year, dim_sample, dim_scenario), -9999.0, prec = "float",
                         compression = 4, longname = "Counterfactual bleaching at annual max DHW")

nc_out <- nc_create(out_file, list(var_lon, var_lat, var_factual, var_cf))
ncatt_put(nc_out, 0, "scenarios",         paste(scenarios, collapse = ","))
ncatt_put(nc_out, 0, "site_ids",          paste(cf_ids, collapse = ","))   # coord_id (grid cell)
ncatt_put(nc_out, 0, "provinces",         paste(cells$province, collapse = ","))
ncatt_put(nc_out, 0, "years",             paste(years_all, collapse = ","))
ncatt_put(nc_out, 0, "model",             OUT_PREFIX)
ncatt_put(nc_out, 0, "batch",             batch)
ncatt_put(nc_out, 0, "cf_file",           basename(cf_file))
ncatt_put(nc_out, 0, "bootstrap_indices", paste(draw_indices, collapse = ","))

ncvar_put(nc_out, var_lon, cells$coord_lon)
ncvar_put(nc_out, var_lat, cells$coord_lat)
pred_factual[is.na(pred_factual)] <- -9999.0
pred_cf[is.na(pred_cf)]           <- -9999.0
ncvar_put(nc_out, var_factual, pred_factual)
ncvar_put(nc_out, var_cf,      pred_cf)
nc_close(nc_out)
cat("Saved:", out_file, "\n")
cat("Batch", batch, "complete\n")
