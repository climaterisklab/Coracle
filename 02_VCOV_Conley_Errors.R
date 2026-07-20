### 07_linear_and_binned_conley_models.R
### Produces all four model output files needed for the panel B plotting script:
###   1. linear_lat_conley200km_coefs.csv   — 2 rows (DHW, DHW:abs_lat)
###   2. linear_lat_conley200km_vcov.rds    — 2×2 matrix
###   3. FE_binned_DHW_conley200km_coefs.csv — one row per bin term
###   4. FE_binned_DHW_conley200km_vcov.rds  — N_bins × N_bins matrix
###
### IMPORTANT: each model's vcov is saved under its OWN variable name to
### prevent the silent cross-contamination bug where running both models in
### the same session overwrites the wrong vcov file.

library(dplyr)
library(fixest)

#### ── DATA ──────────────────────────────────────────────────────────────── ####
All_Bleaching_Events_Data_AllDHW <- read.csv("path/to/Coracle/Datasets/Final_Panel_Data_allDHW_01July.csv")
All_Bleaching_Events_Data_AllDHW$Proportion_Bleached <- All_Bleaching_Events_Data_AllDHW$Percent_Bleached / 100
All_Bleaching_Events_Data_AllDHW$abs_lat             <- abs(All_Bleaching_Events_Data_AllDHW$Latitude_Degrees)
All_Bleaching_Events_Data_AllDHW$mass_bleaching      <- ifelse(All_Bleaching_Events_Data_AllDHW$Percent_Bleached >= 30, 1, 0)
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Ecoregion_Month = interaction(Ecoregion_Name, Date_Month, drop = TRUE)
  )
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Proportion_Bleached     = Percent_Bleached / 100,
    Proportion_Bleached_adj = (Proportion_Bleached * (n() - 1) + 0.5) / n(),
    any_bleaching           = ifelse(Percent_Bleached > 0, 1, 0)
  )


All_Bleaching_Events_Data_AllDHW_Max <- read.csv("path/to/data/Final_Combined_Data_all_levels_01July2026_extravars.csv")
All_Bleaching_Events_Data_AllDHW_Max <- left_join(All_Bleaching_Events_Data_AllDHW, All_Bleaching_Events_Data_AllDHW_Max[, c(1:8, 17:34)], by = c("Site_ID", "Latitude_Degrees", "Longitude_Degrees", "Date_Year", "Date_Month", "Date_Day", "Ecoregion_Name", "Percent_Bleached"))

DF <- All_Bleaching_Events_Data_AllDHW_Max

DF$abs_lat <- abs(DF$Latitude_Degrees)

DF_model <- DF %>%
  filter(!is.na(DHW), !is.na(Percent_Bleached), !is.na(Site_ID), !is.na(Date_Year)) %>%
  mutate(
    Site_ID   = as.factor(Site_ID),
    Date_Year = as.factor(Date_Year)
  )

#### ── HELPER: strip fixest_vcov class and save ───────────────────────────── ####

clean_vcov <- function(vcov_con) {
  m <- unclass(vcov_con)
  class(m) <- "matrix"
  m
}

save_model_outputs <- function(model, vcov_con, coef_file, vcov_file,
                               coef_value_col = "x") {
  coefs    <- coef(model)
  coef_df  <- data.frame(term = names(coefs))
  coef_df[[coef_value_col]] <- as.numeric(coefs)
  write.csv(coef_df, coef_file, row.names = FALSE)
  
  vcov_mat <- clean_vcov(vcov_con)
  saveRDS(vcov_mat, vcov_file)
  
  cat(sprintf("\nSaved: %s\n", coef_file))
  cat(sprintf("Saved: %s\n", vcov_file))
  cat(sprintf("  vcov dim: %d x %d\n", nrow(vcov_mat), ncol(vcov_mat)))
  cat(sprintf("  row/col names match: %s\n",
              identical(rownames(vcov_mat), colnames(vcov_mat))))
  
  invisible(vcov_mat)
}

#### ── (i) LINEAR FE — Conley (200km, spherical) ─────────────────────────── ####

cat("\n==== Fitting linear FE (DHW + DHW:abs_lat) ====\n")

model_linear_lat <- feols(
  Percent_Bleached ~ DHW + DHW:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = DF_model
)
summary(model_linear_lat)

cat("Computing Conley vcov for linear model...\n")
vcov_con_linear <- vcov_conley(
  model_linear_lat,
  lat      = "Latitude_Degrees",
  lon      = "Longitude_Degrees",
  cutoff   = 200,
  distance = "spherical"
)
summary(model_linear_lat, vcov_con_linear)

linear_vcov_mat <- save_model_outputs(
  model         = model_linear_lat,
  vcov_con      = vcov_con_linear,
  coef_file     = "linear_lat_conley200km_coefs.csv",
  vcov_file     = "linear_lat_conley200km_vcov.rds",
  coef_value_col = "x"
)

stopifnot(
  "linear vcov must be 2x2"          = all(dim(linear_vcov_mat) == c(2, 2)),
  "linear vcov must not contain bins" = !any(grepl("DHW_bin", rownames(linear_vcov_mat)))
)

#### ── (ii) FE-BINNED DHW — Conley (200km, spherical) ───────────────────── ####

cat("\n==== Fitting FE-binned DHW model ====\n")

# Bins: 0.5-unit below DHW=6, 1-unit from 6-12, tail [12,Inf)
bin_edges <- c(seq(0, 6, by = 0.5), 8, Inf)
DF_model$DHW_bin <- cut(DF_model$DHW, breaks = bin_edges, right = FALSE)

cat("N per DHW bin:\n")
print(table(DF_model$DHW_bin))

ref_bin <- levels(DF_model$DHW_bin)[1]
cat(sprintf("\nRef bin: %s\n", ref_bin))

m_binned <- feols(
  Percent_Bleached ~ i(DHW_bin, ref = ref_bin) +
    i(DHW_bin, abs_lat, ref = ref_bin) | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = DF_model
)
summary(m_binned)

cat("Computing Conley vcov for binned model...\n")
vcov_con_binned <- vcov_conley(
  m_binned,
  lat      = "Latitude_Degrees",
  lon      = "Longitude_Degrees",
  cutoff   = 200,
  distance = "spherical"
)
summary(m_binned, vcov_con_binned)

cat(sprintf("%d of %d observations retained (%.1f%%)\n",
            nobs(m_binned), nrow(DF_model),
            100 * nobs(m_binned) / nrow(DF_model)))

binned_vcov_mat <- save_model_outputs(
  model          = m_binned,
  vcov_con       = vcov_con_binned,
  coef_file      = "FE_binned_DHW_conley200km_coefs.csv",
  vcov_file      = "FE_binned_DHW_conley200km_vcov.rds",
  coef_value_col = "estimate"
)

stopifnot(
  "binned vcov must contain DHW_bin terms" =
    any(grepl("DHW_bin", rownames(binned_vcov_mat)))
)
