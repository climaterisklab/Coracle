library(fixest)
library(dplyr)

All_Bleaching_Events_Data_AllDHW <- read.csv("~/Library/CloudStorage/Dropbox/Coracle/Merged_Mermaid_Panel_Bleaching_Data.csv")

# Step 1: Get annual max
annual_max_bleaching <- All_Bleaching_Events_Data_AllDHW %>%
  group_by(Site_ID, Date_Year) %>%
  summarise(
    Max_Bleaching = max(Percent_Bleached, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: For each site, explicitly create lags by looking up prior years
annual_max_bleaching <- annual_max_bleaching %>%
  arrange(Site_ID, Date_Year) %>%
  group_by(Site_ID) %>%
  mutate(
    # For Lag1: Look up Max_Bleaching from (current year - 1)
    Lag1_Max_Bleaching = Max_Bleaching[match(Date_Year - 1, Date_Year)],
    # For Lag2: Look up Max_Bleaching from (current year - 2)
    Lag2_Max_Bleaching = Max_Bleaching[match(Date_Year - 2, Date_Year)],
    # For Lag3: Look up Max_Bleaching from (current year - 3)
    Lag3_Max_Bleaching = Max_Bleaching[match(Date_Year - 3, Date_Year)]
  ) %>%
  ungroup()

# Verify this worked
cat("\n=== Lag Verification ===\n")
cat("Site 4 lags (should show NA when prior year doesn't exist):\n")
annual_max_bleaching %>%
  filter(Site_ID == 4) %>%
  arrange(Date_Year) %>%
  print(n = 20)

cat("\nSite 4 lags:\n")
annual_max_bleaching %>%
  filter(Site_ID == 4) %>%
  arrange(Date_Year) %>%
  print(n = 20)

# Step 3: Filter to only keep sites with at least some consecutive years
sites_with_consecutive <- annual_max_bleaching %>%
  filter(!is.na(Lag1_Max_Bleaching)) %>%
  distinct(Site_ID)

cat("\nSites with at least one consecutive year pair:", nrow(sites_with_consecutive), "\n")

annual_max_bleaching_filtered <- annual_max_bleaching %>%
  semi_join(sites_with_consecutive, by = "Site_ID")

cat("Site-years after filtering:", nrow(annual_max_bleaching_filtered), "\n\n")

# Step 4: Merge back to original data
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  select(-starts_with("Lag")) %>%  # Remove any old lag columns
  left_join(
    annual_max_bleaching_filtered %>% select(Site_ID, Date_Year, starts_with("Lag")),
    by = c("Site_ID", "Date_Year")
  )


#### Model Testing ####

# Model 1: Without lag (baseline)
model_no_lag <- fepois(
  Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)

# Model 2: With 1-year lag
model_with_lag <- fepois(
  Percent_Bleached ~ log1p(dhw) + Lag1_Max_Bleaching | Site_ID + Date_Year,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)

# Model 3: With multiple lags
model_multi_lag <- fepois(
  Percent_Bleached ~ log1p(dhw) + Lag1_Max_Bleaching + Lag2_Max_Bleaching + Lag3_Max_Bleaching | Site_ID + Date_Year,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)


#### Results ####

cat("\n\n=== Model Comparison ===\n")
etable(model_no_lag, model_with_lag, model_multi_lag)


#### Interpretation ####

cat("\n\n=== Bleaching Dependence Test ===\n")

lag1_coef <- coef(model_with_lag)["Lag1_Max_Bleaching"]
lag1_se <- se(model_with_lag)["Lag1_Max_Bleaching"]
lag1_pval <- pvalue(model_with_lag)["Lag1_Max_Bleaching"]

cat("Lag 1 coefficient:", round(lag1_coef, 4), "\n")
cat("Standard error:", round(lag1_se, 4), "\n")
cat("P-value:", format.pval(lag1_pval, digits = 3), "\n\n")

if(!is.na(lag1_pval) && lag1_pval < 0.05) {
  if(lag1_coef > 0) {
    cat("Result: POSITIVE dependence - sensitization/damage\n\n")
    cat("Interpretation: Each 1% increase in last year's bleaching\n")
    cat("                → ", round((exp(lag1_coef) - 1) * 100, 2), "% increase in current bleaching\n", sep = "")
  } else {
    cat("Result: NEGATIVE dependence - adaptation/selection\n")
  }
} else {
  cat("Result: NO significant dependence - independent events\n")
}


#### Data Summary ####

cat("\n\n=== Data Summary ===\n")
cat("Total observations:", nrow(All_Bleaching_Events_Data_AllDHW), "\n")
cat("With Lag1 data:", sum(!is.na(All_Bleaching_Events_Data_AllDHW$Lag1_Max_Bleaching)), "\n")
cat("With Lag2 data:", sum(!is.na(All_Bleaching_Events_Data_AllDHW$Lag2_Max_Bleaching)), "\n")
cat("With Lag3 data:", sum(!is.na(All_Bleaching_Events_Data_AllDHW$Lag3_Max_Bleaching)), "\n\n")

# Show example - should now have proper NAs
cat("Example for Site 4 (should show NAs for non-consecutive years):\n")
All_Bleaching_Events_Data_AllDHW %>%
  filter(Site_ID == 4) %>%
  select(Site_ID, Date_Year, Percent_Bleached, 
         Lag1_Max_Bleaching, Lag2_Max_Bleaching, Lag3_Max_Bleaching) %>%
  distinct(Date_Year, .keep_all = TRUE) %>%
  arrange(Date_Year) %>%
  print()

