#### author: Puja Pande
#### date: 20 January 2026
#### description: Integrated uncertainty pipeline for factual/counterfactual coral bleaching attribution - ALL SCENARIOS

#### libraries ####
library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)

#### Load data ####
# Impact model bootstrap results
boot_results_clean <- read.csv("Bootstrap_Coefficients_Impact_Model.csv")

# Bleaching observations with factual DHW
All_Bleaching_Events_Data_AllDHW <- read.csv("Final_Relevant_Scripts/Merged_Mermaid_Panel_Bleaching_Data.csv")

# Load original impact model
library(fixest)
model_negative_binomial <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, 
                                    cluster = ~Ecoregion_Name, 
                                    data = All_Bleaching_Events_Data_AllDHW)

# Get fixed effects from original model
site_fe <- fixef(model_negative_binomial)$Site_ID
year_fe <- fixef(model_negative_binomial)$Date_Year

cat("Original model coefficients:\n")
print(coef(model_negative_binomial))
cat("\n")

# Counterfactual file paths
counterfactual_files <- paste0("/Users/pujapande/Library/CloudStorage/Dropbox/Coracle/Counterfactuals_12_12/DHW_count_source_mm_", 0:9, ".nc")

#### Helper function: Extract counterfactual DHW ####
extract_counterfactual_dhw <- function(nc_file, site_ids, year_month_indices, sample_idx, scenario_name) {
  # Open nc file
  nc <- nc_open(nc_file)
  
  # Get dimension info
  site_dim <- ncvar_get(nc, "Site_ID")
  
  # Get scenario index (firm dimension)
  scenarios <- c("USA", "Europe", "natural", "all", "topten")
  scenario_idx <- which(scenarios == scenario_name)
  
  # Extract full DHW array: [time, Site_ID, sample, firm]
  dhw_full <- ncvar_get(nc, "DHW")
  
  # Close file
  nc_close(nc)
  
  # Get valid time range
  max_time <- dim(dhw_full)[1]
  
  # Extract DHW for specific sites, times, sample, and scenario
  dhw_values <- sapply(1:length(site_ids), function(i) {
    site_idx <- which(site_dim == site_ids[i])
    time_idx <- year_month_indices[i]
    
    # Check bounds
    if(length(site_idx) == 0 || is.na(time_idx) || time_idx > max_time || time_idx < 1) {
      return(NA)
    }
    
    # DHW[time, Site_ID, sample, firm]
    dhw_full[time_idx, site_idx, sample_idx, scenario_idx]
  })
  
  return(dhw_values)
}

#### Create time index for observations ####
time_origin <- as.Date("1985-01-31")
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    obs_date = as.Date(paste(Date_Year, Date_Month, "15", sep = "-")),
    time_idx = as.numeric(floor((obs_date - time_origin) / 30.44)) + 1
  ) #%>%
  #filter(time_idx >= 1 & time_idx <= 480)

cat("Filtered observations:", nrow(All_Bleaching_Events_Data_AllDHW), "\n")
cat("Time index range:", range(All_Bleaching_Events_Data_AllDHW$time_idx, na.rm = TRUE), "\n\n")

#### Monte Carlo uncertainty propagation - ALL SCENARIOS ####
set.seed(123)
n_iterations <- 1000
scenarios <- c("natural", "USA", "Europe", "all", "topten")  # ALL scenarios

n_obs <- nrow(All_Bleaching_Events_Data_AllDHW)
results_list <- vector("list", length = n_iterations)

cat("Starting Monte Carlo uncertainty propagation...\n")
cat("Iterations:", n_iterations, "\n")
cat("Observations:", n_obs, "\n")
cat("Scenarios:", length(scenarios), "\n\n")

pb <- txtProgressBar(min = 0, max = n_iterations, style = 3)

# Monte Carlo loop with fixed effects
for(iter in 1:n_iterations) {
  
  # Step 1: Sample impact model coefficient (only beta)
  boot_idx <- sample(1:nrow(boot_results_clean), 1)
  beta <- boot_results_clean$log1p_dhw[boot_idx]
  
  # Step 2: Sample counterfactual file and uncertainty sample
  file_idx <- sample(1:10, 1)
  sample_idx <- sample(1:10, 1)
  nc_file <- counterfactual_files[file_idx]
  
  # Initialize iteration results
  iter_results <- data.frame(
    iteration = iter,
    Site_ID = All_Bleaching_Events_Data_AllDHW$Site_ID,
    Date_Year = All_Bleaching_Events_Data_AllDHW$Date_Year,
    Date_Month = All_Bleaching_Events_Data_AllDHW$Date_Month,
    Ecoregion_Name = All_Bleaching_Events_Data_AllDHW$Ecoregion_Name,
    observed_dhw = All_Bleaching_Events_Data_AllDHW$dhw,
    observed_bleaching = All_Bleaching_Events_Data_AllDHW$Percent_Bleached,
    pred_factual = NA
  )
  
  # Predict factual bleaching with fixed effects
  for(i in 1:nrow(iter_results)) {
    # Get site and year fixed effects
    site_id_char <- as.character(iter_results$Site_ID[i])
    year_char <- as.character(iter_results$Date_Year[i])
    
    # Get fixed effect values (default to 0 if missing)
    site_effect <- ifelse(site_id_char %in% names(site_fe), site_fe[site_id_char], 0)
    year_effect <- ifelse(year_char %in% names(year_fe), year_fe[year_char], 0)
    
    # Predict factual: exp(site_FE + year_FE + β × log1p(dhw))
    iter_results$pred_factual[i] <- exp(site_effect + year_effect + beta * log1p(iter_results$observed_dhw[i]))
  }
  
  # Step 3: Extract counterfactual DHW and predict for all scenarios
  for(scenario in scenarios) {
    
    # Extract counterfactual DHW
    cf_dhw <- extract_counterfactual_dhw(
      nc_file = nc_file,
      site_ids = All_Bleaching_Events_Data_AllDHW$Site_ID,
      year_month_indices = All_Bleaching_Events_Data_AllDHW$time_idx,
      sample_idx = sample_idx,
      scenario_name = scenario
    )
    
    # Predict counterfactual bleaching with fixed effects
    pred_cf <- rep(NA, length(cf_dhw))
    for(i in 1:length(cf_dhw)) {
      # Get site and year fixed effects (same as factual)
      site_id_char <- as.character(iter_results$Site_ID[i])
      year_char <- as.character(iter_results$Date_Year[i])
      
      site_effect <- ifelse(site_id_char %in% names(site_fe), site_fe[site_id_char], 0)
      year_effect <- ifelse(year_char %in% names(year_fe), year_fe[year_char], 0)
      
      # Predict counterfactual: exp(site_FE + year_FE + β × log1p(cf_dhw))
      pred_cf[i] <- exp(site_effect + year_effect + beta * log1p(cf_dhw[i]))
    }
    
    # Store results
    iter_results[[paste0("cf_dhw_", scenario)]] <- cf_dhw
    iter_results[[paste0("pred_cf_", scenario)]] <- pred_cf
  }
  
  # Store iteration results
  results_list[[iter]] <- iter_results
  
  setTxtProgressBar(pb, iter)
}

close(pb)

# Combine all iterations
all_results <- bind_rows(results_list)

cat("\n\nMonte Carlo complete!\n")
cat("Total results rows:", nrow(all_results), "\n\n")

#### Calculate summary statistics ####
summary_stats <- all_results %>%
  group_by(Site_ID, Date_Year, Date_Month, Ecoregion_Name, observed_bleaching) %>%
  summarise(
    # Observed DHW
    observed_dhw_mean = mean(observed_dhw, na.rm = TRUE),
    
    # Factual predictions
    pred_factual_mean = mean(pred_factual, na.rm = TRUE),
    pred_factual_sd = sd(pred_factual, na.rm = TRUE),
    pred_factual_ci_lower = quantile(pred_factual, 0.025, na.rm = TRUE),
    pred_factual_ci_upper = quantile(pred_factual, 0.975, na.rm = TRUE),
    
    # Counterfactual predictions for each scenario
    across(starts_with("pred_cf_"), 
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE),
                ci_lower = ~quantile(.x, 0.025, na.rm = TRUE),
                ci_upper = ~quantile(.x, 0.975, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    
    .groups = "drop"
  )

cat("Summary statistics calculated\n")
print(head(summary_stats))

#### Save results ####
#write.csv(all_results, "Attribution_Uncertainty_Full_Results.csv", row.names = FALSE)
#write.csv(summary_stats, "Attribution_Uncertainty_Summary.csv", row.names = FALSE)

cat("\nResults saved!\n\n")

#### PLOTTING ####
cat("Creating plots...\n")

#### 1. Time series with all scenarios ####
bleaching_time <- summary_stats %>%
  group_by(Date_Year) %>%
  summarise(
    Factual = mean(pred_factual_mean, na.rm = TRUE),
    Natural = mean(pred_cf_natural_mean, na.rm = TRUE),
    USA = mean(pred_cf_USA_mean, na.rm = TRUE),
    Europe = mean(pred_cf_Europe_mean, na.rm = TRUE),
    All_Majors = mean(pred_cf_all_mean, na.rm = TRUE),
    Top_Ten = mean(pred_cf_topten_mean, na.rm = TRUE),
    
    # CI for natural scenario
    Natural_lower = mean(pred_cf_natural_ci_lower, na.rm = TRUE),
    Natural_upper = mean(pred_cf_natural_ci_upper, na.rm = TRUE)
  )

# Time series plot
p1 <- ggplot(bleaching_time, aes(x = Date_Year)) +
  # Uncertainty band for natural
  #geom_ribbon(aes(ymin = Natural_lower, ymax = Natural_upper), alpha = 0.2, fill = "black") +
  
  # Lines for all scenarios
  geom_line(aes(y = Factual, color = "Factual"), linewidth = 1) +
  geom_line(aes(y = Natural, color = "Natural"), linewidth = 1) +
  geom_line(aes(y = All_Majors, color = "All Majors"), linewidth = 1) +
  geom_line(aes(y = Top_Ten, color = "Top Ten"), linewidth = 1) +
  geom_line(aes(y = USA, color = "USA"), linewidth = 1) +
  geom_line(aes(y = Europe, color = "Europe"), linewidth = 1) +
  
  scale_color_manual(values = c(
    "Factual" = "red",
    "Natural" = "black",
    "All Majors" = "blue",
    "Top Ten" = "orange",
    "USA" = "green",
    "Europe" = "purple"
  )) +
  labs(
    title = "Coral Bleaching: Factual vs Counterfactual Scenarios",
    subtitle = "Shaded area = 95% CI for natural-only forcing",
    x = "Year",
    y = "Mean Predicted Bleaching (%)",
    color = "Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(p1)

#### 2. Faceted by scenario ####
bleaching_facet <- bleaching_time %>%
  select(Date_Year, Factual, Natural, USA, Europe, All_Majors, Top_Ten) %>%
  pivot_longer(cols = -Date_Year, names_to = "Scenario", values_to = "Mean_Bleaching") %>%
  mutate(Scenario = factor(Scenario, 
                           levels = c("Factual", "Natural", "All_Majors", "Top_Ten", "USA", "Europe")))

p2 <- ggplot(bleaching_facet, aes(x = Date_Year, y = Mean_Bleaching, color = Scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  facet_wrap(~Scenario, ncol = 2) +
  scale_color_manual(values = c(
    "Factual" = "red",
    "Natural" = "black",
    "All_Majors" = "blue",
    "Top_Ten" = "orange",
    "USA" = "green",
    "Europe" = "purple"
  )) +
  labs(
    title = "Bleaching Progression by Scenario",
    x = "Year",
    y = "Mean Predicted Bleaching (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p2)

#### 3. Attribution (difference from factual) ####
bleaching_attr <- summary_stats %>%
  group_by(Date_Year) %>%
  summarise(
    Natural = mean(pred_factual_mean - pred_cf_natural_mean, na.rm = TRUE),
    USA = mean(pred_factual_mean - pred_cf_USA_mean, na.rm = TRUE),
    Europe = mean(pred_factual_mean - pred_cf_Europe_mean, na.rm = TRUE),
    All_Majors = mean(pred_factual_mean - pred_cf_all_mean, na.rm = TRUE),
    Top_Ten = mean(pred_factual_mean - pred_cf_topten_mean, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = -Date_Year, names_to = "Scenario", values_to = "Attributable")

p3 <- ggplot(bleaching_attr, aes(x = Date_Year, y = Attributable, fill = Scenario)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(
    "Natural" = "black", 
    "All_Majors" = "blue",
    "Top_Ten" = "orange",
    "USA" = "green", 
    "Europe" = "purple"
  )) +
  labs(
    title = "Bleaching Attributable to Human Forcing",
    subtitle = "Factual - Counterfactual",
    x = "Year",
    y = "Attributable Bleaching (%)"
  ) +
  theme_minimal()

print(p3)

cat("\nPlots saved!\n")
cat("- Bleaching_Timeseries_All_Scenarios.png\n")
cat("- Bleaching_Timeseries_Faceted.png\n")
cat("- Bleaching_Attribution_Barplot.png\n")