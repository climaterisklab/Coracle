#### author: Puja Pande
#### date: 03 November 2025
#### description: block bootstrapping for uncertainty estimation

#### libraries ####
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)

#### mildly obsolete but bootstrapping ####
## may need to use some data from MarineHeatwaveRegressions - GMST data and all_ts_data

# Storage for all bootstrap results
all_bootstrap_results <- list()

# Run all months sequentially
for (target_month in 1:12) {
  cat("STARTING MONTH:", month.name[target_month], "\n")
  
  set.seed(123)
  
  # Parameters
  block_width <- 2
  n_bootstrap <- 100
  
  # Setup parallel backend
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  cat("Processing Month:", month.name[target_month], "\n")
  cat("Using", n_cores, "cores for parallel processing\n")
  
  # Get the range of years
  year_range <- range(all_ts_data$year)
  all_years <- year_range[1]:year_range[2]
  n_years <- length(all_years)
  
  # Create blocks
  blocks <- list()
  for (i in seq(1, n_years - block_width + 1, by = 1)) {
    blocks[[length(blocks) + 1]] <- all_years[i:(i + block_width - 1)]
  }
  
  # Pre-filter data for target month (more efficient)
  month_data <- all_ts_data %>% filter(month == target_month)
  
  cat("Starting", n_bootstrap, "bootstrap iterations...\n")
  
  start_time <- Sys.time()
  
  # Parallel bootstrap
  bootstrap_results <- foreach(
    boot_iter = 1:n_bootstrap,
    .combine = rbind,
    .packages = c('dplyr', 'tidyr'),
    .errorhandling = 'pass'
  ) %dopar% {
    
    # Sample blocks with replacement
    set.seed(123 + target_month * 1000 + boot_iter)  # Reproducible seed
    n_blocks_needed <- ceiling(n_years / block_width)
    sampled_blocks <- sample(blocks, size = n_blocks_needed, replace = TRUE)
    
    # Flatten to get bootstrap sample of years
    boot_years <- unlist(sampled_blocks)[1:n_years]
    
    # Create mapping from original years to bootstrap years
    year_mapping <- data.frame(
      original_year = all_years,
      boot_year = boot_years
    )
    
    # Bootstrap GMST data
    boot_gmst <- gmst_data %>%
      inner_join(year_mapping, by = c("Year" = "boot_year")) %>%
      mutate(Year = original_year) %>%
      select(Year, GMT, date)
    
    # Storage for this bootstrap iteration
    iter_results <- list()
    
    # Loop over each location
    for (loc_id in seq_len(nrow(unique_locations))) {
      
      # Get this location's data (already filtered to target month)
      loc_data <- month_data %>% filter(location_id == loc_id)
      
      if (nrow(loc_data) == 0) next
      
      lon <- loc_data$lon[1]
      lat <- loc_data$lat[1]
      
      # Bootstrap DHW data using the same year mapping
      boot_dhw <- loc_data %>%
        inner_join(year_mapping, by = c("year" = "boot_year")) %>%
        mutate(year = original_year)
      
      # Aggregate to yearly
      ts_yearly <- boot_dhw %>%
        group_by(year) %>%
        summarise(mean_dhw = mean(dhw, na.rm = TRUE), .groups = 'drop')
      
      # Join with bootstrapped GMST
      merged_df <- ts_yearly %>%
        inner_join(boot_gmst, by = c("year" = "Year"))
      
      if (nrow(merged_df) < 2) next
      
      # Remove NAs
      merged_df <- merged_df %>% drop_na(mean_dhw, GMT)
      
      if (nrow(merged_df) < 2) next
      
      # Handle zeros
      merged_df$mean_dhw <- pmax(merged_df$mean_dhw, 0.001)
      
      # Fit Poisson regression
      tryCatch({
        lm_fit <- glm(mean_dhw ~ GMT, family = "poisson", data = merged_df)
        coef_summary <- summary(lm_fit)$coefficients
        
        # Store only essential results for bootstrap
        boot_result <- data.frame(
          boot_iteration = boot_iter,
          location_id = loc_id,
          longitude = lon,
          latitude = lat,
          month = target_month,
          beta = coef_summary[2, 1],
          beta_se = coef_summary[2, 2],
          p_value = coef_summary[2, 4]
        )
        
        iter_results[[length(iter_results) + 1]] <- boot_result
        
      }, error = function(e) {
        # Silently skip failed models
      })
    }
    
    # Return combined results for this iteration
    if (length(iter_results) > 0) {
      bind_rows(iter_results)
    } else {
      NULL
    }
  }
  
  # Stop cluster
  stopCluster(cl)
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat("Bootstrap complete for", month.name[target_month], "\n")
  cat("Time elapsed:", round(elapsed_time, 2), "minutes\n")
  
  # Save bootstrap results
  output_file <- paste0("bootstrap_results_month_", target_month, ".csv")
  write.csv(bootstrap_results, output_file, row.names = FALSE)
  
  cat("Iterations completed:", n_distinct(bootstrap_results$boot_iteration), "\n")
  
  # Store in list for later use
  all_bootstrap_results[[target_month]] <- bootstrap_results
  
  # Also create a named object in the environment for easy access
  assign(paste0("bootstrap_month_", target_month), bootstrap_results, envir = .GlobalEnv)
  
  cat("COMPLETED MONTH:", month.name[target_month], "\n\n")
}

cat("\n====================================\n")
cat("ALL MONTHS COMPLETE!\n")
cat("====================================\n")

# Optionally combine all months into one dataframe
bootstrap_all_months <- bind_rows(all_bootstrap_results)

# Save combined results
write.csv(bootstrap_all_months, "bootstrap_results_all_months_combined.csv", row.names = FALSE)
cat("\nCombined results saved to: bootstrap_results_all_months_combined.csv\n")




#### plots??? ####

#### Step 8: Quantify Uncertainty from Bootstrap Samples - ALL MONTHS ####

cat("\n====================================\n")
cat("STEP 8: CALCULATING UNCERTAINTY\n")
cat("====================================\n\n")

# Storage for uncertainty results
all_uncertainty_results <- list()

# Process each month
for (target_month in 1:12) {
  
  cat("Processing uncertainty for", month.name[target_month], "...\n")
  
  # Load bootstrap results for this month
  bootstrap_df <- read.csv(paste0("bootstrap_results_month_", target_month, ".csv"))
  
  # Calculate uncertainty metrics for each location
  uncertainty_summary <- bootstrap_df %>%
    group_by(location_id, longitude, latitude, month) %>%
    summarise(
      # Number of successful bootstrap iterations
      n_boot = n(),
      
      # Beta coefficient statistics
      beta_mean = mean(beta, na.rm = TRUE),
      beta_median = median(beta, na.rm = TRUE),
      beta_sd = sd(beta, na.rm = TRUE),
      
      # Confidence intervals (95%)
      beta_ci_lower_95 = quantile(beta, probs = 0.025, na.rm = TRUE),
      beta_ci_upper_95 = quantile(beta, probs = 0.975, na.rm = TRUE),
      
      # Confidence intervals (90%)
      beta_ci_lower_90 = quantile(beta, probs = 0.05, na.rm = TRUE),
      beta_ci_upper_90 = quantile(beta, probs = 0.95, na.rm = TRUE),
      
      # Interquartile range
      beta_q25 = quantile(beta, probs = 0.25, na.rm = TRUE),
      beta_q75 = quantile(beta, probs = 0.75, na.rm = TRUE),
      beta_iqr = IQR(beta, na.rm = TRUE),
      
      # Coefficient of variation (relative uncertainty)
      beta_cv = beta_sd / abs(beta_mean),
      
      # Width of confidence interval
      ci_width_95 = beta_ci_upper_95 - beta_ci_lower_95,
      ci_width_90 = beta_ci_upper_90 - beta_ci_lower_90,
      
      # Statistical significance (bootstrap p-value)
      # Proportion of bootstrap samples where beta crosses zero
      prop_positive = mean(beta > 0, na.rm = TRUE),
      prop_negative = mean(beta < 0, na.rm = TRUE),
      boot_pvalue = 2 * pmin(prop_positive, prop_negative),
      
      # Mean standard error from individual models
      mean_beta_se = mean(beta_se, na.rm = TRUE),
      
      .groups = 'drop'
    )
  
  # Add significance flag based on bootstrap CI
  uncertainty_summary <- uncertainty_summary %>%
    mutate(
      significant_95 = !(beta_ci_lower_95 <= 0 & beta_ci_upper_95 >= 0),
      significant_90 = !(beta_ci_lower_90 <= 0 & beta_ci_upper_90 >= 0),
      sign_beta = sign(beta_mean)
    )
  
  # Merge with original model results
  original_results <- final_poisson_results %>%
    filter(month == target_month) %>%
    select(location_id, longitude, latitude, month, 
           beta_original = beta, 
           beta_se_original = beta_se,
           p_value_original = p_value,
           pseudo_r2, dispersion, n_obs)
  
  uncertainty_with_original <- left_join(
    uncertainty_summary,
    original_results,
    by = c("location_id", "longitude", "latitude", "month")
  )
  
  # Calculate bias (difference between bootstrap mean and original estimate)
  uncertainty_with_original <- uncertainty_with_original %>%
    mutate(
      bootstrap_bias = beta_mean - beta_original,
      bias_percent = (bootstrap_bias / beta_original) * 100
    )
  
  # Save uncertainty results for this month
  output_file <- paste0("bootstrap_uncertainty_month_", target_month, ".csv")
  write.csv(uncertainty_with_original, output_file, row.names = FALSE)
  
  # Store in list
  all_uncertainty_results[[target_month]] <- uncertainty_with_original
  
  # Also create named object in environment
  assign(paste0("uncertainty_month_", target_month), uncertainty_with_original, envir = .GlobalEnv)
  
  # Print summary statistics for this month
  cat("  Locations analyzed:", nrow(uncertainty_with_original), "\n")
  cat("  Mean beta:", round(mean(uncertainty_with_original$beta_mean, na.rm = TRUE), 4), "\n")
  cat("  Median CI width (95%):", round(median(uncertainty_with_original$ci_width_95, na.rm = TRUE), 4), "\n")
  cat("  Significant locations (95% CI):", 
      sum(uncertainty_with_original$significant_95, na.rm = TRUE), 
      "(", round(100 * mean(uncertainty_with_original$significant_95, na.rm = TRUE), 1), "%)\n")
  cat("  Saved to:", output_file, "\n\n")
}

cat("====================================\n")
cat("UNCERTAINTY QUANTIFICATION COMPLETE\n")
cat("====================================\n\n")

# Combine all months
uncertainty_all_months <- bind_rows(all_uncertainty_results)

# Save combined uncertainty results
write.csv(uncertainty_all_months, "bootstrap_uncertainty_all_months.csv", row.names = FALSE)

cat("Combined uncertainty results saved to: bootstrap_uncertainty_all_months.csv\n")
cat("Total location-month combinations:", nrow(uncertainty_all_months), "\n\n")

cat("Uncertainty results are available as:\n")
cat("  - Individual objects: uncertainty_month_1, uncertainty_month_2, ..., uncertainty_month_12\n")
cat("  - Combined dataframe: uncertainty_all_months\n")
cat("  - List: all_uncertainty_results\n\n")

#### Summary Statistics Across All Months ####
cat("=== Summary Statistics Across All Months ===\n")

summary_stats <- uncertainty_all_months %>%
  group_by(month) %>%
  summarise(
    n_locations = n(),
    mean_beta = mean(beta_mean, na.rm = TRUE),
    median_ci_width = median(ci_width_95, na.rm = TRUE),
    pct_significant = 100 * mean(significant_95, na.rm = TRUE),
    mean_cv = mean(beta_cv, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(month_name = month.name[month])

print(summary_stats)

write.csv(summary_stats, "uncertainty_summary_by_month.csv", row.names = FALSE)

#### Visualizations ####
cat("\n=== Creating Visualizations ===\n")

# 1. Spatial map of uncertainty for each month
library(ggplot2)

for (m in 1:12) {
  month_data <- uncertainty_all_months %>% filter(month == m)
  
  p <- ggplot(month_data, aes(x = longitude, y = latitude, color = ci_width_95)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_viridis_c(option = "plasma") +
    labs(title = paste("Bootstrap Uncertainty -", month.name[m]),
         subtitle = "Width of 95% Confidence Interval for Beta",
         x = "Longitude",
         y = "Latitude",
         color = "CI Width") +
    theme_minimal() +
    coord_fixed()
  print(p)
  #ggsave(paste0("spatial_uncertainty_month_", m, ".png"), 
  #p, width = 12, height = 8)
}


p_all <- ggplot(uncertainty_all_months, 
                aes(x = longitude, y = latitude, color = ci_width_95)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "Bootstrap Uncertainty - All Months",
       subtitle = "Width of 95% Confidence Interval for Beta",
       x = "Longitude",
       y = "Latitude",
       color = "CI Width") +
  theme_minimal() +
  coord_fixed() +
  facet_wrap(~month, ncol = 4, labeller = labeller(month = month.name))

print(p_all)


cat("Spatial maps saved for all 12 months\n")

# 2. Beta coefficient vs uncertainty across all months
p_beta_vs_uncertainty <- ggplot(uncertainty_all_months, 
                                aes(x = beta_mean, y = beta_sd, color = as.factor(month))) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_viridis_d(option = "turbo", name = "Month") +
  labs(title = "Beta Coefficient vs Bootstrap Uncertainty - All Months",
       x = "Beta (Bootstrap Mean)",
       y = "Beta Standard Deviation") +
  theme_minimal() +
  facet_wrap(~month, ncol = 4, labeller = labeller(month = month.name))

#ggsave("beta_vs_uncertainty_all_months.png", 
p_beta_vs_uncertainty, width = 16, height = 12)

cat("Faceted plot saved: beta_vs_uncertainty_all_months.png\n")

# 3. Comparison with original estimates
p_comparison <- ggplot(uncertainty_all_months, 
                       aes(x = beta_original, y = beta_mean, color = as.factor(month))) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_color_viridis_d(option = "turbo", name = "Month") +
  labs(title = "Bootstrap vs Original Beta Estimates - All Months",
       x = "Original Beta Estimate",
       y = "Bootstrap Mean Beta") +
  theme_minimal()

#ggsave("bootstrap_vs_original_all_months.png", 
p_comparison, width = 12, height = 8)

cat("Comparison plot saved: bootstrap_vs_original_all_months.png\n")

# 4. Distribution of CI widths by month
p_ci_dist <- ggplot(uncertainty_all_months, aes(x = as.factor(month), y = ci_width_95)) +
  geom_boxplot(aes(fill = as.factor(month)), alpha = 0.7) +
  scale_fill_viridis_d(option = "turbo", guide = "none") +
  scale_x_discrete(labels = month.abb) +
  labs(title = "Distribution of Uncertainty (95% CI Width) by Month",
       x = "Month",
       y = "95% Confidence Interval Width") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#ggsave("ci_width_by_month.png", p_ci_dist, width = 12, height = 6)

cat("CI width distribution saved: ci_width_by_month.png\n")

# 5. Significance comparison: original vs bootstrap
significance_comparison <- uncertainty_all_months %>%
  mutate(
    original_sig = p_value_original < 0.05,
    bootstrap_sig = significant_95
  ) %>%
  group_by(month) %>%
  summarise(
    n_total = n(),
    original_sig_count = sum(original_sig, na.rm = TRUE),
    bootstrap_sig_count = sum(bootstrap_sig, na.rm = TRUE),
    both_sig = sum(original_sig & bootstrap_sig, na.rm = TRUE),
    neither_sig = sum(!original_sig & !bootstrap_sig, na.rm = TRUE),
    only_original = sum(original_sig & !bootstrap_sig, na.rm = TRUE),
    only_bootstrap = sum(!original_sig & bootstrap_sig, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(month_name = month.name[month])

print(significance_comparison)
write.csv(significance_comparison, "significance_comparison.csv", row.names = FALSE)

cat("\n=== Step 8 Complete ===\n")
cat("All uncertainty metrics calculated and saved!\n")
