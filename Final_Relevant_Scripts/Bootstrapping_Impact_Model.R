#### author: Puja Pande
#### date: 10 November 2025
#### description: bootstrapping coefficient for impact model

#### libraries ####
library(fixest)
library(dplyr)
library(tidyr)

#### model ####
All_Bleaching_Events_Data_AllDHW <- read.csv("Coracle_backup/Final_Relevant_Scripts/Merged_Mermaid_Panel_Bleaching_Data.csv")
# model_negative_binomial <- fenegbin(Percent_Bleached ~ log1p(dhw) | Site_ID + Date_Year, cluster = ~Ecoregion_Name, 
#                                     data = All_Bleaching_Events_Data_AllDHW)

#### Change to logit model ####
# Need Proportion_Bleached (0-1) for logit - create if not already there
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(Proportion_Bleached = Percent_Bleached / 100)

# Logit model
model_logit <- feglm(Proportion_Bleached ~ log1p(dhw) | Site_ID + Date_Year,
                     family = binomial(link = "logit"),
                     cluster = ~Ecoregion_Name,
                     data = All_Bleaching_Events_Data_AllDHW)

summary(model_logit)

#### spatial bootstrap by ecoregion ####
# Set seed for reproducibility
set.seed(123)

# Number of bootstrap iterations
n_boot <- 1000

# Get unique ecoregions
ecoregions <- unique(All_Bleaching_Events_Data_AllDHW$Ecoregion_Name)
n_ecoregions <- length(ecoregions)

# Initialize storage for bootstrap coefficients
boot_coefs <- matrix(NA, nrow = n_boot, ncol = 2)
colnames(boot_coefs) <- c("Intercept", "log1p_dhw")

# Bootstrap loop
cat("Starting spatial bootstrap by ecoregion...\n")
pb <- txtProgressBar(min = 0, max = n_boot, style = 3)

for(i in 1:n_boot) {
  
  # Resample ecoregions with replacement
  sampled_ecoregions <- sample(ecoregions, size = n_ecoregions, replace = TRUE)
  
  # Create bootstrap sample by combining all observations from sampled ecoregions
  boot_data <- NULL
  for(j in 1:length(sampled_ecoregions)) {
    eco_data <- All_Bleaching_Events_Data_AllDHW %>%
      filter(Ecoregion_Name == sampled_ecoregions[j]) %>%
      mutate(boot_ecoregion_id = j)  # Track which bootstrap sample this came from
    boot_data <- rbind(boot_data, eco_data)
  }
  boot_data <- boot_data %>%
    mutate(Proportion_Bleached = Percent_Bleached / 100)
  
  # Fit model on bootstrap sample
  tryCatch({
    boot_model <- feglm(
      Proportion_Bleached ~ log1p(dhw) | Site_ID + Date_Year,
      family = binomial(link = "logit"),
      cluster = ~boot_ecoregion_id,
      data = boot_data
    )
    # Store coefficients
    boot_coefs[i, ] <- coef(boot_model)
    
  }, error = function(e) {
    # If model fails to converge, store NA
    boot_coefs[i, ] <- NA
  })
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Convert to data frame
boot_results <- as.data.frame(boot_coefs)

# Remove any failed iterations
boot_results_clean <- boot_results %>%
  filter(!is.na(Intercept) & !is.na(log1p_dhw))

cat("\n\nBootstrap complete!\n")
cat("Successful iterations:", nrow(boot_results_clean), "out of", n_boot, "\n\n")

# Calculate bootstrap statistics
boot_summary <- data.frame(
  Coefficient = c("Intercept", "log1p(dhw)"),
  Original = coef(model_logit),
  Boot_Mean = colMeans(boot_results_clean),
  Boot_SE = apply(boot_results_clean, 2, sd),
  Boot_CI_Lower = apply(boot_results_clean, 2, quantile, probs = 0.025),
  Boot_CI_Upper = apply(boot_results_clean, 2, quantile, probs = 0.975)
)

print(boot_summary)

# Plot bootstrap distributions
par(mfrow = c(1, 2))

hist(boot_results_clean$Intercept, 
     main = "Bootstrap Distribution: Intercept",
     xlab = "Coefficient Value", 
     col = "lightblue", 
     breaks = 30)
abline(v = coef(model_logit)[1], col = "red", lwd = 2, lty = 2)
abline(v = boot_summary$Boot_CI_Lower[1], col = "darkblue", lwd = 2, lty = 3)
abline(v = boot_summary$Boot_CI_Upper[1], col = "darkblue", lwd = 2, lty = 3)

hist(boot_results_clean$log1p_dhw, 
     main = "Bootstrap Distribution: log1p(dhw)",
     xlab = "Coefficient Value", 
     col = "lightgreen", 
     breaks = 30)
abline(v = coef(model_negative_binomial)[2], col = "red", lwd = 2, lty = 2)
abline(v = boot_summary$Boot_CI_Lower[2], col = "darkblue", lwd = 2, lty = 3)
abline(v = boot_summary$Boot_CI_Upper[2], col = "darkblue", lwd = 2, lty = 3)

legend("topright", 
       legend = c("Original", "95% CI"), 
       col = c("red", "darkblue"), 
       lty = c(2, 3), 
       lwd = 2)

par(mfrow = c(1, 1))

## # Save bootstrap results
write.csv(boot_results_clean, "Bootstrap_Coefficients_Impact_Model_Logit.csv", row.names = FALSE)



