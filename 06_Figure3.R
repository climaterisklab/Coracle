#### author: Puja Pande 
#### date: 20/05/2026
#### description: code for figure 3!!
#### note: builds main text Figure 3 (coral bleaching attribution ridge plots and bar
####       charts from FaIR counterfactual simulations)
library(ncdf4)
library(dplyr)
library(ggplot2)
library(patchwork)
library(MetBrewer)

#### Read counts ####
# v2_mixed impact model (province-specific linear/sigmoid), all reef locations (~72,500 5 km cells).
# The 106 GB of per-cell output on BSC (/esarchive/Coracle/bleaching_attribution_v2/reef_locations/) is
# reduced on the cluster to per-country event counts by 04_HPC_Scripts/reef_counts_v2.py
# (slurm_reef_counts_v2.sh); this script only reads that small file.
source("reef_counts_v2_helpers.R")
counts_file <- Sys.getenv("REEF_COUNTS", "path/to/reef_counts_v2.nc")
fig_dir     <- Sys.getenv("FIG_DIR", ".")
# CELLSET = "reefs" (all reef cells) or "sites" (only cells with survey sites); goes into the file names
cellset     <- Sys.getenv("CELLSET", "reefs")

rc <- read_reef_counts(counts_file, cellset)
cat("Countries:", length(rc$country), "Years:", length(rc$year), "Scenarios:", length(rc$scenario),
    "Samples:", length(rc$sample), "\n")

# Global counts (all reef cells), conditional counting as before
step2_10_v2 <- group_counts(rc, 10)
step2_30_v2 <- group_counts(rc, 30)
cat("Reef cells:", step2_10_v2$n_cells, "\n")

#### Extract CUMULATIVE attribution (1985-2024) ####
get_total_attribution_df <- function(step2, label_map) {
  lapply(names(label_map), function(scen) {
    sc_idx <- which(step2$scenario_names == scen)
    
    # Sum factual and CF counts across all years
    total_factual <- colSums(step2$factual_count, na.rm = TRUE)  # 840 values
    total_cf <- colSums(step2$cf_count[, , sc_idx], na.rm = TRUE)  # 840 values
    
    # Attribution: (factual - cf) / factual * 100
    value <- (total_factual - total_cf) / total_factual * 100
    
    data.frame(
      scenario = label_map[[scen]], 
      value = value
    )
  }) %>%
    bind_rows() %>%
    group_by(scenario) %>%
    mutate(med = median(value, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(scenario = reorder(scenario, med))
}

#### Extract 2024 ONLY attribution ####
get_2024_attribution_df <- function(step2, label_map) {
  yr_idx <- 40  # 2024 is year 40
  
  lapply(names(label_map), function(scen) {
    sc_idx <- which(step2$scenario_names == scen)
    
    # Get 2024 counts only
    factual_2024 <- step2$factual_count[yr_idx, ]  # 840 values
    cf_2024 <- step2$cf_count[yr_idx, , sc_idx]     # 840 values
    
    # Attribution: (factual - cf) / factual * 100
    value <- (factual_2024 - cf_2024) / factual_2024 * 100
    
    data.frame(
      scenario = label_map[[scen]], 
      value = value
    )
  }) %>%
    bind_rows() %>%
    group_by(scenario) %>%
    mutate(med = median(value, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(scenario = reorder(scenario, med))
}

#### LABELS ####
nat_label <- "All anthropogenic\nemissions"

country_labels <- c(
  "natural" = nat_label, "AOSIS" = "AOSIS", "Indonesia" = "Indonesia",
  "India" = "India", "Japan" = "Japan", "China" = "China", "Russia" = "Russia",
  "EU27UK" = "EU27+UK", "USA" = "USA", "Australia" = "Australia"
)


company_labels <- c(
  "natural" = nat_label, "Saudi Aramco" = "Saudi Aramco", "ExxonMobil" = "ExxonMobil",
  "Chevron" = "Chevron", "topten" = "Top 10 carbon majors", "all" = "All carbon majors", "Gazprom" = "Gazprom", 
  "National Iranian Oil Company" = "National Iranian\nOil Company", 
  "BP" = "BP"
)

#### PLOTTING FUNCTIONS ####
bar_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y        = element_text(size = 6),
    axis.title.y       = element_text(size = 7),
    plot.title         = element_text(size = 7, hjust = 0.5),
    plot.tag           = element_text(size = 8, face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    legend.position    = "none"
  )

make_colours <- function(df) {
  scenarios <- levels(df$scenario)
  n <- length(scenarios)
  cols <- met.brewer("Tam", n = n + 3, type = "continuous", direction = -1)
  cols <- cols[1:n]
  names(cols) <- scenarios
  cols
}

make_bar <- function(df, tag, title, y_limits = c(0, 100), thr = 10) {
  summary_df <- df %>%
    group_by(scenario) %>%
    summarise(
      med   = median(value, na.rm = TRUE),
      ci_lo = quantile(value, 0.025, na.rm = TRUE),
      ci_hi = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  cols <- make_colours(df)
  
  p <- ggplot(summary_df, aes(x = scenario, y = med, fill = scenario)) +
    geom_col(colour = "black", linewidth = 0.3, alpha = 0.85, width = 0.7) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  width = 0.25, linewidth = 0.5, colour = "black") +
    #geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "black") +
    scale_fill_manual(values = cols) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    coord_cartesian(ylim = y_limits) +
    labs(#y = "Percentage of site-level bleaching events (1985-2024)", 
      y = paste0("Percentage of bleaching events greater than ", thr, "%"), 
      #y = "Percentage of site-level bleaching events", 
      x = NULL, tag = tag, title = title) +
    bar_theme
  
  p
}

#### Generate plotting data ####
# Cumulative (1985-2024)
#tot_attr_0_country  <- get_total_attribution_df(step2_0_v2,  country_labels)
#tot_attr_0_company  <- get_total_attribution_df(step2_0_v2,  company_labels)
tot_attr_10_country  <- get_total_attribution_df(step2_10_v2,  country_labels)
tot_attr_10_company  <- get_total_attribution_df(step2_10_v2,  company_labels)
tot_attr_30_country  <- get_total_attribution_df(step2_30_v2,  country_labels)
tot_attr_30_company  <- get_total_attribution_df(step2_30_v2,  company_labels)

# 2024 only
#yr2024_0_country  <- get_2024_attribution_df(step2_0_v2,  country_labels)
#yr2024_0_company  <- get_2024_attribution_df(step2_0_v2,  company_labels)
yr2024_10_country  <- get_2024_attribution_df(step2_10_v2,  country_labels)
yr2024_10_company  <- get_2024_attribution_df(step2_10_v2,  company_labels)
yr2024_30_country  <- get_2024_attribution_df(step2_30_v2,  country_labels)
yr2024_30_company  <- get_2024_attribution_df(step2_30_v2,  company_labels)

#### Build bar charts ####
# 2024 plots
#p_2024_0_country  <- make_bar(yr2024_0_country,  "a", "Country contributions to any bleaching (>0%)", c(0, 100))
#p_2024_0_company  <- make_bar(yr2024_0_company,  "b", "Company contributions to any bleaching (>0%)", c(0, 100))
p_2024_10_country  <- make_bar(yr2024_10_country,  "b", "Country contributions to bleaching (2024)", c(0, 100), thr = 10)
p_2024_10_company  <- make_bar(yr2024_10_company,  "d", "Carbon major contributions to bleaching (2024)", c(0, 100), thr = 10)
p_2024_30_country  <- make_bar(yr2024_30_country,  "b", "Country contributions to bleaching (2024)", c(0, 100), thr = 30)
p_2024_30_company  <- make_bar(yr2024_30_company,  "d", "Carbon major contributions to bleaching (2024)", c(0, 100), thr = 30)

# Cumulative plots
#p_cum_0_country  <- make_bar(tot_attr_0_country,  "a", "Country contributions to any bleaching (>0%)", c(0, 100))
#p_cum_0_company  <- make_bar(tot_attr_0_company,  "b", "Company contributions to any bleaching (>0%)", c(0, 100))
p_cum_10_country  <- make_bar(tot_attr_10_country,  "a", "Country contributions to bleaching (1985-2024)", c(0, 100), thr = 10)
p_cum_10_company  <- make_bar(tot_attr_10_company,  "c", "Carbon major contributions to bleaching (1985-2024)", c(0, 100), thr = 10)
p_cum_30_country  <- make_bar(tot_attr_30_country,  "a", "Country contributions to bleaching (1985-2024)", c(0, 100), thr = 30)
p_cum_30_company  <- make_bar(tot_attr_30_company,  "c", "Carbon major contributions to bleaching (1985-2024)", c(0, 100), thr = 30)

#### Assemble figures ####

fig_cumulative <- (p_cum_10_country | p_cum_30_country) / (p_cum_10_company | p_cum_30_company)
fig_cumulative



fig_2024 <- (p_2024_10_country | p_2024_30_country) / (p_2024_10_company | p_2024_30_company)
fig_2024


fig_10 <- (p_cum_10_country | p_2024_10_country) / (p_cum_10_company | p_2024_10_company)

fig_30 <- (p_cum_30_country | p_2024_30_country) / (p_cum_30_company | p_2024_30_company)



ggsave(file.path(fig_dir, paste0("fig_10_", cellset, ".png")), fig_10,
       width  = 183,
       height = 150,
       units  = "mm",
       dpi    = 600)

ggsave(file.path(fig_dir, paste0("fig_30_", cellset, ".png")), fig_30,
       width  = 183,
       height = 150,
       units  = "mm",
       dpi    = 600)


