#### author: Puja Pande 
#### date: 20/05/2026
#### description: code for figure 3!!
#### note: builds Supplementary Figure S6 (attribution panels by country/company),
####       adapted from the Figure 3 pipeline with added filter support
library(ncdf4)
library(dplyr)
library(ggplot2)
library(patchwork)
library(MetBrewer)

#### Read batch metadata ####
batch_dir <- "path/to/hpc_output/linear_lat_conley_season_nd/"

batch_files <- sort(list.files(batch_dir, pattern = "\\.nc$", full.names = TRUE))

nc1 <- nc_open(batch_files[1])
site_ids_all <- as.integer(strsplit(ncatt_get(nc1, 0, "site_ids")$value, ",")[[1]])
years_all <- as.integer(strsplit(ncatt_get(nc1, 0, "years")$value, ",")[[1]])
scenarios <- strsplit(ncatt_get(nc1, 0, "scenarios")$value, ",")[[1]]
nc_close(nc1)

n_sites <- length(site_ids_all)
n_years <- length(years_all)
n_scenarios <- length(scenarios)
n_samples <- length(batch_files) * 10  # 840 total

cat("Sites:", n_sites, "Years:", n_years, "Scenarios:", n_scenarios, "Samples:", n_samples, "\n")

#### ── ADD FILTER SUPPORT TO PIPELINE ──────────────────────────────────── ####
get_batch_counts_corrected_v2 <- function(batch_dir, threshold = 0,
                                          site_filter = NULL,
                                          year_from   = NULL) {
  batch_files <- sort(list.files(batch_dir, pattern = "\\.nc$", full.names = TRUE))
  
  # Year index filter
  yr_idx <- if (!is.null(year_from)) which(years_all >= year_from) else seq_len(n_years)
  n_years_filt <- length(yr_idx)
  
  # Site index filter (rows of pf/pcf map to site_ids_all)
  site_keep <- if (!is.null(site_filter)) site_ids_all %in% site_filter else rep(TRUE, n_sites)
  
  factual_count <- matrix(NA_real_, nrow = n_years_filt, ncol = n_samples)
  cf_count      <- array(NA_real_,  dim = c(n_years_filt, n_samples, n_scenarios))
  
  for (b in seq_along(batch_files)) {
    cat("Batch", b, "/", length(batch_files), "\n")
    nc  <- nc_open(batch_files[b])
    pf  <- ncvar_get(nc, "pred_factual"); pf[pf == -9999]   <- NA
    pcf <- ncvar_get(nc, "pred_cf");      pcf[pcf == -9999] <- NA
    nc_close(nc)
    
    # Apply site filter
    pf  <- pf[site_keep, , , drop = FALSE]
    pcf <- pcf[site_keep, , , , drop = FALSE]
    
    s_idx <- ((b - 1) * 10 + 1):(b * 10)
    
    for (yi in seq_along(yr_idx)) {
      y <- yr_idx[yi]
      for (s in 1:10) {
        sample_idx <- s_idx[s]
        
        factual_exceeding   <- pf[, y, s] > threshold
        n_factual_exceeding <- sum(factual_exceeding, na.rm = TRUE)
        factual_count[yi, sample_idx] <- n_factual_exceeding
        
        if (n_factual_exceeding == 0) {
          cf_count[yi, sample_idx, ] <- 0
          next
        }
        for (sc in seq_len(n_scenarios)) {
          cf_values <- pcf[factual_exceeding, y, s, sc]
          cf_count[yi, sample_idx, sc] <- sum(cf_values > threshold, na.rm = TRUE)
        }
      }
    }
    rm(pf, pcf); gc()
  }
  
  list(factual_count = factual_count, cf_count = cf_count,
       years = years_all[yr_idx], scenario_names = scenarios, threshold = threshold)
}

#### Run ####
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

# company_labels <- c(
#   "natural" = nat_label, "Saudi Aramco" = "Saudi Aramco", "ExxonMobil" = "ExxonMobil",
#   "Chevron" = "Chevron", "Holcim Group" = "Holcim Group", "All_US" = "All US",
#   "topten" = "Top 10", "all" = "All"
# )


company_labels <- c(
  "natural" = nat_label, "Saudi Aramco" = "Saudi Aramco", "ExxonMobil" = "ExxonMobil",
  "Chevron" = "Chevron", "topten" = "Top 10", "all" = "All", "Gazprom" = "Gazprom", 
  "National Iranian Oil Company" = "National Iranian\nOil Company", 
  "BP" = "BP"
)


#### PLOTTING FUNCTIONS ####
bar_theme <- theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.y       = element_text(size = 11),
    plot.title         = element_text(size = 14, hjust = 0.5),
    plot.tag           = element_text(size = 14, face = "bold"),
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

make_bar <- function(df, tag, title, y_limits = c(0, 100)) {
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
      #y = "Percentage of site-level bleaching events (2024)", 
      y = "Percentage of site-level bleaching events", 
      x = NULL, tag = tag, title = title) +
    bar_theme
  
  p
}

#### Generate plotting data ####
#### ── PANEL SITE IDs ──────────────────────────────────────────────────── ####
# Adjust the data object name if different in your environment
All_Bleaching_Events_Data_AllDHW <- read.csv("path/to/Coracle/Datasets/Final_Panel_Data_allDHW_01July.csv")

panel_site_ids <- All_Bleaching_Events_Data_AllDHW %>%
  group_by(Site_ID) %>% filter(n() > 1) %>% pull(Site_ID) %>% unique()

#### ── RUN 4 VERSIONS (10% threshold) ─────────────────────────────────── ####
step2_all        <- get_batch_counts_corrected_v2(batch_dir, threshold = 20)
step2_panel      <- get_batch_counts_corrected_v2(batch_dir, threshold = 20, site_filter = panel_site_ids)
step2_1998       <- get_batch_counts_corrected_v2(batch_dir, threshold = 20, year_from = 1998)
step2_panel_1998 <- get_batch_counts_corrected_v2(batch_dir, threshold = 20, site_filter = panel_site_ids, year_from = 1998)

#### ── ATTRIBUTION DATAFRAMES ──────────────────────────────────────────── ####
# Countries
tot_country_all        <- get_total_attribution_df(step2_all,        country_labels)
tot_country_panel      <- get_total_attribution_df(step2_panel,      country_labels)
tot_country_1998       <- get_total_attribution_df(step2_1998,       country_labels)
tot_country_panel_1998 <- get_total_attribution_df(step2_panel_1998, country_labels)

# Companies
tot_company_all        <- get_total_attribution_df(step2_all,        company_labels)
tot_company_panel      <- get_total_attribution_df(step2_panel,      company_labels)
tot_company_1998       <- get_total_attribution_df(step2_1998,       company_labels)
tot_company_panel_1998 <- get_total_attribution_df(step2_panel_1998, company_labels)

#### ── COUNTRY FIGURE (2×2) ────────────────────────────────────────────── ####
p_country_all        <- make_bar(tot_country_all,        "a", "Country contributions to bleaching (>20%), all sites 1985–2024")
p_country_panel      <- make_bar(tot_country_panel,      "b", "Country contributions to bleaching (>20%), panel sites 1985–2024")
p_country_1998       <- make_bar(tot_country_1998,       "c", "Country contributions to bleaching (>20%), all sites 1998–2024")
p_country_panel_1998 <- make_bar(tot_country_panel_1998, "d", "Country contributions to bleaching (>20%), panel sites 1998–2024")

fig_country <- (p_country_all | p_country_panel) /
  (p_country_1998 | p_country_panel_1998)
fig_country

#### ── COMPANY FIGURE (2×2) ────────────────────────────────────────────── ####
p_company_all        <- make_bar(tot_company_all,        "a", "All sites (1985–2024)")
p_company_panel      <- make_bar(tot_company_panel,      "b", "Panel sites (1985–2024)")
p_company_1998       <- make_bar(tot_company_1998,       "c", "All sites (1998–2024)")
p_company_panel_1998 <- make_bar(tot_company_panel_1998, "d", "Panel sites (1998–2024)")

fig_company <- (p_company_all | p_company_panel) /
  (p_company_1998 | p_company_panel_1998)
fig_company

#### ── SAVE ────────────────────────────────────────────────────────────── ####
ggsave("fig_country_10pct_sensitivity.pdf", fig_country,
       width = 250, height = 200, units = "mm", dpi = 300)
ggsave("fig_company_10pct_sensitivity.pdf", fig_company,
       width = 250, height = 200, units = "mm", dpi = 300)