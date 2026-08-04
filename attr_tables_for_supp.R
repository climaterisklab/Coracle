
#### attribution tables for supp

library(ncdf4)
library(dplyr)

#### Load site-country mapping ####
site_country <- read.csv("/Users/pujapande/Library/CloudStorage/Dropbox/Coracle/Datasets - 01 July 2026/Final_Panel_Data_allDHW_01July.csv")

site_country_map <- site_country %>%
  select(Site_ID, Country_Name) %>%
  distinct()

#### Define AOSIS members ####
aosis_members <- c(
  'Bahamas', 'Barbados', 'Belize', 'Cuba', 'Dominica', 'Dominican Republic',
  'Fiji', 'Grenada', 'Haiti', 'Jamaica', 'Maldives', 'Marshall Islands',
  'Mauritius', 'Federated States of Micronesia', 'Palau', 'Papua New Guinea',
  'Saint Kitts and Nevis', 'Saint Lucia', 'Saint Vincent and the Grenadines',
  'Sao Tome & Principe', 'Solomon Islands', 'Trinidad and Tobago', 'Vanuatu',
  'Cook Islands'
)

#### Define EU27+UK members ####
eu27_uk_members <- c(
  'Austria', 'Belgium', 'Bulgaria', 'Croatia', 'Cyprus', 'Czech Republic',
  'Denmark', 'Estonia', 'Finland', 'France', 'Germany', 'Greece', 'Hungary',
  'Ireland', 'Italy', 'Latvia', 'Lithuania', 'Luxembourg', 'Malta',
  'Netherlands', 'Poland', 'Portugal', 'Romania', 'Slovakia', 'Slovenia',
  'Spain', 'Sweden', 'United Kingdom'
)

#### Batch metadata ####
#batch_dir <- "/Users/pujapande/Downloads/linear_lat_v2_nd/"
batch_dir <- "/Users/pujapande/Downloads/linear_lat_conley_season_nd/"
batch_files <- sort(list.files(batch_dir, pattern = "\\.nc$", full.names = TRUE))

nc1 <- nc_open(batch_files[1])
site_ids_all <- as.integer(strsplit(ncatt_get(nc1, 0, "site_ids")$value, ",")[[1]])
years_all <- as.integer(strsplit(ncatt_get(nc1, 0, "years")$value, ",")[[1]])
scenarios <- strsplit(ncatt_get(nc1, 0, "scenarios")$value, ",")[[1]]
nc_close(nc1)

n_sites <- length(site_ids_all)
n_years <- length(years_all)
n_scenarios <- length(scenarios)
n_samples <- length(batch_files) * 10

cat("Sites:", n_sites, "Years:", n_years, "Scenarios:", n_scenarios, "Samples:", n_samples, "\n")

#### Get country counts with CONDITIONAL FILTERING ####
get_batch_counts_country <- function(batch_dir, country_site_ids, site_ids_all, threshold = 0) {
  batch_files <- sort(list.files(batch_dir, pattern = "\\.nc$", full.names = TRUE))
  
  country_site_indices <- which(site_ids_all %in% country_site_ids)
  
  if (length(country_site_indices) == 0) {
    return(NULL)
  }
  
  n_years <- 40
  n_samples <- length(batch_files) * 10
  
  nc1 <- nc_open(batch_files[1])
  scenarios <- strsplit(ncatt_get(nc1, 0, "scenarios")$value, ",")[[1]]
  nc_close(nc1)
  n_scenarios <- length(scenarios)
  
  # Store COUNTS (conditional filtering)
  factual_count <- matrix(NA_real_, nrow = n_years, ncol = n_samples)
  cf_count <- array(NA_real_, dim = c(n_years, n_samples, n_scenarios))
  
  for (b in seq_along(batch_files)) {
    cat("  Batch", b, "/", length(batch_files), "\r")
    nc <- nc_open(batch_files[b])
    
    pf <- ncvar_get(nc, "pred_factual")
    pf[pf == -9999] <- NA
    
    pcf <- ncvar_get(nc, "pred_cf")
    pcf[pcf == -9999] <- NA
    
    nc_close(nc)
    
    s_idx <- ((b - 1) * 10 + 1):(b * 10)
    
    # Filter to country sites
    pf_country <- pf[country_site_indices, , , drop = FALSE]
    pcf_country <- pcf[country_site_indices, , , , drop = FALSE]
    
    # CONDITIONAL FILTERING
    for (y in seq_len(n_years)) {
      for (s in 1:10) {
        sample_idx <- s_idx[s]
        
        factual_exceeding <- pf_country[, y, s] > threshold
        n_factual_exceeding <- sum(factual_exceeding, na.rm = TRUE)
        
        factual_count[y, sample_idx] <- n_factual_exceeding
        
        if (n_factual_exceeding == 0) {
          cf_count[y, sample_idx, ] <- 0
          next
        }
        
        for (sc in seq_len(n_scenarios)) {
          cf_values <- pcf_country[factual_exceeding, y, s, sc]
          n_cf_exceeding <- sum(cf_values > threshold, na.rm = TRUE)
          
          cf_count[y, sample_idx, sc] <- n_cf_exceeding
        }
      }
    }
    
    rm(pf, pcf, pf_country, pcf_country)
    gc()
  }
  
  cat("\n")
  
  list(
    factual_count = factual_count,
    cf_count = cf_count,
    scenario_names = scenarios
  )
}

#### Create country table - 2024 ONLY ####
create_country_table_2024 <- function(batch_dir, site_ids_all, site_country_map, aosis_members, eu27_uk_members, threshold = 0) {
  
  row_countries <- c("Global", "AOSIS", "Indonesia", "Japan", "China", "USA", "Australia")
  
  col_scenarios <- c("natural", "all", "Saudi Aramco", "ExxonMobil", 
                     "Chevron", "Holcim Group", "topten", "Gazprom", "National Iranian Oil Company", "BP", 
                     "China", "USA", "AOSIS", "Indonesia", "India", "Japan", 
                     "Russia", "EU27UK", "Australia")
  col_labels <- c("All anthropogenic", "All carbon majors", 
                  "Saudi Aramco", "ExxonMobil", "Chevron", "Holcim Group", "Top 10", "Gazprom", "National Iranian Oil Company", "BP", 
                  "China", "USA", "AOSIS", "Indonesia", "India", "Japan",
                  "Russia", "EU27+UK", "Australia")
  
  results <- data.frame(
    Country = row_countries,
    N_sites = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  for (r in seq_along(row_countries)) {
    country <- row_countries[r]
    
    cat("\n========================================\n")
    cat("Processing:", country, "\n")
    cat("========================================\n")
    
    # Get site IDs for this country
    if (country == "Global") {
      country_site_ids <- site_ids_all
    } else if (country == "AOSIS") {
      country_site_ids <- site_country_map %>%
        filter(Country_Name %in% aosis_members) %>%
        pull(Site_ID) %>%
        unique()
    } else if (country == "USA") {
      country_site_ids <- site_country_map %>%
        filter(Country_Name == "United States") %>%
        pull(Site_ID) %>%
        unique()
    } else {
      country_site_ids <- site_country_map %>%
        filter(Country_Name == country) %>%
        pull(Site_ID) %>%
        unique()
    }
    
    if (length(country_site_ids) == 0) {
      cat("  No sites found\n")
      results[r, "N_sites"] <- 0
      for (c in seq_along(col_labels)) {
        results[r, col_labels[c]] <- NA
      }
      next
    }
    
    cat("  Found", length(country_site_ids), "sites\n")
    results[r, "N_sites"] <- length(country_site_ids)
    
    # Get country data with conditional filtering
    country_data <- get_batch_counts_country(batch_dir, country_site_ids, site_ids_all, threshold)
    
    if (is.null(country_data)) {
      for (c in seq_along(col_labels)) {
        results[r, col_labels[c]] <- NA
      }
      next
    }
    
    # For each column scenario
    for (c in seq_along(col_scenarios)) {
      sc_idx <- which(country_data$scenario_names == col_scenarios[c])
      
      if (length(sc_idx) == 0) {
        results[r, col_labels[c]] <- NA
        next
      }
      
      # Extract 2024 (year 40)
      yr_idx <- 40
      factual_2024 <- country_data$factual_count[yr_idx, ]
      cf_2024 <- country_data$cf_count[yr_idx, , sc_idx]
      
      # Attribution: (factual - cf) / factual * 100
      attribution <- (factual_2024 - cf_2024) / factual_2024 * 100
      
      # Calculate statistics across 840 samples
      median_attr <- median(attribution, na.rm = TRUE)
      ci_lo_95 <- quantile(attribution, (1 - 0.95) / 2, na.rm = TRUE)  # 2.5%
      ci_hi_95 <- quantile(attribution, 1 - (1 - 0.95) / 2, na.rm = TRUE)  # 97.5%
      p1 <- quantile(attribution, 0.01, na.rm = TRUE)  # 1st percentile
      
      # Add star if 1st percentile > 0 (virtually certain contribution)
      star <- ifelse(p1 > 0, "*", "")
      
      # Format as: median* on one line, (ci_lo - ci_hi) on the next
      results[r, col_labels[c]] <- sprintf("%.1f%%%s\n(%.1f%% - %.1f%%)",
                                           median_attr, star, ci_lo_95, ci_hi_95)
    }
  }
  
  results
}

#### Create country table - CUMULATIVE (1985-2024) ####
create_country_table_cumulative <- function(batch_dir, site_ids_all, site_country_map, aosis_members, eu27_uk_members, threshold = 0) {
  
  row_countries <- c("Global", "AOSIS", "Indonesia", "Japan", "China", "USA", "Australia")
  
  col_scenarios <- c("natural", "all", "Saudi Aramco", "ExxonMobil", 
                     "Chevron", "Holcim Group", "topten", "Gazprom", "National Iranian Oil Company", "BP", 
                     "China", "USA", "AOSIS", "Indonesia", "India", "Japan", 
                     "Russia", "EU27UK", "Australia")
  col_labels <- c("All anthropogenic", "All carbon majors", 
                  "Saudi Aramco", "ExxonMobil", "Chevron", "Holcim Group", "Top 10", "Gazprom", "National Iranian Oil Company", "BP", 
                  "China", "USA", "AOSIS", "Indonesia", "India", "Japan",
                  "Russia", "EU27+UK", "Australia")
  
  results <- data.frame(
    Country = row_countries,
    N_sites = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  for (r in seq_along(row_countries)) {
    country <- row_countries[r]
    
    cat("\n========================================\n")
    cat("Processing:", country, "\n")
    cat("========================================\n")
    
    # Get site IDs for this country
    if (country == "Global") {
      country_site_ids <- site_ids_all
    } else if (country == "AOSIS") {
      country_site_ids <- site_country_map %>%
        filter(Country_Name %in% aosis_members) %>%
        pull(Site_ID) %>%
        unique()
    } else if (country == "USA") {
      country_site_ids <- site_country_map %>%
        filter(Country_Name == "United States") %>%
        pull(Site_ID) %>%
        unique()
    } else {
      country_site_ids <- site_country_map %>%
        filter(Country_Name == country) %>%
        pull(Site_ID) %>%
        unique()
    }
    
    if (length(country_site_ids) == 0) {
      cat("  No sites found\n")
      results[r, "N_sites"] <- 0
      for (c in seq_along(col_labels)) {
        results[r, col_labels[c]] <- NA
      }
      next
    }
    
    cat("  Found", length(country_site_ids), "sites\n")
    results[r, "N_sites"] <- length(country_site_ids)
    
    # Get country data with conditional filtering
    country_data <- get_batch_counts_country(batch_dir, country_site_ids, site_ids_all, threshold)
    
    if (is.null(country_data)) {
      for (c in seq_along(col_labels)) {
        results[r, col_labels[c]] <- NA
      }
      next
    }
    
    # For each column scenario
    for (c in seq_along(col_scenarios)) {
      sc_idx <- which(country_data$scenario_names == col_scenarios[c])
      
      if (length(sc_idx) == 0) {
        results[r, col_labels[c]] <- NA
        next
      }
      
      # Sum counts across all 40 years
      total_factual <- colSums(country_data$factual_count, na.rm = TRUE)
      total_cf <- colSums(country_data$cf_count[, , sc_idx], na.rm = TRUE)
      
      # Attribution: (factual - cf) / factual * 100
      attribution <- (total_factual - total_cf) / total_factual * 100
      
      # Calculate statistics across 840 samples
      median_attr <- median(attribution, na.rm = TRUE)
      ci_lo_95 <- quantile(attribution, (1 - 0.95) / 2, na.rm = TRUE)  # 2.5%
      ci_hi_95 <- quantile(attribution, 1 - (1 - 0.95) / 2, na.rm = TRUE)  # 97.5%
      p1 <- quantile(attribution, 0.01, na.rm = TRUE)  # 1st percentile
      
      # Add star if 1st percentile > 0 (virtually certain contribution)
      star <- ifelse(p1 > 0, "*", "")
      
      # Format as: median* on one line, (ci_lo - ci_hi) on the next
      results[r, col_labels[c]] <- sprintf("%.1f%%%s\n(%.1f%% - %.1f%%)",
                                           median_attr, star, ci_lo_95, ci_hi_95)
    }
  }
  
  results
}

#### Run for 2024 only ####
cat("\n========== GENERATING 2024 TABLES ==========\n")
country_table_10_2024 <- create_country_table_2024(batch_dir, site_ids_all, site_country_map, aosis_members, eu27_uk_members, threshold = 10)
country_table_20_2024 <- create_country_table_2024(batch_dir, site_ids_all, site_country_map, aosis_members, eu27_uk_members, threshold = 20)

#### Run for cumulative (1985-2024) ####
cat("\n========== GENERATING CUMULATIVE TABLES ==========\n")
country_table_10_cumulative <- create_country_table_cumulative(batch_dir, site_ids_all, site_country_map, aosis_members, eu27_uk_members, threshold = 10)
country_table_20_cumulative <- create_country_table_cumulative(batch_dir, site_ids_all, site_country_map, aosis_members, eu27_uk_members, threshold = 20)





library(kableExtra)
library(tinytex)

create_table_pdf <- function(df, output_file, caption = NULL) {

  # Replace NA with em-dash and fix column names for LaTeX
  df[is.na(df)] <- "---"
  colnames(df) <- gsub("_", " ", colnames(df))  # Replace underscores in headers

  bold_row <- which(df$Country == "All anthropogenic\nemissions")

  # Fold N sites into the Country cell (country name, then N = ... below it),
  # then drop the now-redundant N sites column.
  if ("N sites" %in% colnames(df)) {
    df$Country <- paste0(df$Country, "\nN = ", df[["N sites"]])
    df[["N sites"]] <- NULL
  }

  # Convert the embedded "\n" (median on one line, CI on the next, from the
  # sprintf in create_country_table_*) into a real LaTeX line break via
  # \makecell (loaded in the preamble below). linebreak() does NOT escape
  # special characters (e.g. "%", which starts a LaTeX comment and truncates
  # the rest of the row) -- escape_latex() handles that first, leaving the
  # embedded newline intact for linebreak() to convert. kbl() below must NOT
  # escape again or it will mangle the \makecell commands just inserted.
  df[] <- lapply(seq_along(df), function(j) {
    col <- df[[j]]
    if (!is.character(col)) return(col)
    col <- kableExtra:::escape_latex(col)
    kableExtra::linebreak(col, align = if (j == 1) "l" else "c")
  })

  # Create LaTeX table
  latex_table <- df %>%
    kbl(
      format   = "latex",
      booktabs = TRUE,
      caption  = caption,
      escape   = FALSE,   # linebreak() above already escaped + inserted \makecell
      align    = c("l", rep("c", ncol(df) - 1)),
      linesep  = ""
    ) %>%
    kable_styling(
      latex_options = c("HOLD_position"),
      font_size     = 9.5,
      position      = "center"
    ) %>%
    #row_spec(0, bold = TRUE) %>%
    row_spec(bold_row, bold = TRUE, hline_after = TRUE) %>%
    row_spec(1, bold = TRUE, hline_after = TRUE) %>%
    column_spec(1, bold = TRUE, width = "3.5cm") %>%
    footnote(
      #general           = "* virtually certain contribution (1st percentile > 0%). Values: median (95% likely range).",
      general_title     = "",
      footnote_as_chunk = TRUE,
      escape            = TRUE
    )
  
  tex_file <- gsub("\\.pdf$", ".tex", output_file)
  
  latex_doc <- paste0(
    "\\documentclass[10pt]{article}\n",
    "\\usepackage[landscape, margin=0.8cm]{geometry}\n",
    "\\usepackage{float}\n",
    "\\usepackage{graphicx}\n",
    "\\usepackage{booktabs}\n",
    "\\usepackage{array}\n",
    "\\usepackage{threeparttable}\n",
    "\\usepackage{makecell}\n",
    "\\usepackage{helvet}\n",
    "\\renewcommand{\\familydefault}{\\sfdefault}\n",
    "\\usepackage{xcolor}\n",
    "\\usepackage{colortbl}\n",
    "\\usepackage{caption}\n",
    "\\captionsetup[table]{font=small, labelfont=bf}\n",
    "\\begin{document}\n",
    "\\pagestyle{empty}\n",
    latex_table,
    "\n\\end{document}"
  )
  
  writeLines(latex_doc, tex_file)
  tinytex::pdflatex(tex_file)
  
  cat("Created:", output_file, "\n")
  
  # Clean up
  file.remove(tex_file)
  aux_files <- list.files(
    pattern = paste0(gsub("\\.tex$", "", basename(tex_file)), "\\.(aux|log)$")
  )
  if (length(aux_files) > 0) file.remove(aux_files)
}














transpose_table <- function(df) {
  # Save country names and N_sites
  countries <- df$Country
  n_sites   <- df$N_sites
  
  # Create header: "Country (N sites)"
  col_headers <- paste0(countries, "\n(N=", n_sites, ")")
  
  # Transpose the attribution columns only
  df_t <- df %>%
    select(-Country, -N_sites) %>%
    t() %>%
    as.data.frame()
  
  colnames(df_t) <- col_headers
  df_t <- cbind(Scenario = rownames(df_t), df_t)
  rownames(df_t) <- NULL
  
  df_t
}







create_table_pdf(
  transpose_table(country_table_10_2024),
  "table_country_10_2024.pdf"
  #caption = "Country-specific bleaching attribution for 2024 ($>$10\\% threshold)."
)

create_table_pdf(
  transpose_table(country_table_20_2024),
  "table_country_20_2024.pdf"
  #caption = "Country-specific bleaching attribution for 2024 ($>$20\\% threshold)."
)

create_table_pdf(
  transpose_table(country_table_10_cumulative),
  "table_country_10_cumulative.pdf"
  #caption = "Country-specific cumulative bleaching attribution 1985--2024 ($>$10\\% threshold)."
)

create_table_pdf(
  transpose_table(country_table_20_cumulative),
  "table_country_20_cumulative.pdf"
  #caption = "Country-specific cumulative bleaching attribution 1985--2024 ($>$20\\% threshold)."
)

