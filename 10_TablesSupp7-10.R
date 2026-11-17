
#### attribution tables for supp (Tables S7-S10)

library(dplyr)

#### Load counts ####
# v2_mixed impact model, all reef locations (~72,500 5 km cells). Per-country event counts are made on BSC
# from /esarchive/Coracle/bleaching_attribution_v2/reef_locations/ by 04_HPC_Scripts/reef_counts_v2.py.
# Every reef cell has a country (`admin` in /esarchive/Coracle/climate_data/reef_locations/cells.csv);
# the table rows are groups of those countries, defined in reef_counts_v2_helpers.R
# (reef_country_groups: full 39-member AOSIS, USA incl. territories, Australia incl. external territories).
source("reef_counts_v2_helpers.R")
counts_file <- Sys.getenv("REEF_COUNTS", "path/to/reef_counts_v2.nc")
# CELLSET = "sites_area" (survey-site grid cells, cos-latitude area weighted; default), "reefs_area" (all reef
# cells, area weighted), or unweighted "sites" / "reefs"; goes into the file names
cellset     <- Sys.getenv("CELLSET", "sites_area")
table_dir   <- Sys.getenv("TABLE_DIR", ".")
rc <- read_reef_counts(counts_file, cellset)

cat("Countries:", length(rc$country), "Years:", length(rc$year), "Scenarios:", length(rc$scenario),
    "Samples:", length(rc$sample), "\n")

#### Create country table - 2024 ONLY ####
create_country_table_2024 <- function(rc, threshold) {
  
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
    
    # Sum the counts of this row's countries (Global = all reef cells)
    country_data <- if (country == "Global") group_counts(rc, threshold) else
      group_counts(rc, threshold, reef_country_groups[[country]])
    results[r, "N_sites"] <- if (is.null(country_data)) 0 else country_data$n_cells   # number of reef cells
    cat("  Reef cells:", results[r, "N_sites"], "\n")
    
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
create_country_table_cumulative <- function(rc, threshold) {
  
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
    
    # Sum the counts of this row's countries (Global = all reef cells)
    country_data <- if (country == "Global") group_counts(rc, threshold) else
      group_counts(rc, threshold, reef_country_groups[[country]])
    results[r, "N_sites"] <- if (is.null(country_data)) 0 else country_data$n_cells   # number of reef cells
    cat("  Reef cells:", results[r, "N_sites"], "\n")
    
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
country_table_10_2024 <- create_country_table_2024(rc, threshold = 10)
country_table_30_2024 <- create_country_table_2024(rc, threshold = 30)

#### Run for cumulative (1985-2024) ####
cat("\n========== GENERATING CUMULATIVE TABLES ==========\n")
country_table_10_cumulative <- create_country_table_cumulative(rc, threshold = 10)
country_table_30_cumulative <- create_country_table_cumulative(rc, threshold = 30)







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
  file.path(table_dir, paste0("table_country_10_2024_", cellset, ".pdf"))
  #caption = "Country-specific bleaching attribution for 2024 ($>$10\\% threshold)."
)

create_table_pdf(
  transpose_table(country_table_30_2024),
  file.path(table_dir, paste0("table_country_30_2024_", cellset, ".pdf"))
  #caption = "Country-specific bleaching attribution for 2024 ($>$30\\% threshold)."
)

create_table_pdf(
  transpose_table(country_table_10_cumulative),
  file.path(table_dir, paste0("table_country_10_cumulative_", cellset, ".pdf"))
  #caption = "Country-specific cumulative bleaching attribution 1985--2024 ($>$10\\% threshold)."
)

create_table_pdf(
  transpose_table(country_table_30_cumulative),
  file.path(table_dir, paste0("table_country_30_cumulative_", cellset, ".pdf"))
  #caption = "Country-specific cumulative bleaching attribution 1985--2024 ($>$30\\% threshold)."
)




