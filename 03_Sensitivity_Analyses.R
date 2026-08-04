### author: Puja Pande
### date: 12 May 2026 - updated and new
### description: the script runs a variety of sensitivity test to see how the results change with different parameters and assumptions
### covers: model type comparisons, MPA interactions, cluster/lag/polynomial (quadratic,
###         cubic) sensitivity, Conley (200km spherical) SE checks, seasonality FE,
###         and DHW specification (DHW vs DHW_adj) sensitivity, with summary tables



#### load libraries ####
if (!requireNamespace("kableExtra", quietly = TRUE)) install.packages("kableExtra")
if (!requireNamespace("tinytex", quietly = TRUE)) install.packages("tinytex")
if (tinytex::is_tinytex() && !("siunitx" %in% tinytex::tl_pkgs())) tinytex::tlmgr_install("siunitx")
if (tinytex::is_tinytex() && !("adjustbox" %in% tinytex::tl_pkgs())) tinytex::tlmgr_install(c("adjustbox", "collectbox"))

library(dplyr)
library(fixest)
library(ggplot2)
library(tidyr)
library(patchwork)
library(plotly)
library(sf)
library(glmmTMB)
library(modelsummary)
library(flextable)
library(officer)
library(kableExtra)

# modelsummary defaults to the 'tinytable' backend for LaTeX; make_latex_table()
# below relies on kableExtra::add_header_above()/footnote(), so force that backend.
options(modelsummary_factory_latex = "kableExtra")


#### load datasets ####
All_Bleaching_Events_Data_AllDHW <- read.csv("/Users/pujapande/Library/CloudStorage/Dropbox/Coracle/Datasets - 01 July 2026/Final_Panel_Data_allDHW_01July.csv")
All_Bleaching_Events_Data_AllDHW$Proportion_Bleached <- All_Bleaching_Events_Data_AllDHW$Percent_Bleached / 100
All_Bleaching_Events_Data_AllDHW$abs_lat             <- abs(All_Bleaching_Events_Data_AllDHW$Latitude_Degrees)
All_Bleaching_Events_Data_AllDHW$mass_bleaching      <- ifelse(All_Bleaching_Events_Data_AllDHW$Percent_Bleached >= 30, 1, 0)
All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
  mutate(
    Ecoregion_Month = interaction(Ecoregion_Name, Date_Month, drop = TRUE)
  )

All_Bleaching_Events_Data_AllDHW_Max <- read.csv("/Users/pujapande/Downloads/Final_Combined_Data_all_levels_01July2026_extravars.csv")
All_Bleaching_Events_Data_AllDHW_Max <- left_join(All_Bleaching_Events_Data_AllDHW, All_Bleaching_Events_Data_AllDHW_Max[, c(1:8, 17:34)], by = c("Site_ID", "Latitude_Degrees", "Longitude_Degrees", "Date_Year", "Date_Month", "Date_Day", "Ecoregion_Name", "Percent_Bleached"))


#### different model types ####
model_linear <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  cluster = ~Ecoregion_Name, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_lat <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

model_linear_turbidity <- feols(
  Percent_Bleached ~ dhw + dhw:Turbidity | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)

model_linear_depth <- feols(
  Percent_Bleached ~ dhw + dhw:Depth_m | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW
)

model_linear_dhw_std <- feols(
  Percent_Bleached ~ dhw + dhw:DHW_1985.2005_std | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW_Max
)

model_linear_hotspot_std <- feols(
  Percent_Bleached ~ dhw + dhw:hotspot_1985.2005_std | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW_Max
)

model_linear_hotspot_warming <- feols(
  Percent_Bleached ~ dhw + dhw:hotspot_warming | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW_Max
)

model_linear_lat_depth <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat + dhw:Depth_m | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW_Max
)

model_linear_lat_hotspot <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat + dhw:hotspot_1985.2005_std | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data = All_Bleaching_Events_Data_AllDHW_Max
)

#### Adding MPA interactions to bleaching events data ####
# Read all files
poly_files <- c(
  "~/Library/CloudStorage/Dropbox/Coracle/Maps and Spatial Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_0/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp",
  "~/Library/CloudStorage/Dropbox/Coracle/Maps and Spatial Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_1/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp",
  "~/Library/CloudStorage/Dropbox/Coracle/Maps and Spatial Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_2/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp"
)

point_files <- c(
  "~/Library/CloudStorage/Dropbox/Coracle/Maps and Spatial Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_0/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp",
  "~/Library/CloudStorage/Dropbox/Coracle/Maps and Spatial Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_1/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp",
  "~/Library/CloudStorage/Dropbox/Coracle/Maps and Spatial Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_2/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp"
)

# Merge polygons
poly_list <- lapply(poly_files, st_read)
merged_poly <- do.call(rbind, poly_list)

# Merge points
point_list <- lapply(point_files, st_read)
merged_points <- do.call(rbind, point_list)

# Get all unique column names
all_cols <- union(names(merged_poly), names(merged_points))

# Add missing columns to each dataset
for (col in all_cols) {
  if (!col %in% names(merged_poly)) {
    merged_poly[[col]] <- NA
  }
  if (!col %in% names(merged_points)) {
    merged_points[[col]] <- NA
  }
}

# Reorder columns to match
merged_poly <- merged_poly[, all_cols]
merged_points <- merged_points[, all_cols]

# Now combine
all_marine <- rbind(merged_poly, merged_points)

# Convert your data to sf object (points)
merged_sf <- st_as_sf(All_Bleaching_Events_Data_AllDHW, 
                      coords = c("Longitude_Degrees", "Latitude_Degrees"),  # Adjust column names
                      crs = 4326)  # WGS84

# Fix invalid geometries in MPA data
all_marine <- st_make_valid(all_marine)

# Spatial join
sf_use_s2(FALSE)
bleaching_with_mpa <- st_join(merged_sf, all_marine["SITE_ID"], left = TRUE)
sf_use_s2(TRUE)

# Add binary MPA indicator
bleaching_with_mpa$in_mpa <- !is.na(bleaching_with_mpa$SITE_ID)

# Drop geometry and deduplicate
# st_join can create duplicates if a point falls in multiple MPAs
All_Bleaching_Events_Data_AllDHW_MPA <- bleaching_with_mpa %>%
  st_drop_geometry() %>%
  group_by(across(-c(SITE_ID, in_mpa))) %>%
  summarise(in_mpa = any(in_mpa), .groups = "drop")

#### Run MPA model ####
model_linear_mpa <- feols(
  Percent_Bleached ~ dhw + dhw:in_mpa | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW_MPA
)


#### cluster sensitivity test ####
model_linear_nocluster <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_ecoregion <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  cluster = ~Ecoregion_Name, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_site <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
  cluster = ~Site_ID, 
  data = All_Bleaching_Events_Data_AllDHW)


#### lag sensitivity test ####
model_lagminus3 <- feols(Percent_Bleached ~ dhw_lag.3 + dhw_lag.3:abs_lat | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                         data = All_Bleaching_Events_Data_AllDHW)
model_lagminus2 <- feols(Percent_Bleached ~ dhw_lag.2 + dhw_lag.2:abs_lat | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                         data = All_Bleaching_Events_Data_AllDHW)
model_lagminus1 <- feols(Percent_Bleached ~ dhw_lag.1 + dhw_lag.1:abs_lat | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                         data = All_Bleaching_Events_Data_AllDHW)
model_lag0 <- feols(Percent_Bleached ~ dhw_lag0 + dhw_lag0:abs_lat | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)
model_lag1 <- feols(Percent_Bleached ~ dhw_lag1 + dhw_lag1:abs_lat | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)
model_lag2 <- feols(Percent_Bleached ~ dhw_lag2 + dhw_lag2:abs_lat | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)
model_lag3 <- feols(Percent_Bleached ~ dhw_lag3 + dhw_lag3:abs_lat| Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)
model_lag4 <- feols(Percent_Bleached ~ dhw_lag4 + dhw_lag4:abs_lat| Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)
model_lag5 <- feols(Percent_Bleached ~ dhw_lag5 + dhw_lag5:abs_lat| Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
                    data = All_Bleaching_Events_Data_AllDHW)



#### Quadratic DHW models ####
model_quadratic <- feols(
  Percent_Bleached ~ dhw + I(dhw^2) | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

model_quadratic_lat <- feols(
  Percent_Bleached ~ dhw + I(dhw^2) + dhw:abs_lat + I(dhw^2):abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

#### Cubic DHW models ####
model_cubic <- feols(
  Percent_Bleached ~ dhw + I(dhw^2) + I(dhw^3) | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

model_cubic_lat <- feols(
  Percent_Bleached ~ dhw + I(dhw^2) + I(dhw^3) +
    dhw:abs_lat + I(dhw^2):abs_lat + I(dhw^3):abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)






################## conley errors ################
#### Compute Conley (200km, spherical) vcov for each model ####
conley_vcov <- function(model) {
  vcov_conley(model,
              lat      = "Latitude_Degrees",
              lon      = "Longitude_Degrees",
              cutoff   = 200,
              distance = "spherical")
}

#### Different model types ####
model_linear_vcov <- conley_vcov(model_linear)
model_linear_lat_vcov <- conley_vcov(model_linear_lat)
model_linear_turbidity_vcov <- conley_vcov(model_linear_turbidity)
model_linear_depth_vcov <- conley_vcov(model_linear_depth)
model_linear_dhw_std_vcov <- conley_vcov(model_linear_dhw_std)
model_linear_hotspot_std_vcov <- conley_vcov(model_linear_hotspot_std)
model_linear_hotspot_warming_vcov <- conley_vcov(model_linear_hotspot_warming)
model_linear_lat_depth_vcov <- conley_vcov(model_linear_lat_depth)
model_linear_lat_hotspot_vcov <- conley_vcov(model_linear_lat_hotspot)

All_Bleaching_Events_Data_AllDHW_MPA <- All_Bleaching_Events_Data_AllDHW_MPA %>%
  left_join(
    All_Bleaching_Events_Data_AllDHW %>% 
      select(Site_ID, Date_Year, Latitude_Degrees, Longitude_Degrees),
    by = c("Site_ID", "Date_Year")
  )
model_linear_mpa_vcov <- conley_vcov(model_linear_mpa)


#### Lag models ####
model_lagminus3_vcov <- conley_vcov(model_lagminus3)
model_lagminus2_vcov <- conley_vcov(model_lagminus2)
model_lagminus1_vcov <- conley_vcov(model_lagminus1)
model_lag0_vcov <- conley_vcov(model_lag0)
model_lag1_vcov <- conley_vcov(model_lag1)
model_lag2_vcov <- conley_vcov(model_lag2)
model_lag3_vcov <- conley_vcov(model_lag3)
model_lag4_vcov <- conley_vcov(model_lag4)
model_lag5_vcov <- conley_vcov(model_lag5)


#### Quadratic models ####
model_quadratic_vcov <- conley_vcov(model_quadratic)
model_quadratic_lat_vcov <- conley_vcov(model_quadratic_lat)

#### Cubic models ####
model_cubic_vcov <- conley_vcov(model_cubic)
model_cubic_lat_vcov <- conley_vcov(model_cubic_lat)



#### Custom styling function ####
# Mimics the classic 3-line journal regression table (e.g. Table S2 style):
# italicized "Dependent variable:" spanning header, rule under headers,
# rule separating coefficients from Observations/R2 block, no vertical rules.
nature_table_style <- function(ft, dv_label = "Dependent variable: Percent Bleached (%)") {
  ft <- ft %>%
    add_header_row(values = c("", dv_label), colwidths = c(1, ncol_keys(ft) - 1), top = TRUE)

  header_rows <- nrow_part(ft, "header")

  gof_labels <- c("Num.Obs.", "R2", "R2 Adj.", "AIC", "BIC", "RMSE", "Std.Errors")
  body_col1  <- ft$body$dataset[[1]]
  gof_row    <- which(body_col1 %in% gof_labels)

  ft <- ft %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 11, part = "all") %>%
    italic(part = "header") %>%
    align(align = "center", part = "header") %>%
    align(j = 1, align = "left", part = "body") %>%
    align(j = 2:ncol_keys(ft), align = "center", part = "body") %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 1.5), part = "all") %>%
    hline_bottom(border = fp_border(width = 1.5), part = "all") %>%
    hline(i = header_rows, border = fp_border(width = 1), part = "header")

  if (length(gof_row) > 0) {
    ft <- ft %>% hline(i = min(gof_row) - 1, border = fp_border(width = 1), part = "body")
  }

  ft %>%
    autofit() %>%
    width(width = 1.2, unit = "in")
}

#### Add a small gap between different coefficients' rows ####
# Each retained term occupies 2 lines (estimate, then std. error). Coefficient
# rows sit between the first \midrule (below the header) and the second
# \midrule (above Num.Obs./R2/etc.). This inserts a booktabs \addlinespace
# after every term's pair of rows (except the last), so different
# explanatory variables are visually separated without a full rule.
add_coef_row_spacing <- function(tex_text) {
  lines <- strsplit(tex_text, "\n")[[1]]
  midrule_idx <- grep("^\\\\midrule\\s*$", lines)
  if (length(midrule_idx) < 2) return(tex_text)

  coef_lines <- (midrule_idx[1] + 1):(midrule_idx[2] - 1)
  n <- length(coef_lines)
  if (n < 4) return(tex_text)

  insert_after <- coef_lines[seq(2, n - 1, by = 2)]
  for (idx in rev(insert_after)) {
    lines <- append(lines, "\\addlinespace", after = idx)
  }
  paste(lines, collapse = "\n")
}

#### Relabel GOF rows (Num.Obs. -> Observations, R2 -> R^2) ####
# Applied as a text substitution after modelsummary builds the table, rather
# than via gof_map, so the gof_omit regexes each call already uses (which
# match on the default clean names like "R2 Adj") keep working unchanged.
rename_gof_labels <- function(tex_text) {
  tex_text <- gsub("Num\\.Obs\\.", "Observations", tex_text)
  tex_text <- gsub("\\bR2\\b", "$R^2$", tex_text)
  tex_text
}

#### Auto-scale wide tables to fit the page ####
# Wraps just the tabular in an adjustbox that shrinks it to \textwidth only
# if it would otherwise overflow (never enlarges a table that already fits),
# so wide tables (many model columns, long coefficient labels) stay on the
# page without having to hand-tune font size or column widths per table.
wrap_table_scale <- function(tex_text) {
  tex_text <- sub("\\\\begin\\{tabular\\}", "\\\\begin{adjustbox}{max width=\\\\textwidth}\n\\\\begin{tabular}", tex_text)
  tex_text <- sub("\\\\end\\{tabular\\}", "\\\\end{tabular}\n\\\\end{adjustbox}", tex_text)
  tex_text
}

#### Add a (1) (2) (3)... row under the model-name header row ####
# Matches the journal example, which numbers each column just above the rule
# that separates the header from the coefficient rows.
add_column_numbers <- function(tex_text, n_models) {
  lines <- strsplit(tex_text, "\n")[[1]]
  midrule_idx <- grep("^\\\\midrule\\s*$", lines)[1]
  if (is.na(midrule_idx)) return(tex_text)

  num_row <- paste0("  & ", paste(paste0("(", seq_len(n_models), ")"), collapse = " & "), "\\\\")
  lines <- append(lines, num_row, after = midrule_idx - 1)
  paste(lines, collapse = "\n")
}

#### Custom LaTeX styling function ####
# Builds a modelsummary table directly to booktabs LaTeX, styled like the
# journal example (Table S2): italic "Dependent variable:" header spanning
# all model columns, rule under the header, no vertical rules, notes as a
# threeparttable footnote below the table.
make_latex_table <- function(models, vcov_list, coef_rename, gof_omit,
                              notes, title, file,
                              dv_label = "Percent Bleached (\\%)",
                              stars = c('*' = 0.05, '**' = 0.01)) {

  tab <- modelsummary(
    models,
    vcov        = vcov_list,
    stars       = stars,
    fmt         = 2,
    statistic   = "({std.error})",
    coef_omit   = "Intercept",
    coef_rename = coef_rename,
    gof_omit    = gof_omit,
    title       = title,
    output      = "latex",
    escape      = FALSE
  )

  n_models <- length(models)

  tab <- tab %>%
    kableExtra::add_header_above(
      c(" " = 1, setNames(n_models, paste0("\\textit{Dependent variable: ", dv_label, "}"))),
      escape = FALSE, line = TRUE
    ) %>%
    kableExtra::kable_styling(latex_options = "hold_position") %>%
    kableExtra::footnote(
      general         = unlist(notes),
      general_title   = "",
      threeparttable  = TRUE,
      escape          = FALSE
    )

  tab_text <- as.character(tab)
  tab_text <- add_coef_row_spacing(tab_text)
  tab_text <- rename_gof_labels(tab_text)
  tab_text <- add_column_numbers(tab_text, n_models)
  tab_text <- wrap_table_scale(tab_text)
  writeLines(tab_text, file)

  compile_latex_table_pdf(file)

  invisible(tab)
}

#### Compile a standalone table .tex (as written by make_latex_table) to PDF ####
# Wraps the raw \begin{table}...\end{table} fragment in a minimal document with
# the packages modelsummary/kableExtra need (booktabs, threeparttable, siunitx),
# plus T1 fontenc so glyphs like "|" render correctly, then renders via TinyTeX.
compile_latex_table_pdf <- function(tex_file) {
  if (!requireNamespace("tinytex", quietly = TRUE)) {
    warning("Package 'tinytex' not installed; skipping PDF compilation for ", tex_file)
    return(invisible(NULL))
  }

  wrapper <- tempfile(fileext = ".tex")
  writeLines(c(
    "\\documentclass[12pt]{article}",
    "\\usepackage[T1]{fontenc}",
    "\\usepackage[a4paper, margin=1in]{geometry}",
    "\\usepackage{booktabs}",
    "\\usepackage{threeparttable}",
    "\\usepackage{siunitx}",
    "\\usepackage[normalem]{ulem}",
    "\\usepackage{adjustbox}",
    "\\pagestyle{empty}",
    "\\begin{document}",
    paste0("\\input{", normalizePath(tex_file, mustWork = TRUE), "}"),
    "\\end{document}"
  ), wrapper)

  pdf_out <- tinytex::pdflatex(wrapper, clean = TRUE)
  file.copy(pdf_out, sub("\\.tex$", ".pdf", tex_file), overwrite = TRUE)
  file.remove(wrapper, pdf_out)
  invisible(NULL)
}

#### Interaction models ####
make_latex_table(
  models = list(
    "Latitude Interaction"      = model_linear_lat,
    "Turbidity Interaction"     = model_linear_turbidity,
    "Depth Interaction"         = model_linear_depth,
    "\\shortstack{DHW Anomaly \\\\ Variability Interaction}"                  = model_linear_dhw_std,
    "\\shortstack{SST Anomaly \\\\ Variability Interaction}"              = model_linear_hotspot_std,
    "\\shortstack{SST Anomaly\\\\ Warming Interaction}"      = model_linear_hotspot_warming,
    "\\shortstack{Latitude and SST \\\\ Anomaly Variability Interaction}" = model_linear_lat_hotspot,
    "\\shortstack{Marine Protected Area \\\\ Interaction}"                       = model_linear_mpa
  ),
  vcov_list = list(
    model_linear_lat_vcov,
    model_linear_turbidity_vcov,
    model_linear_depth_vcov,
    model_linear_dhw_std_vcov,
    model_linear_hotspot_std_vcov,
    model_linear_hotspot_warming_vcov,
    model_linear_lat_hotspot_vcov,
    model_linear_mpa_vcov
  ),
  coef_rename = c(
    "dhw"                       = "DHW",
    "dhw:abs_lat"               = "\\shortstack[l]{DHW × \\\\ |Latitude|}",
    "dhw:Turbidity"             = "\\shortstack[l]{DHW × \\\\ Turbidity}",
    "dhw:Depth_m"               = "\\shortstack[l]{DHW × \\\\ Depth}",
    "dhw:DHW_1985.2005_std"     = "\\shortstack[l]{DHW × \\\\ DHW Anomaly}",
    "dhw:hotspot_1985.2005_std" = "\\shortstack[l]{DHW × \\\\ SST Anomaly Variability}",
    "dhw:hotspot_warming"       = "\\shortstack[l]{DHW × \\\\ SST Anomaly Warming}",
    "dhw:in_mpaTRUE"            = "\\shortstack[l]{DHW × \\\\ Marine Protected Area}"
  ),
  gof_omit = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes = list(),
  title = "",
  file  = "table_interactions_nature.tex"
)


#### Lag sensitivity ####
make_latex_table(
  models = list(
    "Lead 3" = model_lagminus3, "Lead 2" = model_lagminus2, "Lead 1" = model_lagminus1,
    "No Lag or Lead"  = model_lag0,
    "Lag 1" = model_lag1, "Lag 2" = model_lag2, "Lag 3" = model_lag3,
    "Lag 4" = model_lag4, "Lag 5" = model_lag5
  ),
  vcov_list = list(
    model_lagminus3_vcov, model_lagminus2_vcov, model_lagminus1_vcov,
    model_lag0_vcov,
    model_lag1_vcov, model_lag2_vcov, model_lag3_vcov,
    model_lag4_vcov, model_lag5_vcov
  ),
  coef_rename = c(
    "dhw_lag.3" = "DHW Lead 3", "dhw_lag.2" = "DHW Lead 2", "dhw_lag.1" = "DHW Lead 1",
    "dhw_lag0"  = "DHW",
    "dhw_lag1"  = "DHW Lag 1", "dhw_lag2"  = "DHW Lag 2", "dhw_lag3"  = "DHW Lag 3",
    "dhw_lag4"  = "DHW Lag 4", "dhw_lag5"  = "DHW Lag 5",
    "dhw_lag.3:abs_lat" = "DHW Lead 3 × |Latitude|",
    "dhw_lag.2:abs_lat" = "DHW Lead 2 × |Latitude|",
    "dhw_lag.1:abs_lat" = "DHW Lead 1 × |Latitude|",
    "dhw_lag0:abs_lat"  = "DHW × |Latitude|",
    "dhw_lag1:abs_lat"  = "DHW Lag 1 × |Latitude|",
    "dhw_lag2:abs_lat"  = "DHW Lag 2 × |Latitude|",
    "dhw_lag3:abs_lat"  = "DHW Lag 3 × |Latitude|",
    "dhw_lag4:abs_lat"  = "DHW Lag 4 × |Latitude|",
    "dhw_lag5:abs_lat"  = "DHW Lag 5 × |Latitude|"
  ),
  gof_omit = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes = list(),
  title = "",
  file  = "table_lags_nature.tex"
)

#### Polynomial sensitivity ####
make_latex_table(
  models = list(
    "Linear with Latitude Interaction"        = model_linear_lat,
    "Linear"        = model_linear,
    "Quadratic with Latitude Interaction" = model_quadratic_lat,
    "Quadratic"     = model_quadratic,
    "Cubic with Latitude Interaction"     = model_cubic_lat,
    "Cubic"         = model_cubic
  ),
  vcov_list = list(
    model_linear_lat_vcov,
    model_linear_vcov,
    model_quadratic_lat_vcov,
    model_quadratic_vcov,
    model_cubic_lat_vcov,
    model_cubic_vcov
  ),
  coef_rename = c(
    "dhw"                 = "DHW",
    "I(I(dhw^2))"         = "DHW²",
    "I(I(dhw^3))"         = "DHW³",
    "dhw:abs_lat"         = "DHW × |Latitude|",
    "I(I(dhw^2)):abs_lat" = "DHW² × |Latitude|",
    "I(I(dhw^3)):abs_lat" = "DHW³ × |Latitude|"
  ),
  gof_omit = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes = list(),
  title = "",
  file  = "table_polynomial_sensitivity.tex"
)






model_linear_nocluster <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat| Site_ID + Date_Year + Ecoregion_Month, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_ecoregion <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat  | Site_ID + Date_Year + Ecoregion_Month, 
  cluster = ~Ecoregion_Name, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_site <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month, 
  cluster = ~Site_ID, 
  data = All_Bleaching_Events_Data_AllDHW)

model_linear_nocluster_vcov <- conley_vcov(model_linear_nocluster)
model_linear_ecoregion_vcov <- conley_vcov(model_linear_ecoregion)
model_linear_site_vcov <- conley_vcov(model_linear_site)


make_latex_table(
  models = list(
    "Conley 200km Clustering" = model_linear_lat,
    "Ecoregion Clustering"    = model_linear_ecoregion,
    "Site Clustering"         = model_linear_site,
    "No Clustering"           = model_linear_nocluster
  ),
  vcov_list = list(
    model_linear_lat_vcov,
    ~Ecoregion_Name,
    ~Site_ID,
    "iid"
  ),
  coef_rename = c("dhw" = "DHW", "dhw:abs_lat" = "DHW × |Latitude|"),
  gof_omit = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  stars = c('*' = 0.05, '**' = 0.01),
  notes = "",
  title = "",
  file  = "table_clustering_nature.tex"
)






#### seasonality #### 
#### Model with seasonality FE ####
model_linear_lat_season <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

model_linear_lat_2FE <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year,
  cluster = ~Ecoregion_Name,
  data    = All_Bleaching_Events_Data_AllDHW
)

model_linear_lat_vcov_2FE        <- conley_vcov(model_linear_lat_2FE)
model_linear_lat_season_vcov <- conley_vcov(model_linear_lat_season)

#### Comparison table ####
make_latex_table(
  models = list(
    "Site + Year + Seasonality Fixed Effects" = model_linear_lat_season,
    "Site + Year Fixed Effects"              = model_linear_lat_2FE
  ),
  vcov_list = list(
    model_linear_lat_season_vcov,
    model_linear_lat_vcov_2FE
  ),
  coef_rename = c(
    "dhw"         = "DHW",
    "dhw:abs_lat" = "DHW × |Latitude|"
  ),
  gof_omit = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes = list(
  ),
  title = "",
  file  = "table_seasonality_sensitivity.tex"
)









#### DHW specification sensitivity (DHW vs DHW_adj) ####

All_Bleaching_Events_Data_AllDHW_test <- All_Bleaching_Events_Data_AllDHW
All_Bleaching_Events_Data_AllDHW_test <- left_join(All_Bleaching_Events_Data_AllDHW_test, DF[, c(1, 2, 3, 4, 5, 6, 18)])
#### Fit base + adjusted DHW models (Conley SEs computed post-estimation) ####
model_linear_dhwspec <- feols(
  Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month,
  data = All_Bleaching_Events_Data_AllDHW_Max)
model_linear_dhwspec_vcov <- conley_vcov(model_linear_dhwspec)

model_linear_lat_dhwspec <- feols(
  Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  data = All_Bleaching_Events_Data_AllDHW_Max)
model_linear_lat_dhwspec_vcov <- conley_vcov(model_linear_lat_dhwspec)

model_linear_adj <- feols(
  Percent_Bleached ~ DHW_adj | Site_ID + Date_Year + Ecoregion_Month,
  data = All_Bleaching_Events_Data_AllDHW_Max)
model_linear_adj_vcov <- conley_vcov(model_linear_adj)

model_linear_adj_lat <- feols(
  Percent_Bleached ~ DHW_adj + DHW_adj:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
  data = All_Bleaching_Events_Data_AllDHW_Max)
model_linear_adj_lat_vcov <- conley_vcov(model_linear_adj_lat)

#### Table: DHW specification sensitivity ####
make_latex_table(
  models = list(
    "\\shortstack{DHW with\\\\Latitude Interaction}" = model_linear_lat_dhwspec,
    "DHW" = model_linear_dhwspec,
    "Adjusted DHW" = model_linear_adj,
    "\\shortstack{Adjusted DHW with\\\\Latitude Interaction}" = model_linear_adj_lat
  ),
  vcov_list = list(
    model_linear_lat_dhwspec_vcov,
    model_linear_dhwspec_vcov,
    model_linear_adj_vcov,
    model_linear_adj_lat_vcov
  ),
  coef_rename = c(
    "dhw"             = "DHW",
    "DHW_adj"         = "Adjusted DHW",
    "dhw:abs_lat"     = "DHW × |Latitude|",
    "DHW_adj:abs_lat" = "Adjusted DHW × |Latitude|"
  ),
  gof_omit = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
  notes = list(),
  title = "",
  file  = "table_sensitivity_dhw_spec.tex"
)



# #### load libraries ####
# library(dplyr)
# library(fixest)
# library(ggplot2)
# library(tidyr)
# library(patchwork)
# library(plotly)
# library(sf)
# library(glmmTMB)
# library(modelsummary)
# library(flextable)
# library(officer)


# #### load datasets ####
# All_Bleaching_Events_Data_AllDHW <- read.csv("path/to/Coracle/Datasets/Final_Panel_Data_allDHW_01July.csv")
# All_Bleaching_Events_Data_AllDHW$Proportion_Bleached <- All_Bleaching_Events_Data_AllDHW$Percent_Bleached / 100
# All_Bleaching_Events_Data_AllDHW$abs_lat             <- abs(All_Bleaching_Events_Data_AllDHW$Latitude_Degrees)
# All_Bleaching_Events_Data_AllDHW$mass_bleaching      <- ifelse(All_Bleaching_Events_Data_AllDHW$Percent_Bleached >= 30, 1, 0)
# All_Bleaching_Events_Data_AllDHW <- All_Bleaching_Events_Data_AllDHW %>%
#   mutate(
#     Ecoregion_Month = interaction(Ecoregion_Name, Date_Month, drop = TRUE)
#   )

# All_Bleaching_Events_Data_AllDHW_Max <- read.csv("path/to/data/Final_Combined_Data_all_levels_01July2026_extravars.csv")
# All_Bleaching_Events_Data_AllDHW_Max <- left_join(All_Bleaching_Events_Data_AllDHW, All_Bleaching_Events_Data_AllDHW_Max[, c(1:8, 17:34)], by = c("Site_ID", "Latitude_Degrees", "Longitude_Degrees", "Date_Year", "Date_Month", "Date_Day", "Ecoregion_Name", "Percent_Bleached"))


# #### different model types ####
# model_linear <- feols(
#   Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
#   cluster = ~Ecoregion_Name, 
#   data = All_Bleaching_Events_Data_AllDHW)

# model_linear_turbidity <- feols(
#   Percent_Bleached ~ dhw + dhw:Turbidity | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data = All_Bleaching_Events_Data_AllDHW
# )

# model_linear_depth <- feols(
#   Percent_Bleached ~ dhw + dhw:Depth_m | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data = All_Bleaching_Events_Data_AllDHW
# )

# model_linear_dhw_std <- feols(
#   Percent_Bleached ~ dhw + dhw:DHW_1985.2005_std | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data = All_Bleaching_Events_Data_AllDHW_Max
# )

# model_linear_hotspot_std <- feols(
#   Percent_Bleached ~ dhw + dhw:hotspot_1985.2005_std | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data = All_Bleaching_Events_Data_AllDHW_Max
# )

# model_linear_hotspot_warming <- feols(
#   Percent_Bleached ~ dhw + dhw:hotspot_warming | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data = All_Bleaching_Events_Data_AllDHW_Max
# )

# model_linear_lat_depth <- feols(
#   Percent_Bleached ~ dhw + dhw:abs_lat + dhw:Depth_m | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data = All_Bleaching_Events_Data_AllDHW_Max
# )

# model_linear_lat_hotspot <- feols(
#   Percent_Bleached ~ dhw + dhw:abs_lat + dhw:hotspot_1985.2005_std | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data = All_Bleaching_Events_Data_AllDHW_Max
# )

# #### Adding MPA interactions to bleaching events data ####
# # Read all files
# poly_files <- c(
#   "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_0/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp",
#   "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_1/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp",
#   "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_2/WDPA_WDOECM_Jan2026_Public_marine_shp-polygons.shp"
# )

# point_files <- c(
#   "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_0/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp",
#   "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_1/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp",
#   "path/to/Coracle/Maps_and_Spatial_Data/WDPA_WDOECM_Jan2026_Public_marine_shp/WDPA_WDOECM_Jan2026_Public_marine_shp_2/WDPA_WDOECM_Jan2026_Public_marine_shp-points.shp"
# )

# # Merge polygons
# poly_list <- lapply(poly_files, st_read)
# merged_poly <- do.call(rbind, poly_list)

# # Merge points
# point_list <- lapply(point_files, st_read)
# merged_points <- do.call(rbind, point_list)

# # Get all unique column names
# all_cols <- union(names(merged_poly), names(merged_points))

# # Add missing columns to each dataset
# for (col in all_cols) {
#   if (!col %in% names(merged_poly)) {
#     merged_poly[[col]] <- NA
#   }
#   if (!col %in% names(merged_points)) {
#     merged_points[[col]] <- NA
#   }
# }

# # Reorder columns to match
# merged_poly <- merged_poly[, all_cols]
# merged_points <- merged_points[, all_cols]

# # Now combine
# all_marine <- rbind(merged_poly, merged_points)

# # Convert your data to sf object (points)
# merged_sf <- st_as_sf(All_Bleaching_Events_Data_AllDHW, 
#                       coords = c("Longitude_Degrees", "Latitude_Degrees"),  # Adjust column names
#                       crs = 4326)  # WGS84

# # Fix invalid geometries in MPA data
# all_marine <- st_make_valid(all_marine)

# # Spatial join
# sf_use_s2(FALSE)
# bleaching_with_mpa <- st_join(merged_sf, all_marine["SITE_ID"], left = TRUE)
# sf_use_s2(TRUE)

# # Add binary MPA indicator
# bleaching_with_mpa$in_mpa <- !is.na(bleaching_with_mpa$SITE_ID)

# # Drop geometry and deduplicate
# # st_join can create duplicates if a point falls in multiple MPAs
# All_Bleaching_Events_Data_AllDHW_MPA <- bleaching_with_mpa %>%
#   st_drop_geometry() %>%
#   group_by(across(-c(SITE_ID, in_mpa))) %>%
#   summarise(in_mpa = any(in_mpa), .groups = "drop")

# #### Run MPA model ####
# model_linear_mpa <- feols(
#   Percent_Bleached ~ dhw + dhw:in_mpa | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data    = All_Bleaching_Events_Data_AllDHW_MPA
# )


# #### cluster sensitivity test ####
# model_linear_nocluster <- feols(
#   Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
#   data = All_Bleaching_Events_Data_AllDHW)

# model_linear_ecoregion <- feols(
#   Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
#   cluster = ~Ecoregion_Name, 
#   data = All_Bleaching_Events_Data_AllDHW)

# model_linear_site <- feols(
#   Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
#   cluster = ~Site_ID, 
#   data = All_Bleaching_Events_Data_AllDHW)


# #### lag sensitivity test ####
# model_lagminus3 <- feols(Percent_Bleached ~ dhw_lag.3 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
#                          data = All_Bleaching_Events_Data_AllDHW)
# model_lagminus2 <- feols(Percent_Bleached ~ dhw_lag.2 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
#                          data = All_Bleaching_Events_Data_AllDHW)
# model_lagminus1 <- feols(Percent_Bleached ~ dhw_lag.1 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
#                          data = All_Bleaching_Events_Data_AllDHW)
# model_lag0 <- feols(Percent_Bleached ~ dhw_lag0 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
#                     data = All_Bleaching_Events_Data_AllDHW)
# model_lag1 <- feols(Percent_Bleached ~ dhw_lag1 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
#                     data = All_Bleaching_Events_Data_AllDHW)
# model_lag2 <- feols(Percent_Bleached ~ dhw_lag2 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
#                     data = All_Bleaching_Events_Data_AllDHW)
# model_lag3 <- feols(Percent_Bleached ~ dhw_lag3 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
#                     data = All_Bleaching_Events_Data_AllDHW)
# model_lag4 <- feols(Percent_Bleached ~ dhw_lag4 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
#                     data = All_Bleaching_Events_Data_AllDHW)
# model_lag5 <- feols(Percent_Bleached ~ dhw_lag5 | Site_ID + Date_Year + Ecoregion_Month, cluster = ~Ecoregion_Name, 
#                     data = All_Bleaching_Events_Data_AllDHW)



# #### Quadratic DHW models ####
# model_quadratic <- feols(
#   Percent_Bleached ~ dhw + I(dhw^2) | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data    = All_Bleaching_Events_Data_AllDHW
# )

# model_quadratic_lat <- feols(
#   Percent_Bleached ~ dhw + I(dhw^2) + dhw:abs_lat + I(dhw^2):abs_lat | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data    = All_Bleaching_Events_Data_AllDHW
# )

# #### Cubic DHW models ####
# model_cubic <- feols(
#   Percent_Bleached ~ dhw + I(dhw^2) + I(dhw^3) | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data    = All_Bleaching_Events_Data_AllDHW
# )

# model_cubic_lat <- feols(
#   Percent_Bleached ~ dhw + I(dhw^2) + I(dhw^3) +
#     dhw:abs_lat + I(dhw^2):abs_lat + I(dhw^3):abs_lat | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data    = All_Bleaching_Events_Data_AllDHW
# )






# ################## conley errors ################
# #### Compute Conley (200km, spherical) vcov for each model ####
# conley_vcov <- function(model) {
#   vcov_conley(model,
#               lat      = "Latitude_Degrees",
#               lon      = "Longitude_Degrees",
#               cutoff   = 200,
#               distance = "spherical")
# }

# #### Different model types ####
# model_linear_vcov <- conley_vcov(model_linear)
# model_linear_lat_vcov <- conley_vcov(model_linear_lat)
# model_linear_turbidity_vcov <- conley_vcov(model_linear_turbidity)
# model_linear_depth_vcov <- conley_vcov(model_linear_depth)
# model_linear_dhw_std_vcov <- conley_vcov(model_linear_dhw_std)
# model_linear_hotspot_std_vcov <- conley_vcov(model_linear_hotspot_std)
# model_linear_hotspot_warming_vcov <- conley_vcov(model_linear_hotspot_warming)
# model_linear_lat_depth_vcov <- conley_vcov(model_linear_lat_depth)
# model_linear_lat_hotspot_vcov <- conley_vcov(model_linear_lat_hotspot)

# All_Bleaching_Events_Data_AllDHW_MPA <- All_Bleaching_Events_Data_AllDHW_MPA %>%
#   left_join(
#     All_Bleaching_Events_Data_AllDHW %>% 
#       select(Site_ID, Date_Year, Latitude_Degrees, Longitude_Degrees),
#     by = c("Site_ID", "Date_Year")
#   )
# model_linear_mpa_vcov <- conley_vcov(model_linear_mpa)


# #### Lag models ####
# model_lagminus3_vcov <- conley_vcov(model_lagminus3)
# model_lagminus2_vcov <- conley_vcov(model_lagminus2)
# model_lagminus1_vcov <- conley_vcov(model_lagminus1)
# model_lag0_vcov <- conley_vcov(model_lag0)
# model_lag1_vcov <- conley_vcov(model_lag1)
# model_lag2_vcov <- conley_vcov(model_lag2)
# model_lag3_vcov <- conley_vcov(model_lag3)
# model_lag4_vcov <- conley_vcov(model_lag4)
# model_lag5_vcov <- conley_vcov(model_lag5)


# #### Quadratic models ####
# model_quadratic_vcov <- conley_vcov(model_quadratic)
# model_quadratic_lat_vcov <- conley_vcov(model_quadratic_lat)

# #### Cubic models ####
# model_cubic_vcov <- conley_vcov(model_cubic)
# model_cubic_lat_vcov <- conley_vcov(model_cubic_lat)



# #### Custom styling function ####
# nature_table_style <- function(ft) {
#   ft %>%
#     font(fontname = "Arial", part = "all") %>%
#     fontsize(size = 7, part = "all") %>%
#     align(align = "center", part = "header") %>%
#     align(j = 1, align = "left", part = "body") %>%
#     align(j = 2:ncol_keys(ft), align = "center", part = "body") %>%
#     bold(part = "header") %>%
#     border_remove() %>%
#     hline_top(border = fp_border(width = 1.5), part = "all") %>%
#     hline_bottom(border = fp_border(width = 1.5), part = "all") %>%
#     hline(i = 1, border = fp_border(width = 1), part = "header") %>%
#     autofit() %>%
#     width(width = 1.2, unit = "in")
# }

# #### Interaction models ####
# modelsummary(
#   list(
#     "Latitude"      = model_linear_lat,
#     "Turbidity"     = model_linear_turbidity,
#     "Depth"         = model_linear_depth,
#     "DHW var."      = model_linear_dhw_std,
#     "SST var."  = model_linear_hotspot_std,
#     "SST warm." = model_linear_hotspot_warming,
#     #"Lat+Depth"     = model_linear_lat_depth,
#     "Lat+SST var."   = model_linear_lat_hotspot,
#     "MPA"           = model_linear_mpa
#   ),
#   vcov = list(      # Pass vcov objects separately here
#     model_linear_lat_vcov,
#     model_linear_turbidity_vcov,
#     model_linear_depth_vcov,
#     model_linear_dhw_std_vcov,
#     model_linear_hotspot_std_vcov,
#     model_linear_hotspot_warming_vcov,
#     #model_linear_lat_depth_vcov,
#     model_linear_lat_hotspot_vcov,
#     model_linear_mpa_vcov
#   ),
#   stars     = c('*' = 0.05, '**' = 0.01),
#   fmt       = 2,
#   statistic = "({std.error})",
#   coef_omit = "Intercept",
#   coef_rename = c(
#     "dhw"                       = "DHW",
#     "dhw:abs_lat"               = "DHW × Latitude",
#     "dhw:Turbidity"             = "DHW × Turbidity",
#     "dhw:Depth_m"               = "DHW × Depth",
#     "dhw:DHW_1985.2005_std"     = "DHW × DHW var",
#     "dhw:hotspot_1985.2005_std" = "DHW × SST var",
#     "dhw:hotspot_warming"       = "DHW × SST warm",
#     "dhw:in_mpaTRUE"            = "DHW × MPA"
#   ),
#   gof_omit  = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
#   notes     = list(
#     "Standard errors in parentheses, using Conley 200km errors.",
#     "All models include site and year fixed effects."
#   ),
#   title  = "Interaction models (DHW with site characteristics)",
#   output = "flextable"
# ) %>%
#   nature_table_style() %>%
#   save_as_docx(path = "table_interactions_nature.docx")


# #### Lag sensitivity ####
# modelsummary(
#   list(
#     "−3" = model_lagminus3, "−2" = model_lagminus2, "−1" = model_lagminus1,
#     "0"  = model_lag0,
#     "+1" = model_lag1, "+2" = model_lag2, "+3" = model_lag3,
#     "+4" = model_lag4, "+5" = model_lag5
#   ),
#   vcov = list(
#     model_lagminus3_vcov, model_lagminus2_vcov, model_lagminus1_vcov,
#     model_lag0_vcov,
#     model_lag1_vcov, model_lag2_vcov, model_lag3_vcov,
#     model_lag4_vcov, model_lag5_vcov
#   ),
#   stars     = c('*' = 0.05, '**' = 0.01),
#   fmt       = 2,
#   statistic = "({std.error})",
#   coef_omit = "Intercept",
#   coef_rename = c(
#     "dhw_lag.3" = "DHW", "dhw_lag.2" = "DHW", "dhw_lag.1" = "DHW",
#     "dhw_lag0"  = "DHW",
#     "dhw_lag1"  = "DHW", "dhw_lag2"  = "DHW", "dhw_lag3"  = "DHW",
#     "dhw_lag4"  = "DHW", "dhw_lag5"  = "DHW"
#   ),
#   gof_omit  = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
#   notes     = list(
#     "Standard errors in parentheses, using Conley 200km errors.",
#     "All models include site and year fixed effects.",
#     "Column headers show lag in months relative to bleaching observation."
#   ),
#   title  = "Temporal lag sensitivity analysis",
#   output = "flextable"
# ) %>%
#   nature_table_style() %>%
#   save_as_docx(path = "table_lags_nature.docx")

# #### Polynomial sensitivity ####
# modelsummary(
#   list(
#     "Linear"        = model_linear_lat,
#     "Quadratic"     = model_quadratic,
#     "Quadratic+Lat" = model_quadratic_lat,
#     "Cubic"         = model_cubic,
#     "Cubic+Lat"     = model_cubic_lat
#   ),
#   vcov = list(
#     model_linear_lat_vcov,
#     model_quadratic_vcov,
#     model_quadratic_lat_vcov,
#     model_cubic_vcov,
#     model_cubic_lat_vcov
#   ),
#   stars     = c('*' = 0.05, '**' = 0.01),
#   fmt       = 2,
#   statistic = "({std.error})",
#   coef_omit = "Intercept",
#   coef_rename = c(
#     "dhw"              = "DHW",
#     "I(dhw^2)"         = "DHW²",
#     "I(dhw^3)"         = "DHW³",
#     "dhw:abs_lat"      = "DHW × Lat",
#     "I(dhw^2):abs_lat" = "DHW² × Lat",
#     "I(dhw^3):abs_lat" = "DHW³ × Lat"
#   ),
#   gof_omit = "R2$|R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
#   notes    = list(
#     "* p < 0.05, ** p < 0.01",
#     "Standard errors in parentheses, using Conley 200km errors.",
#     "All models include site and year fixed effects."
#   ),
#   title  = "Polynomial DHW sensitivity analysis",
#   output = "flextable"
# ) %>%
#   nature_table_style() %>%
#   save_as_docx(path = "table_polynomial_sensitivity.docx")






# model_linear_nocluster <- feols(
#   Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
#   data = All_Bleaching_Events_Data_AllDHW)

# model_linear_ecoregion <- feols(
#   Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
#   cluster = ~Ecoregion_Name, 
#   data = All_Bleaching_Events_Data_AllDHW)

# model_linear_site <- feols(
#   Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month, 
#   cluster = ~Site_ID, 
#   data = All_Bleaching_Events_Data_AllDHW)

# model_linear_nocluster_vcov <- conley_vcov(model_linear_nocluster)
# model_linear_ecoregion_vcov <- conley_vcov(model_linear_ecoregion)
# model_linear_site_vcov <- conley_vcov(model_linear_site)


# modelsummary(
#   list(
#     "Conley 200km"  = model_linear, 
#     "Ecoregion" = model_linear_ecoregion,
#     "Site"      = model_linear_site,
#     "None"      = model_linear_nocluster
    
#   ),
#   vcov = list(
#     model_linear_vcov,
#     ~Ecoregion_Name,           # use model's own ecoregion clustering
#     ~Site_ID,                  # use model's own site clustering
#     "iid"                      # no clustering
#   ),
#   stars     = c('*' = 0.05, '**' = 0.01),
#   fmt       = 2,
#   statistic = "({std.error})",
#   coef_omit = "Intercept",
#   coef_rename = c("dhw" = "DHW"),
#   gof_omit  = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
#   notes     = list(
#     "Standard errors in parentheses.",
#     "All models include site and year fixed effects.",
#     "Clustering level varies by column as indicated."
#   ),
#   title  = "Standard error clustering sensitivity",
#   output = "flextable"
# ) %>%
#   nature_table_style() %>%
#   save_as_docx(path = "table_clustering_nature.docx")






# #### seasonality #### 
# #### Model with seasonality FE ####
# model_linear_lat_season <- feols(
#   Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
#   cluster = ~Ecoregion_Name,
#   data    = All_Bleaching_Events_Data_AllDHW
# )


# #### Compare to baseline model ####
# model_linear_lat <- feols(
#   Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year,
#   cluster = ~Ecoregion_Name,
#   data    = All_Bleaching_Events_Data_AllDHW
# )


# model_linear_lat_vcov        <- conley_vcov(model_linear_lat)
# model_linear_lat_season_vcov <- conley_vcov(model_linear_lat_season)

# #### Comparison table ####
# modelsummary(
#   list(
#     "Baseline (Site + Year FE)"            = model_linear_lat,
#     "+ Seasonality (Ecoregion×Month FE)"   = model_linear_lat_season
#   ),
#   vcov = list(
#     model_linear_lat_vcov,
#     model_linear_lat_season_vcov
#   ),
#   stars     = c('*' = 0.05, '**' = 0.01),
#   fmt       = 2,
#   statistic = "({std.error})",
#   coef_omit = "Intercept",
#   coef_rename = c(
#     "dhw"         = "DHW",
#     "dhw:abs_lat" = "DHW × |Latitude|"
#   ),
#   gof_omit  = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
#   notes     = list(
#     "Standard errors in parentheses, Conley (200km) spatial correction.",
#     "Both models include Site_ID and Date_Year fixed effects.",
#     "Ecoregion×Month FE absorbs seasonal differences in bleaching timing between hemispheres."
#   ),
#   title  = "Sensitivity to seasonality (ecoregion-month fixed effects)",
#   output = "flextable"
# ) %>%
#   nature_table_style() %>%
#   save_as_docx(path = "table_seasonality_sensitivity.docx")









# #### DHW specification sensitivity (DHW vs DHW_adj) ####

# All_Bleaching_Events_Data_AllDHW_test <- All_Bleaching_Events_Data_AllDHW
# All_Bleaching_Events_Data_AllDHW_test <- left_join(All_Bleaching_Events_Data_AllDHW_test, DF[, c(1, 2, 3, 4, 5, 6, 18)])
# #### Fit base + adjusted DHW models (Conley SEs computed post-estimation) ####
# model_linear_dhwspec <- feols(
#   Percent_Bleached ~ dhw | Site_ID + Date_Year + Ecoregion_Month,
#   data = All_Bleaching_Events_Data_AllDHW_Max)
# model_linear_dhwspec_vcov <- conley_vcov(model_linear_dhwspec)

# model_linear_lat_dhwspec <- feols(
#   Percent_Bleached ~ dhw + dhw:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
#   data = All_Bleaching_Events_Data_AllDHW_Max)
# model_linear_lat_dhwspec_vcov <- conley_vcov(model_linear_lat_dhwspec)

# model_linear_adj <- feols(
#   Percent_Bleached ~ DHW_adj | Site_ID + Date_Year + Ecoregion_Month,
#   data = All_Bleaching_Events_Data_AllDHW_Max)
# model_linear_adj_vcov <- conley_vcov(model_linear_adj)

# model_linear_adj_lat <- feols(
#   Percent_Bleached ~ DHW_adj + DHW_adj:abs_lat | Site_ID + Date_Year + Ecoregion_Month,
#   data = All_Bleaching_Events_Data_AllDHW_Max)
# model_linear_adj_lat_vcov <- conley_vcov(model_linear_adj_lat)

# #### Table: DHW specification sensitivity ####
# modelsummary(
#   list(
#     "DHW"               = model_linear_dhwspec,
#     "DHW + Lat."        = model_linear_lat_dhwspec,
#     "DHW (adj.)"        = model_linear_adj,
#     "DHW (adj.) + Lat." = model_linear_adj_lat
#   ),
#   vcov = list(
#     model_linear_dhwspec_vcov,
#     model_linear_lat_dhwspec_vcov,
#     model_linear_adj_vcov,
#     model_linear_adj_lat_vcov
#   ),
#   stars     = c('*' = 0.05, '**' = 0.01),
#   fmt       = 2,
#   statistic = "({std.error})",
#   coef_omit = "Intercept",
#   coef_rename = c(
#     "DHW"             = "DHW",
#     "DHW_adj"         = "DHW (adjusted)",
#     "DHW:abs_lat"     = "DHW × |Latitude|",
#     "DHW_adj:abs_lat" = "DHW (adj.) × |Latitude|"
#   ),
#   gof_omit  = "R2 Adj|R2 Within Adj|AIC|RMSE|FE|Std",
#   notes     = list(
#     "Standard errors in parentheses, Conley (200km) spatial correction.",
#     "All models include site and year fixed effects.",
#     "DHW (adjusted) uses the thermal anomaly adjustment from Ainsworth et al. (2016)."
#   ),
#   title  = "Sensitivity to DHW specification",
#   output = "flextable"
# ) %>%
#   nature_table_style() %>%
#   save_as_docx(path = "table_sensitivity_dhw_spec.docx")
