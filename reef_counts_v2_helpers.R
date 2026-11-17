#### helpers for the v2 all-reef-location attribution (Figure 3, Tables S7-S10)
#### The BSC output (/esarchive/Coracle/bleaching_attribution_v2/reef_locations/, 106 GB) is reduced on
#### the cluster by 04_HPC_Scripts/reef_counts_v2.py to reef_counts_v2.nc: per country (cells.csv `admin`),
#### threshold, year, scenario and sample, the number of factual bleaching events (cells > thr) and the
#### number of those same cells that also exceed thr in the counterfactual. Counts add over cells, so
#### any group of countries is a sum over countries.
library(ncdf4)

# cellset: "reefs" = all modelled reef cells (~72,500); "sites" = only the 4,698 cells holding survey sites
read_reef_counts <- function(path, cellset = c("reefs", "sites")) {
  cellset <- match.arg(cellset)
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  sets <- as.vector(ncvar_get(nc, "cellset"))
  k <- match(cellset, sets)
  if (is.na(k)) stop("cellset ", cellset, " not in ", path)
  # ncdf4 reverses the (C-order) dimensions:
  # factual_count -> [sample, year, country, cellset, threshold]
  # cf_count      -> [sample, scenario, year, country, cellset, threshold]
  # n_cells       -> [country, cellset]
  fact <- ncvar_get(nc, "factual_count", collapse_degen = FALSE)
  cf   <- ncvar_get(nc, "cf_count", collapse_degen = FALSE)
  nce  <- ncvar_get(nc, "n_cells", collapse_degen = FALSE)
  country <- as.vector(ncvar_get(nc, "country"))
  list(
    cellset   = cellset,
    # drop the (now length-1) cellset dimension; removing a length-1 dim keeps the element order
    factual   = array(fact[, , , k, , drop = FALSE], dim = dim(fact)[-4]),
    cf        = array(cf[, , , , k, , drop = FALSE], dim = dim(cf)[-5]),
    n_cells   = setNames(as.vector(nce[, k]), country),
    threshold = as.vector(ncvar_get(nc, "threshold")),
    country   = country,
    year      = as.vector(ncvar_get(nc, "year")),
    scenario  = as.vector(ncvar_get(nc, "scenario")),
    sample    = as.vector(ncvar_get(nc, "sample"))
  )
}

#### Country groups for the table rows, in cells.csv `admin` names ####
# AOSIS: all 39 members (Guinea-Bissau, Guyana and Suriname have no reef cells, so they add nothing)
aosis_members <- c(
  "Antigua and Barbuda", "The Bahamas", "Barbados", "Belize", "Cape Verde", "Comoros", "Cook Islands",
  "Cuba", "Dominica", "Dominican Republic", "Fiji", "Federated States of Micronesia", "Grenada",
  "Guinea-Bissau", "Guyana", "Haiti", "Jamaica", "Kiribati", "Maldives", "Marshall Islands", "Mauritius",
  "Nauru", "Niue", "Palau", "Papua New Guinea", "Samoa", "Sao Tome and Principe", "Seychelles",
  "Singapore", "Saint Kitts and Nevis", "Saint Lucia", "Saint Vincent and the Grenadines",
  "Solomon Islands", "Suriname", "East Timor", "Tonga", "Trinidad and Tobago", "Tuvalu", "Vanuatu"
)

reef_country_groups <- list(
  AOSIS     = aosis_members,
  Indonesia = "Indonesia",
  Japan     = "Japan",
  # survey data labelled Hong Kong sites as China; Taiwan is kept separate as before
  China     = c("China", "Hong Kong S.A.R.", "Macau S.A.R"),
  # USA including its territories (the survey-site version mixed "United States" with PR, Guam, USVI, USMOI)
  USA       = c("United States of America", "Puerto Rico", "United States Virgin Islands", "Guam",
                "Northern Mariana Islands", "American Samoa", "United States Minor Outlying Islands",
                "US Naval Base Guantanamo Bay"),
  # as in the survey-site version, external territories count as Australia
  Australia = c("Australia", "Coral Sea Islands", "Indian Ocean Territories", "Ashmore and Cartier Islands")
)

#### Sum counts over a set of countries -> same structure as get_batch_counts_corrected_v2() returned ####
# countries = NULL means Global (all reef cells)
group_counts <- function(rc, threshold, countries = NULL) {
  t_idx <- match(threshold, rc$threshold)
  if (is.na(t_idx)) stop("threshold ", threshold, " not in counts file (", paste(rc$threshold, collapse = ", "), ")")
  c_idx <- if (is.null(countries)) seq_along(rc$country) else which(rc$country %in% countries)
  missing <- setdiff(countries, rc$country)
  if (length(missing)) message("  no reef cells for: ", paste(missing, collapse = ", "))
  if (length(c_idx) == 0) return(NULL)

  f  <- rc$factual[, , c_idx, t_idx, drop = FALSE]           # sample, year, country
  cf <- rc$cf[, , , c_idx, t_idx, drop = FALSE]              # sample, scenario, year, country
  f  <- apply(f, c(2, 1), sum)                               # year, sample
  cf <- aperm(apply(cf, c(1, 2, 3), sum), c(3, 1, 2))        # year, sample, scenario

  list(
    factual_count  = f,
    cf_count       = cf,
    years          = rc$year,
    scenario_names = rc$scenario,
    threshold      = threshold,
    n_cells        = sum(rc$n_cells[c_idx])
  )
}
