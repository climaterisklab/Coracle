# Calculate Monthly Maximum DHW from Daily NOAA Data
# This script downloads daily DHW data and calculates monthly maximum values
# to extend the NOAA monthly aggregated DHW dataset beyond October 2023

# Required packages
library(ncdf4)
library(httr)
library(lubridate)

# Function to download a file if it doesn't exist
download_if_needed <- function(url, destfile) {
  if (!file.exists(destfile)) {
    cat("Downloading:", basename(destfile), "\n")
    tryCatch({
      response <- GET(url, write_disk(destfile, overwrite = TRUE), progress())
      
      # Check if download was successful
      if (status_code(response) == 200) {
        # Verify file size is reasonable (> 1KB)
        file_size <- file.info(destfile)$size
        if (is.na(file_size) || file_size < 1000) {
          cat("  WARNING: Downloaded file is too small or empty\n")
          file.remove(destfile)
          return(FALSE)
        }
        cat("  Success! File size:", round(file_size/1024/1024, 2), "MB\n")
        return(TRUE)
      } else {
        cat("  ERROR: HTTP status", status_code(response), "\n")
        if (file.exists(destfile)) file.remove(destfile)
        return(FALSE)
      }
    }, error = function(e) {
      cat("  ERROR downloading:", e$message, "\n")
      if (file.exists(destfile)) file.remove(destfile)
      return(FALSE)
    })
  } else {
    # Check existing file is valid
    file_size <- file.info(destfile)$size
    if (is.na(file_size) || file_size < 1000) {
      cat("Existing file is invalid, re-downloading:", basename(destfile), "\n")
      file.remove(destfile)
      return(download_if_needed(url, destfile))
    }
    cat("File already exists:", basename(destfile), "(", round(file_size/1024/1024, 2), "MB )\n")
    return(TRUE)
  }
}

# Function to get all days in a month
get_days_in_month <- function(year, month) {
  start_date <- as.Date(paste0(year, "-", sprintf("%02d", month), "-01"))
  end_date <- ceiling_date(start_date, "month") - days(1)
  seq(start_date, end_date, by = "day")
}

# Function to download daily DHW files for a specific month
download_daily_dhw_month <- function(year, month, output_dir = "daily_dhw") {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get all dates in the month
  dates <- get_days_in_month(year, month)
  
  # Base URL for daily DHW files (includes year subdirectory)
  base_url <- paste0("https://www.ncei.noaa.gov/data/oceans/crw/5km/v3.1/nc/v1.0/daily/dhw/", year, "/")
  
  downloaded_files <- c()
  
  for (date in as.character(dates)) {
    date_str <- gsub("-", "", date)
    filename <- paste0("ct5km_dhw_v3.1_", date_str, ".nc")
    url <- paste0(base_url, filename)
    destfile <- file.path(output_dir, filename)
    
    if (download_if_needed(url, destfile)) {
      downloaded_files <- c(downloaded_files, destfile)
    }
  }
  
  return(downloaded_files)
}

# Function to calculate monthly maximum DHW from daily files
calculate_monthly_max_dhw <- function(daily_files, output_file) {
  cat("\nCalculating monthly maximum DHW from", length(daily_files), "daily files...\n")
  
  if (length(daily_files) == 0) {
    stop("No daily files provided")
  }
  
  # Validate files can be opened
  cat("Validating NetCDF files...\n")
  valid_files <- c()
  for (f in daily_files) {
    tryCatch({
      test_nc <- nc_open(f)
      nc_close(test_nc)
      valid_files <- c(valid_files, f)
    }, error = function(e) {
      cat("  WARNING: Cannot open", basename(f), "-", e$message, "\n")
    })
  }
  
  if (length(valid_files) == 0) {
    stop("No valid NetCDF files could be opened!")
  }
  
  cat("Valid files:", length(valid_files), "out of", length(daily_files), "\n")
  daily_files <- valid_files
  
  # Open first file to get dimensions and coordinates
  nc_first <- nc_open(daily_files[1])
  
  # Get dimensions
  lon <- ncvar_get(nc_first, "lon")
  lat <- ncvar_get(nc_first, "lat")
  
  # Get dimension info
  lon_dim <- nc_first$dim$lon
  lat_dim <- nc_first$dim$lat
  
  # Initialize array for maximum values
  dhw_max <- array(NA, dim = c(length(lon), length(lat)))
  
  nc_close(nc_first)
  
  # Process each daily file
  for (i in seq_along(daily_files)) {
    cat("Processing file", i, "of", length(daily_files), ":", basename(daily_files[i]), "\n")
    
    nc <- nc_open(daily_files[i])
    dhw_daily <- ncvar_get(nc, "degree_heating_week")
    nc_close(nc)
    
    # Update maximum (handling NA values)
    if (i == 1) {
      dhw_max <- dhw_daily
    } else {
      dhw_max <- pmax(dhw_max, dhw_daily, na.rm = TRUE)
    }
  }
  
  # Create output NetCDF file
  cat("\nCreating output file:", output_file, "\n")
  
  # Reopen first file to copy attributes
  nc_first <- nc_open(daily_files[1])
  
  # Define dimensions
  lon_dim_out <- ncdim_def("lon", "degrees_east", lon)
  lat_dim_out <- ncdim_def("lat", "degrees_north", lat)
  
  # Define variable
  dhw_var <- ncvar_def(
    name = "degree_heating_week",
    units = "degree_C weeks",
    dim = list(lon_dim_out, lat_dim_out),
    missval = -9999,
    longname = "Monthly Maximum Degree Heating Week",
    prec = "float"
  )
  
  # Create NetCDF file
  nc_out <- nc_create(output_file, dhw_var)
  
  # Write data
  ncvar_put(nc_out, dhw_var, dhw_max)
  
  # Add global attributes
  ncatt_put(nc_out, 0, "title", "Monthly Maximum Degree Heating Week")
  ncatt_put(nc_out, 0, "source", "Calculated from NOAA Coral Reef Watch daily DHW data")
  ncatt_put(nc_out, 0, "institution", "NOAA Coral Reef Watch")
  ncatt_put(nc_out, 0, "Conventions", "CF-1.6")
  ncatt_put(nc_out, 0, "creation_date", as.character(Sys.time()))
  
  nc_close(nc_out)
  nc_close(nc_first)
  
  cat("Output file created successfully!\n")
  
  return(dhw_max)
}

# Function to verify against official NOAA monthly max file
verify_calculation <- function(calculated_file, official_file) {
  cat("\n=== VERIFICATION ===\n")
  cat("Comparing calculated file with official NOAA file...\n")
  
  # Open both files
  nc_calc <- nc_open(calculated_file)
  nc_official <- nc_open(official_file)
  
  # Read DHW data
  dhw_calc <- ncvar_get(nc_calc, "degree_heating_week")
  dhw_official <- ncvar_get(nc_official, "degree_heating_week")
  
  nc_close(nc_calc)
  nc_close(nc_official)
  
  # Compare dimensions
  if (!all(dim(dhw_calc) == dim(dhw_official))) {
    cat("WARNING: Dimensions don't match!\n")
    cat("Calculated:", dim(dhw_calc), "\n")
    cat("Official:", dim(dhw_official), "\n")
    return(FALSE)
  }
  
  # Calculate differences
  diff <- dhw_calc - dhw_official
  
  # Statistics (ignoring NA values)
  valid_mask <- !is.na(dhw_calc) & !is.na(dhw_official)
  diff_valid <- diff[valid_mask]
  
  cat("\nComparison Statistics:\n")
  cat("  Total pixels:", length(dhw_calc), "\n")
  cat("  Valid pixels in both:", sum(valid_mask), "\n")
  cat("  Mean difference:", mean(diff_valid, na.rm = TRUE), "\n")
  cat("  Max absolute difference:", max(abs(diff_valid), na.rm = TRUE), "\n")
  cat("  RMS difference:", sqrt(mean(diff_valid^2, na.rm = TRUE)), "\n")
  
  # Check if values are essentially identical
  max_diff <- max(abs(diff_valid), na.rm = TRUE)
  if (max_diff < 0.01) {
    cat("\n✓ VERIFICATION PASSED: Values match official NOAA data!\n")
    return(TRUE)
  } else {
    cat("\n✗ VERIFICATION FAILED: Significant differences found.\n")
    cat("  Check your calculation methodology.\n")
    return(FALSE)
  }
}

# Main workflow function
process_month <- function(year, month, verify = FALSE) {
  cat("\n=================================================\n")
  cat("Processing:", year, "-", sprintf("%02d", month), "\n")
  cat("=================================================\n")
  
  # Download daily files
  daily_files <- download_daily_dhw_month(year, month)
  
  if (length(daily_files) == 0) {
    cat("No files downloaded. Exiting.\n")
    return(NULL)
  }
  
  # Calculate monthly maximum
  output_file <- sprintf("ct5km_dhw-max_v3.1_%04d%02d.nc", year, month)
  dhw_max <- calculate_monthly_max_dhw(daily_files, output_file)
  
  # Verification (optional)
  if (verify) {
    # Download official file for comparison
    official_url <- sprintf(
      "https://www.ncei.noaa.gov/data/oceans/crw/5km/v3.1/nc/v1.0/monthly/%04d/ct5km_dhw-max_v3.1_%04d%02d.nc",
      year, year, month
    )
    official_file <- sprintf("official_ct5km_dhw-max_v3.1_%04d%02d.nc", year, month)
    
    if (download_if_needed(official_url, official_file)) {
      verify_calculation(output_file, official_file)
    } else {
      cat("Could not download official file for verification.\n")
    }
  }
  
  cat("\nDone! Output saved to:", output_file, "\n")
  return(output_file)
}

# ========================================
# EXAMPLE USAGE
# ========================================

# Example 1: Verify your calculation method using October 2023 (last available official data)
cat("EXAMPLE 1: Verification using October 2023\n")
cat("This will test your methodology against official NOAA data\n\n")

# Uncomment to run verification:
process_month(2025, 12, verify = TRUE)

# Example 2: Calculate monthly max for November 2023 onwards
cat("\nEXAMPLE 2: Calculate for months after October 2023\n\n")

# Uncomment to process November 2023:
# process_month(2023, 11, verify = FALSE)

# Example 3: Process multiple months
cat("\nEXAMPLE 3: Process multiple months (Nov 2023 - Dec 2024)\n\n")

# Uncomment to process all months from Nov 2023 to Dec 2024:
for (year in 2025:2025) {
  start_month <- ifelse(year == 2024, 11, 1)
  end_month <- ifelse(year == 2025, 1, 1)

  for (month in start_month:end_month) {
    process_month(year, month, verify = FALSE)
  }
}

cat("\n=================================================\n")
cat("Script loaded successfully!\n")
cat("=================================================\n")
cat("\nTo use this script:\n")
cat("1. First verify your method: process_month(2023, 10, verify = TRUE)\n")
cat("2. Then calculate new months: process_month(2023, 11, verify = FALSE)\n")
cat("3. Or batch process: see Example 3 in the script\n")