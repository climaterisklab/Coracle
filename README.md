# Coracle

This repository contains code necessary to reproduce results from Kotz et al 2026: "Anthropogenic contributions to global coral reef bleaching are overwhelming and inequitable".

Note that this repository contains the code necessary to:
1. Merge and clean site-level bleaching records from GCBD and MERMAID, and extract DHW for each observation.
2. Estimate the primary linear and FE-binned DHW models with Conley (200km) standard errors, and run associated sensitivity/robustness analyses.
3. All scripts necessary to run the pipeline on an HPC cluster.
4. Generate Figures 1 and 3, Supplementary Figures S1, S6, and S12 and, Extended Data Figure 6.

Other scripts, including those used to generate GMST counterfactuals using the FaIR simple climate model, DHW counterfactuals using pattern scaling and block-bootstrap uncertainty estimates, and plotting scripts for Figures 2 and 4, all other Extended Data Figures, and the remaining SI Figures, can be found at the other GitHub repository. https://github.com/maxkotz17/Coracle_maxkotz.git, maintained by maximilian.kotz@bsc.es.


## Data
All data used in this analysis draw on two coral bleaching survey databases and one climate data product:

- **Global Coral Bleaching Database (GCBD):** [Global Coral Bleaching Database SQLite v11_24_21](https://springernature.figshare.com/articles/dataset/Global_Coral_Bleaching_Database/17076287?file=31573421) — site, sample, cover, and bleaching prevalence records.
- **MERMAID:** [datamermaid.org](https://datamermaid.org/) — coral reef monitoring data, accessed via the `mermaidr` R package.
- **NOAA Coral Reef Watch Degree Heating Week (DHW):** [NOAA CRW DHW v3.1](https://coralreefwatch.noaa.gov/product/5km/index_5km_dhw.php) — monthly 5km-resolution thermal stress data, 1985–2025.

Site-level bleaching records from GCBD and MERMAID are deduplicated and merged into a single site-year panel, with DHW extracted per observation.

Other data referenced in the scripts can be found here: https://www.dropbox.com/scl/fo/otfp6bzipjcjoc8iafyhl/ADucfteBxwoVWzrd7e8emHU?rlkey=2o35hszvnkwctdo5lc1jg3z79&st=0adhkm37&dl=0. 

## Code

All R scripts to execute the data merging, modelling, and visualisation are stored in the `main` directory:

```
Pipeline
├── A - Data preparation
│   └── 01_Merging_and_Cleaning_GCBD_MERMAID.R                   # merge GCBD + MERMAID into a site-year panel; extract DHW
├── B - Model estimation
│   └── 02_VCOV_Conley_Errors.R                     # primary linear FE + FE-binned DHW models with Conley (200km) SEs
├── C - Sensitivity analyses
│   └── 03_Sensitivity_Analyses.R                   # robustness checks: model type, clustering, lags, polynomial DHW,
│                                                 #   MPA interactions, seasonality FE, DHW specification
├── D - HPC pipeline
│   ├── 04_HPC_Scripts/setup_linear_lat_conley_season.R         # prepares inputs/environment for the batch run
│   ├── 04_HPC_Scripts/mc_batch_linear_lat_conley_season.R      # Monte Carlo batch estimation script run on each HPC node
│   ├── 04_HPC_Scripts/slurm_setup_linear_lat_conley_season.sh  # SLURM submission script for the setup step
│   └── 04_HPC_Scripts/slurm_mc_linear_lat_conley_season.sh     # SLURM submission script for the batch MC run
│
├── E - Figure generation for main text
│   ├── 05_Figure1.R                                # emissions/GMT time series + DHW response + coefficient panels
│   └── 06_Figure3.R                                # coral bleaching source attribution bar charts
└── F - Figure generation for supplement
    ├── 07_FigureSupp1.R                                 # attribution panels by country/company (Fig. S6)
    ├── 08_FigureSupp6.R                              # model response curves with FE-binned estimates (Fig. S1)
    └── 09_FigureSupp12.R                             # global map of panel bleaching survey sites (Fig. S12)
```
    
## Statistical framework

The primary impact model is a linear fixed-effects specification (`feols`, `fixest` package):

```
Percent_Bleached ~ DHW + DHW:abs_lat | Site_ID + Date_Year + Ecoregion_Month
```

with Conley (200km, spherical) standard errors clustered by ecoregion.

## How to replicate results

1. Clone this repository.
2. Obtain the GCBD SQLite database and NOAA CRW DHW `.nc` files (see Data section above), and a MERMAID account/API access for `mermaidr`.
3. Edit the file paths at the top of each script (currently set to local Dropbox paths) to point to your own data directory.
4. Run scripts in the order shown in the Code section: A. Data preparation → B. Model estimation → C. Sensitivity analyses → Scripts 01-08 of Max's repository → D. HPC pipeline → E and F. Figure generation.

## Use of code and data

This code is provided to accompany the manuscript for review and reproducibility purposes.
