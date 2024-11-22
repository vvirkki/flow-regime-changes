#!/bin/bash

# query catchments based on streamflow data availability
Rscript R/01_select_catchments.R

# get attribute information for selected catchments
Rscript R/02_get_metadata.R

# get geometries of selected catchments
Rscript R/03_get_geometries.R

# resolve overlapping and duplicate catchments
Rscript R/04_resolve_geometries.R

# preprocess and fill missing streamflow in selected catchments
Rscript R/05_preprocess_streamflow.R

# prepare drivers
Rscript R/06_preprocess_drivers.R

# extract driver values in selected catchments
Rscript R/07_extract_drivers.R

# build trends
Rscript R/08_determine_trends.R

# produce output results
Rscript R/res01_assign_analyse_frcs.R
Rscript R/res02_export_release_data.R
