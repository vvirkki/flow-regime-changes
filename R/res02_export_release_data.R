# ------------------------------------------------------------------------------
# 
# res02_export_release_data
# 
# Export csv files of main results for supporting data Zenodo release.
# 
# Analysis script for publication:
# 
# "Archetypal flow regime change classes and their associations with
# anthropogenic drivers of global streamflow alterations"
#
# Vili Virkki, Reetik Kumar Sahu, Mikhail Smilovic, Josias Láng-Ritter,
# Miina Porkka, Matti Kummu
#
# Published in Environmental Research Communications
# https://doi.org/10.1088/2515-7620/ad9439
# 
# Data availability:
# https://doi.org/10.5281/zenodo.11102422
# 
# Code availability:
# https://github.com/vvirkki/flow-regime-changes
# 
# Corresponding authors of the article:
# Vili Virkki (vili.virkki@uef.fi)
# Matti Kummu (matti.kummu@aalto.fi)
# 
# Script author:
# Vili Virkki (vili.virkki@uef.fi)
#
# Date:
# 14.11.2024
# 
# ------------------------------------------------------------------------------

library(tidyverse)
library(sf)
options(dplyr.summarise.inform = FALSE)

out_dir <- "output/release_datafiles/"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# CATCHMENT GEOMETRIES ---------------------------------------------------------

geoms <- read_sf("output/figures_tables/fig3_frc_assignments.gpkg") %>%
  select(id) %>%
  st_write(paste0(out_dir, "catchment_geometries.gpkg"), delete_dsn = TRUE)

geom_meta <- read_sf("output/figures_tables/fig3_frc_assignments.gpkg") %>%
  st_drop_geometry() %>%
  select(id, area_km2, parents, children)

# CATCHMENT ATTRIBUTES ---------------------------------------------------------

ctm_records <- readRDS("output/temp_rds_gpkg/01_catchment_records.rds") %>%
  select(id, yr_begin, yr_end, missing_obs_fraction)

ctm_attributes <- readRDS("output/temp_rds_gpkg/02_catchment_attributes.rds") %>%
  select(id, river, station, country, reported_area_km2, GSIM_area_km2, GSIM_geom_qflag)

ctm_export <- ctm_attributes %>%
  left_join(ctm_records, by = "id") %>%
  left_join(geom_meta, by = "id") %>%
  rename(GSIM_estimated_area_km2 = GSIM_area_km2,
         GSIM_shapefile_area_km2 = area_km2,
         streamflow_data_begin = yr_begin,
         streamflow_data_end = yr_end,
         streamflow_missing_obs_fraction = missing_obs_fraction,
         catchment_is_subcatchment_of = parents,
         catchment_has_subcatchments = children) %>%
  select(id, river, station, country, ends_with("km2"), GSIM_geom_qflag, everything())

write_csv(ctm_export, paste0(out_dir, "catchment_attributes.csv"))

# STREAMFLOW TRENDS ------------------------------------------------------------

streamflow_trends <- read_sf("output/figures_tables/fig2_streamflow_trend_components.gpkg") %>%
  st_drop_geometry() %>%
  mutate(trend_unit = "mm_per_month",
         trend_method = "theil_sen_slope") %>%
  select(id, trend_method, trend_unit, mean, sd, quantile05, quantile95) %>%
  rename_with(~ paste0("streamflow_trend_", .x), -c(id, trend_method, trend_unit))

write_csv(streamflow_trends, paste0(out_dir, "streamflow_trends.csv"))

# FRC ASSIGNMENTS --------------------------------------------------------------

frc_assignments <- read_sf("output/figures_tables/fig3_frc_assignments.gpkg") %>%
  st_drop_geometry() %>%
  select(id, case) %>%
  mutate(case = str_replace_all(case, "aa_|ab_|ac_|ad_|ae_", "")) %>%
  rename(frc_assignment = case)

write_csv(frc_assignments, paste0(out_dir, "frc_assignments.csv"))

# DRIVER TRENDS AND CHANGES ----------------------------------------------------

driver_trends <- read_sf("output/figures_tables/fig4_driver_trend_components.gpkg") %>%
  st_drop_geometry() %>%
  select(id, matches("_mean|_sd")) %>%
  mutate(trend_unit = "mm_per_month",
         trend_method = "theil_sen_slope") %>%
  select(id, trend_method, trend_unit, everything()) %>%
  rename_with(~ paste0(.x, "_trend"), -c(id, trend_method, trend_unit))

write_csv(driver_trends, paste0(out_dir, "et_prec_wateruse_trends.csv"))

dor_change <- read_sf("output/figures_tables/fig4_driver_trend_components.gpkg") %>%
  st_drop_geometry() %>%
  select(id, degree_of_regulation_abs_change) %>%
  rename(dor_change = degree_of_regulation_abs_change) %>%
  mutate(dor_change = ifelse(dor_change > 0, "increase", "no_increase"))

write_csv(dor_change, paste0(out_dir, "dor_increases_nonincreases.csv"))

