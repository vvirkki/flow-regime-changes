# ------------------------------------------------------------------------------
# 
# 02_get_metadata
# 
# Fetch attribute and metadata for the selected catchment sample.
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
# LINK TO ADD
# 
# Data availability:
# https://doi.org/10.5281/zenodo.11102423
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

# GSIM attributes & metadata ---------------------------------------------------

message("fetching catchment attributes...")

gsim_selected_data <- readRDS("output/temp_rds_gpkg/01_catchment_records.rds")
catchment_attr <- read_csv("Data/GSIM/GSIM_metadata/GSIM_catalog/GSIM_catchment_characteristics.csv",
                           show_col_types = FALSE)

sel_catchment_attr <- catchment_attr %>%
  select(gsim.no, long.new, lat.new, area.meta, area.est, quality,
         climate.type, landcover.type, ele.mean, slp.mean) %>%
  setNames(c("id", "x", "y", "reported_area_km2", "GSIM_area_km2", "GSIM_geom_qflag",
             "climateZone", "landcoverClass", "meanElevation", "meanSlope"))

catchment_meta <- read_csv("Data/GSIM/GSIM_metadata/GSIM_catalog/GSIM_metadata.csv",
                           show_col_types = FALSE)

sel_catchment_meta <- catchment_meta %>%
  mutate(name = paste0(river, " river ", station, " station")) %>%
  select(gsim.no, name, river, station, country) %>%
  rename(id = gsim.no)

gsim_attr <- gsim_selected_data %>%
  select(id) %>%
  left_join(sel_catchment_meta, by = "id") %>%
  left_join(sel_catchment_attr, by = "id") %>%
  mutate(GSIM_geom_qflag = tolower(GSIM_geom_qflag))

saveRDS(gsim_attr, "output/temp_rds_gpkg/02_catchment_attributes.rds")


