# ------------------------------------------------------------------------------
# 
# 03_get_geometries
# 
# Fetch geometries for the selected catchment sample.
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
library(sf)

# Get geometries for selected GSIM catchments ----------------------------------

message("fetching catchment geometries...")

gsim_selected_data <- readRDS("output/temp_rds_gpkg/01_catchment_records.rds")

geom_files <- gsim_selected_data %>%
  pull(file) %>%
  str_split("/") %>%
  lapply(function (x) {tail(x, 1) %>%
      str_replace(".mon", ".shp") %>%
      paste0("Data/GSIM/GSIM_metadata/GSIM_catchments/", .)}) %>%
  unlist() %>%
  tolower()

gsim_missing_geoms <- geom_files[which(!file.exists(geom_files))] # suspect stations
geom_files <- geom_files[!geom_files %in% gsim_missing_geoms]

gsim_sf <- vector(mode = "list", length = length(geom_files))
for (i in 1:length(geom_files)) {
  gsim_sf[[i]] <- read_sf(geom_files[i]) %>%
    st_cast("MULTIPOLYGON")
}

gsim_geoms <- gsim_sf %>%
  bind_rows() %>%
  st_make_valid() %>%
  rename(id = FILENAME) %>%
  mutate(id = toupper(id),
         area_km2 = as.numeric(st_area(.) / 10^6)) %>%
  select(id, area_km2, geometry) %>%
  arrange(id)

st_write(gsim_geoms, "output/temp_rds_gpkg/03_catchment_geoms.gpkg", delete_dsn = TRUE)



