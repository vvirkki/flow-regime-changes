# ------------------------------------------------------------------------------
# 
# 07_extract_drivers
# 
# Extract zonal statistics of all drivers within all catchments of the selected
# catchment sample.
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
library(terra)
library(sf)
library(exactextractr)
options(dplyr.summarise.inform = FALSE)

# Extract functions ------------------------------------------------------------

.nmfun <- function(values, ...) {return (values)}

.ee_area_weighted_mean <- function(sel_rast, sel_geom,
                                   yr_begin = NULL, yr_end = NULL) {
  
  if (!is.null(yr_begin) & !is.null(yr_end)) {
    sel_yrs <- yr_begin:yr_end %>% paste(collapse = "|")
    sel_rast <- sel_rast[[grepl(sel_yrs, names(sel_rast))]]
  }
  ret <- exact_extract(sel_rast,
                       sel_geom,
                       fun = "weighted_mean",
                       weights = "area",
                       append_cols = "id",
                       colname_fun = .nmfun,
                       force_df = TRUE,
                       max_cells_in_memory = 2^31-1)
  return (ret)
  
}

.ee_to_long <- function(ee, streamflow_period) {

  ee_long <- ee %>%
    pivot_longer(cols = -id,
                 names_to = "variable") %>%
    separate_wider_delim(variable, delim = "_", names = c("variable", "timestep")) %>%
    separate_wider_delim(timestep, delim = ".", names = c("year", "month")) %>%
    mutate(year = as.numeric(year),
           month = as.numeric(month),
           timestep = make_date(year, month)) %>%
    left_join(streamflow_period, by = "id") %>%
    filter(timestep >= min_timestep & timestep <= max_timestep) %>%
    select(-ends_with("timestep"))

  return (ee_long)

}


# Read data --------------------------------------------------------------------

params <- readRDS("output/temp_rds_gpkg/params.rds")
geoms <- read_sf("output/temp_rds_gpkg/03_catchment_geoms.gpkg")

streamflow_period <- readRDS("output/temp_rds_gpkg/05_monthly_streamflow.rds") %>%
  mutate(timestep = make_date(year, month)) %>%
  group_by(id) %>%
  summarise(min_timestep = min(timestep),
            max_timestep = max(timestep))

total_prec <- list.files("output/temp_drivers_tif_gpkg/total-prec", full.names = TRUE) %>%
  rast()
e_total <- list.files("output/temp_drivers_tif_gpkg/e-total", full.names = TRUE) %>%
  rast()
wateruse <- list.files("output/temp_drivers_tif_gpkg/wateruse", full.names = TRUE)

# Extract ----------------------------------------------------------------------

message("extracting zonal statistics...")

# precipitation, evaporation: mm/month -----------------------------------------

tbl_total_prec <- .ee_area_weighted_mean(total_prec, geoms,
                                         yr_begin = params$yr_begin,
                                         yr_end = params$yr_end) %>%
  .ee_to_long(streamflow_period)

saveRDS(tbl_total_prec, "output/temp_rds_gpkg/07_total_prec.rds")
message(paste0("saved total_prec at ", Sys.time()))

tbl_e_total <- .ee_area_weighted_mean(e_total, geoms,
                                      yr_begin = params$yr_begin,
                                      yr_end = params$yr_end) %>%
  .ee_to_long(streamflow_period)

saveRDS(tbl_e_total, "output/temp_rds_gpkg/07_e_total.rds")
message(paste0("saved e_total at ", Sys.time()))

# water use: mm/month ----------------------------------------------------------

tbl_wateruse <- wateruse %>%
  lapply(rast) %>%
  lapply(.ee_area_weighted_mean, sel_geom = geoms) %>%
  lapply(.ee_to_long, streamflow_period = streamflow_period) %>%
  bind_rows()

wu_all <- tbl_wateruse %>%
  group_by(id, year, month) %>%
  summarise(value = sum(value)) %>%
  mutate(variable = "totalWu") %>%
  select(id, variable, year, month, value) %>%
  ungroup()

saveRDS(wu_all, "output/temp_rds_gpkg/07_wateruse_all_sectors.rds")
message(paste0("saved wateruse at ", Sys.time()))

# dams: point-in-polygon -------------------------------------------------------

dams <- read_sf("output/temp_drivers_tif_gpkg/dams.gpkg")
message("processing dams...")

dam_years <- vector(mode = "list", length = length(unique(geoms$id)))
for (i in 1:nrow(geoms)) {

  if (i %% 500 == 0) { message(paste0("dams in ", i, " catchments processed...")) }

  dams_df <- dams[which(st_within(dams, geoms[i,], sparse = FALSE)),]

  duplicate_geometries <- dams_df %>%
    group_by(geom) %>%
    summarise(nn = n(), id = paste(id, collapse = "|")) %>%
    filter(nn > 1)
  
  duplicate_name_year <- dams_df %>%
    group_by(name, year) %>%
    summarise(nn = n(), id = paste(id, collapse = "|")) %>%
    filter(nn > 1)
  
  dams_df_out <- dams_df %>%
    st_drop_geometry() %>%
    arrange(year) %>%
    mutate(reservoir_capacity_cumulative = cumsum(reservoir_capacity),
           ctm_id = geoms[i,]$id) %>%
    rename(dam_id = id) %>%
    select(ctm_id, dam_id, name, country, river, year, starts_with("reservoir"),
           everything())

  dam_years[[i]] <- dams_df_out

}

dams_captured <- dam_years %>%
  bind_rows() %>%
  select(dam_id) %>%
  distinct() %>%
  nrow()

message(paste0(dams_captured, " dams captured within all catchments..."))

ctm_areas <- readRDS("output/temp_rds_gpkg/02_catchment_attributes.rds") %>%
  mutate(GSIM_area_m2 = GSIM_area_km2 * 10^6) %>%
  select(id, GSIM_area_m2)

maf <- readRDS("output/temp_rds_gpkg/05_monthly_streamflow.rds") %>%
  group_by(id, year) %>%
  summarise(annual_streamflow_m = sum(streamflow) / 1000) %>%
  summarise(annual_streamflow_mean_m = mean(annual_streamflow_m)) %>%
  ungroup() %>%
  left_join(ctm_areas, by = "id") %>%
  mutate(maf_m3 = annual_streamflow_mean_m * GSIM_area_m2) %>%
  select(id, maf_m3)

dor <- dam_years %>%
  bind_rows() %>%
  left_join(maf, by = c("ctm_id" = "id")) %>%
  mutate(degree_of_regulation = ((reservoir_capacity_cumulative * 10^3) / maf_m3) * 100) %>%
  select(ctm_id, year, degree_of_regulation) %>%
  rename(id = ctm_id,
         value = degree_of_regulation) %>%
  group_by(id, year) %>%
  summarise(value = max(value)) %>%
  ungroup()

dor_long <- expand_grid(id = unique(streamflow_period$id), variable = "degree_of_regulation",
                        year = 1971:2010, month = 1:12) %>%
  left_join(dor, by = c("year", "id")) %>%
  arrange(id, variable, year, month) %>%
  group_by(id) %>%
  fill(value, .direction = "down") %>%
  ungroup() %>%
  mutate(value = ifelse(is.na(value), 0, value),
         timestep = make_date(year, month)) %>%
  left_join(streamflow_period, by = "id") %>%
  filter(timestep >= min_timestep & timestep <= max_timestep) %>%
  select(-ends_with("timestep"))

saveRDS(dor_long, "output/temp_rds_gpkg/07_dor.rds")



