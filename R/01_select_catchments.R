# ------------------------------------------------------------------------------
# 
# 01_select_catchments
# 
# Query the GSIM database and select catchments that are large enough and have
# a long enough streamflow record.
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
options(dplyr.summarise.inform = FALSE)

# Parameters -------------------------------------------------------------------

params_in <- readLines("params") %>%
  lapply(str_split, ";") %>%
  lapply(unlist)
params <- params_in %>%
  map(2) %>%
  lapply(as.numeric) %>%
  setNames(params_in %>% map(1))

if (!dir.exists("output/temp_rds_gpkg")) {
  dir.create("output/temp_rds_gpkg", recursive = TRUE)
}

saveRDS(params, "output/temp_rds_gpkg/params.rds")

# Functions --------------------------------------------------------------------

.check_gsim_cover <- function(file_in, yr_begin = -Inf, yr_end = Inf) {

  nrow_skip <- which(str_sub(read_lines(file_in), end = 1) == "#") %>%
    max()

  data_in <- read_csv(file_in, skip = nrow_skip, show_col_types = FALSE) %>%
    as_tibble() %>%
    select(date, MEAN, n.available) %>%
    mutate(month = month(date),
           year = year(date)) %>%
    filter(year >= yr_begin & year <= yr_end)

  if (nrow(data_in) == 0) { return (NULL) }

  res <- c(file_in,
           min(data_in$year),
           max(data_in$year),
           data_in %>% filter(is.na(MEAN)) %>% nrow()) %>%
    setNames(c("file", "yr_begin", "yr_end", "missing_obs"))

  return (res)

}

# GSIM data queries ------------------------------------------------------------

message("querying GSIM catchments...")

gsim_ts_files <- list.files("Data/GSIM/GSIM_indices/TIMESERIES/monthly",
                            recursive = TRUE, full.names = TRUE)

gsim_data_cover <- vector(mode = "list", length = length(gsim_ts_files))
for (i in 1:length(gsim_ts_files)) {

  if (i %% 1000 == 0) {message(paste0(i, " GSIM catchments checked at ", Sys.time()))}
  gsim_data_cover[[i]] <- .check_gsim_cover(gsim_ts_files[i],
                                            yr_begin = params$yr_begin,
                                            yr_end = params$yr_end)

}

gsim_selected_data <- gsim_data_cover %>%
  bind_rows() %>%
  mutate(across(-file, ~ as.numeric(.x)),
         record_max_length = yr_end - yr_begin + 1,
         missing_obs_fraction = missing_obs / (12 * record_max_length),
         monthly_obs_available = (1 - missing_obs_fraction) * record_max_length * 12,
         id = file %>% str_split("/") %>% lapply(tail, 1) %>% unlist() %>% str_replace(".mon", "")) %>%
  filter(missing_obs_fraction < params$max_missing_obs_fraction &
         record_max_length >= params$min_record_length_yr)

# exclude suspect stations
suspect_stations <- read_csv("Data/GSIM/GSIM_metadata/GSIM_catalog/GSIM_suspect_coordinates_stations.csv",
                             show_col_types = FALSE)

# exclude small catchments
small_catchments <- read_csv("Data/GSIM/GSIM_metadata/GSIM_catalog/GSIM_catchment_characteristics.csv",
                             show_col_types = FALSE) %>%
  filter(area.meta < params$min_catchment_area | area.est < params$min_catchment_area)

gsim_selected_data <- gsim_selected_data %>%
  filter(!id %in% suspect_stations$gsim.no & !id %in% small_catchments$gsim.no)

saveRDS(gsim_selected_data, "output/temp_rds_gpkg/01_catchment_records.rds")


