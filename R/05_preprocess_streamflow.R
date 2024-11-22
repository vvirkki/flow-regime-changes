# ------------------------------------------------------------------------------
# 
# 05_preprocess_streamflow
# 
# Fetch streamflow records and fill missing monthly values for the selected
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
library(doParallel)
options(dplyr.summarise.inform = FALSE)
params <- readRDS("output/temp_rds_gpkg/params.rds")

# Functions --------------------------------------------------------------------

.process_gsim_streamflow <- function(file_in, params, yr_begin = -Inf, yr_end = Inf) {

  id <- file_in %>%
    str_split("/") %>%
    unlist() %>%
    tail(1) %>%
    str_replace(".mon", "")

  lines_in <- read_lines(file_in)
  nrow_meta <- which(str_sub(lines_in, end = 1) == "#") %>%
    max()

  Q_monthly <- read_csv(file_in, skip = nrow_meta, show_col_types = FALSE) %>%
    select(date, MEAN, n.available) %>%
    mutate(id = id,
           year = year(date),
           month = month(date),
           n.available = ifelse(MEAN < 0, 0, n.available),
           MEAN = ifelse(MEAN < 0, NA, MEAN)) %>% # negative Q --> "missing"
    filter(year >= yr_begin & year <= yr_end) %>%
    mutate(streamflow_qflag = case_when(
      is.na(n.available) ~ "missing",
      n.available == 0 ~ "missing",
      n.available < 10 ~ "low",
      n.available < 20 ~ "medium",
      TRUE ~ "high"
    )) %>%
    rename(Q_m3s = MEAN) %>%
    select(id, year, month, Q_m3s, streamflow_qflag) %>%
    arrange(year, month)

  if (nrow(Q_monthly) != (max(Q_monthly$year) - min(Q_monthly$year) + 1) * 12) {
    message(paste0("mismatch in data counts in ", id))
  }
  return (Q_monthly)

}

.convert_m3s_to_mmPerMonth <- function(Q_m3s, A_km2, month) {

  sec_in_month <- days_in_month(month) * 24 * 60 * 60
  area_m2 <- A_km2 * 10^6 # km2 to m2
  Q_ms <- Q_m3s / area_m2 # m3/s to m/s
  Q_mmPerMonth <- Q_ms * sec_in_month * 1000 # m/s to mm/month
  return (Q_mmPerMonth)

}


.fill_missing_months <- function(ctm, params) {

  y1 <- min(ctm$year)
  yi <- max(ctm$year)
  month_id <- y1:yi %>%
    lapply(function(x) { paste0(x, "-", str_pad(seq(1,12,1), 2, pad = "0")) }) %>%
    unlist()

  ctm_id <- unique(ctm$id)
  ctm <- ctm %>%
    mutate(month_id = paste0(year, "-", str_pad(month, 2, pad = "0")))

  all_months <- tibble(month_id = month_id) %>%
    left_join(ctm, by = "month_id") %>%
    mutate(date = parse_date(month_id, format = "%Y-%m"),
           year = year(date),
           month = month(date)) %>%
    select(-date)

  missing_months <- all_months %>%
    filter(is.na(streamflow))

  if (nrow(missing_months) == 0) {
    return (ctm %>% select(-month_id))
  }

  filled_months <- tibble()
  for (i in 1:nrow(missing_months)) {
    row <- missing_months[i,]
    month_id <- row$month_id
    fill_Q <- ctm %>%
      filter(year <= row$year + params$incomplete_fill_nyrs &
             year >= row$year - params$incomplete_fill_nyrs &
             month == row$month & !is.na(streamflow)) %>%
      pull(streamflow) %>%
      mean()
    if (is.nan(fill_Q)) {
      fill_Q <- ctm %>%
        filter(month == row$month & !is.na(streamflow)) %>%
        pull(streamflow) %>%
        mean()
    }
    filled_months <- filled_months %>%
      bind_rows(tibble(month_id, fill_Q))
  }

  ctm_filled <- all_months %>%
    left_join(filled_months, by = "month_id") %>%
    mutate(streamflow = ifelse(is.na(streamflow), fill_Q, streamflow),
           id = ctm_id) %>%
    select(-c(month_id, fill_Q))

  return (ctm_filled)

}

# GSIM data preprocessing ------------------------------------------------------

gsim_selected_data <- readRDS("output/temp_rds_gpkg/01_catchment_records.rds")

message("fetching streamflow records...")

gsim_mmf_m3s <- vector(mode = "list", length = nrow(gsim_selected_data))
for (i in 1:nrow(gsim_selected_data)) {

  if (i %% 500 == 0) {message(paste0(i, " catchments processed at ", Sys.time()))}
  gsim_mmf_m3s[[i]] <- .process_gsim_streamflow(gsim_selected_data$file[i],
                                                params = params,
                                                yr_begin = params$yr_begin - params$incomplete_fill_nyrs,
                                                yr_end = params$yr_end + params$incomplete_fill_nyrs)

}

ctm_areas <- sf::read_sf("output/temp_rds_gpkg/03_catchment_geoms.gpkg") %>%
  select(id, area_km2) %>%
  sf::st_drop_geometry()

streamflow <- gsim_mmf_m3s %>%
  bind_rows() %>%
  left_join(ctm_areas, by = "id") %>%
  mutate(streamflow = .convert_m3s_to_mmPerMonth(Q_m3s, area_km2, month)) %>%
  select(id, year, month, streamflow, streamflow_qflag)

# Fill missing months ----------------------------------------------------------

ids <- streamflow$id %>%
  unique()

registerDoParallel(cores = detectCores() / 2)
message("doing filling in parallel...")
catchments_filled <- foreach(i = 1:length(ids),
                             .errorhandling = "pass",
                             .verbose = FALSE) %dopar%
{
  if (i %% 500 == 0) {message(paste0(i, " catchments processed at ", Sys.time()))}
  ctm <- streamflow %>%
    filter(id == ids[i])
  return (.fill_missing_months(ctm, params))
}
stopImplicitCluster()

streamflow_filled <- catchments_filled %>%
  bind_rows() %>%
  filter(year >= params$yr_begin & year <= params$yr_end)

# catchments that couldn't be filled (one or more months missing throughout record)
nan_streamflow_catchments <- streamflow_filled %>%
  filter(is.nan(streamflow)) %>%
  pull(id) %>%
  unique()

message(paste0(length(nan_streamflow_catchments), " catchments discarded due to empty months..."))

streamflow_filled <- streamflow_filled %>%
  filter(!id %in% nan_streamflow_catchments)

saveRDS(streamflow_filled, "output/temp_rds_gpkg/05_monthly_streamflow.rds")

readRDS("output/temp_rds_gpkg/01_catchment_records.rds") %>%
  filter(!id %in% nan_streamflow_catchments) %>%
  saveRDS("output/temp_rds_gpkg/01_catchment_records.rds")

readRDS("output/temp_rds_gpkg/02_catchment_attributes.rds") %>%
  filter(!id %in% nan_streamflow_catchments) %>%
  saveRDS("output/temp_rds_gpkg/02_catchment_attributes.rds")

sf::st_read("output/temp_rds_gpkg/03_catchment_geoms.gpkg") %>%
  filter(!id %in% nan_streamflow_catchments) %>%
  sf::st_write("output/temp_rds_gpkg/03_catchment_geoms.gpkg",
           delete_dsn = TRUE)

# Test that catchment ids are the same in all data files so far ----------------

ids1 <- readRDS("output/temp_rds_gpkg/01_catchment_records.rds") %>% pull(id)
ids2 <- readRDS("output/temp_rds_gpkg/02_catchment_attributes.rds") %>% pull(id)
ids3 <- sf::read_sf("output/temp_rds_gpkg/03_catchment_geoms.gpkg") %>% pull(id)
ids4 <- readRDS("output/temp_rds_gpkg/05_monthly_streamflow.rds") %>% pull(id) %>% unique()

test_matching_ids <- all.equal(ids1, ids2, ids3, ids4)

if (test_matching_ids) {
  message(paste0("sample size: ", length(ids1)))
  message("data preparation until streamflow ok!")
} else {
  stop("error in data preparation before streamflow!")
}


