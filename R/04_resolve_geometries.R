# ------------------------------------------------------------------------------
# 
# 04_resolve_geometries
# 
# Overlay the geometries of the selected catchment sample, determine nested
# catchments and remove duplicated catchments from the sample.
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

params <- readRDS("output/temp_rds_gpkg/params.rds")
ctm_attr <- readRDS("output/temp_rds_gpkg/02_catchment_attributes.rds")
geoms_in <- read_sf("output/temp_rds_gpkg/03_catchment_geoms.gpkg")

# Check up-to-date GRDC meta for including some caution catchments -------------

message("namematching with GRDC database...")

ctm_caution_river_station <- ctm_attr %>%
  filter(GSIM_geom_qflag == "caution") %>%
  select(id, country, river, station, GSIM_area_km2) %>%
  mutate(river = tolower(river),
         river = str_replace_all(river, "-", " "),
         river = str_replace_all(river, ",|\\(|\\)| river| creek| wash|gawa|rio ", ""),
         station = tolower(station),
         station = str_replace_all(station, "-", " "),
         station = str_replace_all(station, "\\(|\\)|at |near |below |above |hrs | downs|,", "")) %>%
  rowwise() %>%
  mutate(river =  str_split(river, " ") %>% unlist() %>% paste0(collapse = "|"),
         station = str_split(station, " ") %>% unlist() %>% paste0(collapse = "|")) %>%
  ungroup()

grdc_stations <- readxl::read_xlsx("Data/GRDC_Stations.xlsx") %>%
  select(grdc_no, country, river, station, area) %>%
  mutate(river = tolower(river),
         river = str_replace_all(river, "-", " "),
         river = str_replace_all(river, ",|\\(|\\)| river| creek| wash|gawa| gawa|rio ", ""),
         station = tolower(station),
         station = str_replace_all(station, "-", " "),
         station = str_replace_all(station, "\\(|\\)|at |near |below |above |hrs | downs|,", "")) %>%
  rowwise() %>%
  mutate(river =  str_split(river, " ") %>% unlist() %>% paste0(collapse = "|"),
         station = str_split(station, " ") %>% unlist() %>% paste0(collapse = "|")) %>%
  ungroup()

lst <- vector(mode = "list", length = nrow(ctm_caution_river_station))
for (i in 1:nrow(ctm_caution_river_station)) {

  ctm <- ctm_caution_river_station[i,]
  grdc_refs <- grdc_stations %>%
    rowwise() %>%
    filter(country == ctm$country & grepl(ctm$river, river)) %>% # same country, caution river is substring of GRDC river
    mutate(station_match = grepl(ctm$station, station), # caution station is substring of GRDC station
           river_inverse_match = ifelse(nrow(.) > 0, grepl(river, ctm$river), NA)) %>% # GRDC river is substring of caution river
    ungroup()

  out <- grdc_refs %>%
    mutate(gsim_id = ctm$id,
           gsim_river = ctm$river,
           gsim_station = ctm$station,
           gsim_area = ctm$GSIM_area_km2,
           area_diff = (gsim_area - area) / area) %>%
    arrange(abs(area_diff))

  lst[[i]] <- out

}

matched_ctm <- lst %>%
  bind_rows()

# same river and station, close enough catchment area --> match
full_match <- matched_ctm %>%
  filter(river_inverse_match & station_match & abs(area_diff) < 0.5)

GRDC_matches <- full_match %>%
  group_by(grdc_no) %>% # many GSIM for one GRDC
  arrange(abs(area_diff), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(gsim_id) %>% # many GRDC for one GSIM
  arrange(abs(area_diff), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(grdc_no, gsim_id, country, area, gsim_area, area_diff,
         river, gsim_river, station, gsim_station)

message(paste0(nrow(GRDC_matches), " catchments matched with GRDC data..."))

# include non-caution quality catchments and GRDC matched catchments
GRDC_matched_stations <- GRDC_matches %>%
  mutate(GRDC_matched = TRUE) %>%
  select(gsim_id, GRDC_matched)

sel <- geoms_in %>%
  left_join(ctm_attr %>% select(id, GSIM_geom_qflag), by = "id") %>%
  left_join(GRDC_matched_stations, by = c("id" = "gsim_id")) %>%
  filter(GSIM_geom_qflag != "caution" | GRDC_matched) %>%
  select(-GRDC_matched)

# Preapre for intersections ----------------------------------------------------

# build centroid index
sel_centroids <- suppressWarnings(
  sel %>%
    select(id) %>%
    st_centroid())

# check for catchments that don't intersect any other catchment
is_mat <- st_intersects(sel, sel, sparse = FALSE)
count_intersections <- is_mat %>%
  rowSums()

fully_detached_geoms <- sel[which(count_intersections == 1),] %>%
  pull(id)
remove(is_mat)

# catchments that need resolving
sel_resolve <- sel %>%
  filter(!id %in% fully_detached_geoms) %>%
  select(id, area_km2)

overlaps_remaining <- st_intersects(sel_resolve, sel_resolve)
if (min(rowSums(as.matrix(overlaps_remaining))) < 2) { stop("intersection doesn't work") }

# Perform intersections --------------------------------------------------------

intersections <- vector(length = length(overlaps_remaining), mode = "list")
message(paste0("resolving ", length(overlaps_remaining), " catchment geometries..."))
message("performing intersections...")

for (i in 1:length(overlaps_remaining)) {

  if (i %% 500 == 0) {message(paste0(i, " catchments processed..."))}
  is_to <- overlaps_remaining[[i]]
  intersecting_catchments <- sel_resolve[is_to[-which(is_to == i)],] # exclude self-intersection

  source_area <- sel_resolve[i,]$area_km2
  ist <- suppressWarnings(
    st_intersection(sel_resolve[i,], intersecting_catchments) %>%
      rename(ist_id = id.1) %>%
      select(id, ist_id) %>%
      mutate(ist_area_km2 = as.numeric(st_area(.)) / 10^6,
             ist_area_frac = ist_area_km2 / source_area) %>%
      select(id, ist_id, ist_area_frac) %>%
      st_drop_geometry()
  )

  centroid_within_ist <- sel_centroids %>%
    filter(id %in% ist$ist_id) %>%
    st_within(sel_resolve[i,], sparse = FALSE) %>%
    as.vector()

  ist <- ist %>%
    mutate(centroid_within = centroid_within_ist)

  if (nrow(ist) == 0) {
    # only touches with perfect line, no intersection geometries
    fully_detached_geoms <- c(fully_detached_geoms, sel_resolve[i,]$id)
  } else {
    intersections[[i]] <- ist
  }

}

# Determine parents and children -----------------------------------------------

intersections_df <- intersections %>% bind_rows()
ids <- unique(intersections_df$id)
ctm_resolved_list <- vector(mode = "list", length = length(ids))
message("determining parents and children...")

for (i in 1:length(ids)) {

  if (i %% 500 == 0) {message(paste0(i, " catchments processed..."))}
  ctm <- intersections_df %>%
    filter(id == ids[i] &
           (ist_area_frac > params$border_intersection_threshold | centroid_within)) # exclude sliver intersections

  # no intersections except for possible slivers --> detached catchment
  detached <- ifelse(nrow(ctm) == 0, TRUE, FALSE)
  child_ctms <- c()
  parent_ctms <- c()

  if (!detached) {

    parent_ctms <- ctm %>%
      filter(ist_area_frac > params$nested_intersection_threshold) %>%
      pull(ist_id)

    # check if a catchment is duplicate (child of one of its parent catchments)
    duplicate_to_ctms <- c()
    if (length(parent_ctms) > 0) {
      for (j in 1:length(parent_ctms)) {
        parent <- intersections_df %>%
          filter(id == parent_ctms[j] & ist_id == ids[i])
        if (parent$ist_area_frac > params$nested_intersection_threshold) {
          duplicate_to_ctms <- c(duplicate_to_ctms, parent_ctms[j])
        }
      }
    }

    duplicate <- ifelse(length(duplicate_to_ctms) > 0, TRUE, FALSE)
    parent_ctms <- parent_ctms[!parent_ctms %in% duplicate_to_ctms]

    child_ctms <- ctm %>%
      filter(ist_area_frac <= params$nested_intersection_threshold) %>%
      pull(ist_id)

  }

  ctm_out <- tibble(id = ids[i],
                    duplicate = duplicate,
                    parents = ifelse(length(parent_ctms) != 0, paste(parent_ctms, collapse = "|"), NA),
                    children = ifelse(length(child_ctms) != 0, paste(child_ctms, collapse = "|"), NA),
                    duplicate_to = ifelse(duplicate, paste(duplicate_to_ctms, collapse = "|"), NA))
  ctm_resolved_list[[i]] <- ctm_out

}

ctm_fully_detached <- fully_detached_geoms %>%
  tibble(id = .) %>%
  mutate(duplicate = FALSE, parents = NA, children = NA, duplicate_to = NA)

ctm_resolved <- ctm_resolved_list %>%
  bind_rows() %>%
  bind_rows(ctm_fully_detached) %>%
  arrange(id)

# Resolve duplicates by picking one with the most observations -----------------

ctm_with_duplicates <- ctm_resolved %>%
  filter(duplicate) %>%
  pull(id)

monthly_obs_counts <- readRDS("output/temp_rds_gpkg/01_catchment_records.rds") %>%
  select(id, monthly_obs_available)

duplicate_selections <- vector(mode = "list", length = length(ctm_with_duplicates))
message("removing duplicates...")
for (i in 1:length(ctm_with_duplicates)) {

  ctm <- ctm_resolved %>%
    filter(id == ctm_with_duplicates[i])

  duplicates <- ctm_resolved %>%
    filter(id == ctm_with_duplicates[i] | grepl(ctm$duplicate_to, id)) %>%
    pull(id)

  duplicate_ctms <- sel %>%
    filter(id %in% duplicates) %>%
    st_drop_geometry() %>%
    left_join(monthly_obs_counts, by = "id")

  # select catchment with the most observations out of all duplicate records
  selected_duplicate <- duplicate_ctms %>%
    arrange(-monthly_obs_available) %>%
    slice(1) %>%
    pull(id)

  duplicate_selections[[i]] <- tibble(id = ctm$id,
                                      selected_duplicate = selected_duplicate)

}

discarded_duplicates <- duplicate_selections %>%
  bind_rows() %>%
  filter(!id %in% .$selected_duplicate) %>%
  pull(id)

duplicates_preserved <- duplicate_selections %>%
  bind_rows() %>%
  pull(selected_duplicate) %>%
  unique() %>%
  length()

message(paste0(duplicates_preserved, " duplicates preserved..."))
message(paste0(length(discarded_duplicates), " redundant duplicates discarded..."))

# Update parents and children after duplicate removal --------------------------

message("updating parents and children after duplicate removal...")
for (i in 1:length(discarded_duplicates)) {

  ctm <- ctm_resolved %>%
    filter(id == discarded_duplicates[i])

  parents_ids <- ctm$parents %>%
    str_split("\\|") %>%
    unlist()

  # remove duplicate children from all parents
  if (!any(is.na(parents_ids))) {

    for (j in 1:length(parents_ids)) {

      correct_rowno <- which(parents_ids[j] == ctm_resolved$id)
      parent_ctm_resolved <- ctm_resolved[correct_rowno,]
      old_children <- parent_ctm_resolved$children %>%
        str_split("\\|") %>%
        unlist()
      parent_ctm_resolved$children <- paste(old_children[ctm$id != old_children], collapse = "|")
      ctm_resolved[correct_rowno,] <- parent_ctm_resolved

    }

  }

  children_ids <- ctm$children %>%
    str_split("\\|") %>%
    unlist()

  # remove duplicate parents from all children
  if (!any(is.na(children_ids))) {

    for (j in 1:length(children_ids)) {

      correct_rowno <- which(children_ids[j] == ctm_resolved$id)
      child_ctm_resolved <- ctm_resolved[correct_rowno,]
      old_parents <- child_ctm_resolved$parents %>%
        str_split("\\|") %>%
        unlist()
      child_ctm_resolved$parents <- paste(old_parents[ctm$id != old_parents], collapse = "|")
      ctm_resolved[correct_rowno,] <- child_ctm_resolved

    }
  }
}

ctm_resolved_no_duplicates <- ctm_resolved %>%
  mutate(parents = ifelse(parents == "", NA, parents), # all parents got removed
         children = ifelse(children == "", NA, children)) %>% # all children got removed
  filter(!id %in% discarded_duplicates) %>%
  select(-c(duplicate, duplicate_to)) %>%
  mutate(parents = ifelse(parents == "NA", NA, parents))

# check that children are correctly marked to all their parents
message("fixing parents that miss children...")
for (i in 1:nrow(ctm_resolved_no_duplicates)) {

  if (i %% 500 == 0) {message(paste0(i, " catchments checked..."))}
  ctm <- ctm_resolved_no_duplicates[i,]
  parents_missing_child <- ctm_resolved_no_duplicates %>%
    filter(grepl(ctm$parents, id) & !grepl(ctm$id, children))

  if (nrow(parents_missing_child) > 0) {
    message(paste0("adding missing children to ", ctm$id, "..."))
    for (j in 1:nrow(parents_missing_child)) {

      idx <- which(ctm_resolved_no_duplicates$id == parents_missing_child[j,]$id)
      updated_parent <- ctm_resolved_no_duplicates[idx,] %>%
        mutate(children = paste0(children, "|", ctm$id))
      ctm_resolved_no_duplicates[idx,] <- updated_parent

    }
  }
}

# Update records, attributes and geometries to resolved non-duplicates ---------

ctm_selection_updated <- readRDS("output/temp_rds_gpkg/01_catchment_records.rds") %>%
  filter(id %in% ctm_resolved_no_duplicates$id) %>%
  saveRDS("output/temp_rds_gpkg/01_catchment_records.rds")

attributes_updated <- readRDS("output/temp_rds_gpkg/02_catchment_attributes.rds") %>%
  filter(id %in% ctm_resolved_no_duplicates$id) %>%
  saveRDS("output/temp_rds_gpkg/02_catchment_attributes.rds")

sel_updated <- read_sf("output/temp_rds_gpkg/03_catchment_geoms.gpkg") %>%
  filter(id %in% ctm_resolved_no_duplicates$id) %>%
  left_join(ctm_resolved_no_duplicates, by = "id") %>%
  relocate(geom, .after = last_col()) %>%
  st_write("output/temp_rds_gpkg/03_catchment_geoms.gpkg", delete_dsn = TRUE)



