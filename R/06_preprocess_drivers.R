# ------------------------------------------------------------------------------
# 
# 06_preprocess_drivers
# 
# Pre-process driver data from raw data files to tif rasters (climate & water
# use) or a point data frame (dams).
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
library(ncdf4)
options(dplyr.summarise.inform = FALSE)

if (!dir.exists("output/temp_drivers_tif_gpkg")) {
  dir.create("output/temp_drivers_tif_gpkg", recursive = TRUE)
}

# Preprocess ERA5-Land ---------------------------------------------------------

.process_era5_land <- function(file_in) {

  r <- rast(file_in)
  r <- r[[names(r)[grepl("expver=1", names(r))]]]

  source_var <- file_in %>%
    str_split("/") %>%
    unlist()
  source_var <- source_var[grepl(".nc", source_var)] %>%
    str_replace(".nc", "") %>%
    str_replace_all("_", "-")

  nmr <- time(r) %>%
    format("%Y.%m") %>%
    paste0(source_var, "_", .)
  names(r) <- nmr
  crs(r) <- "+proj=lonlat"

  tmpl <- rast(nrows = nrow(r), ncols = ncol(r),
               xmin = xmin(r) - 180, xmax = xmax(r) - 180,
               ymin = ymin(r), ymax = ymax(r),
               crs = "+proj=lonlat")

  out_dir <- paste0("output/temp_drivers_tif_gpkg/", source_var)
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }

  nyrs <- floor(nlyr(r) / 12)
  land_sea_mask <- rast("data/land_sea_mask_for_era5.tif")

  for (i in 1:nyrs) {

    sel_layers <- ((i-1)*12+1):(i*12)
    rslice <- r[[sel_layers]]

    # hydrological variables from averaged m/day to summed mm/month
    # see: https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation#heading-Monthlymeans
    if (!grepl("temp_2m", file_in)) {
      rslice <- rslice * days_in_month(seq(1, 12, 1)) * 1000
      if (grepl("e_total", file_in)) {
        rslice <- -1 * rslice # evaporative fluxes to up-positive
      }
    }

    rslice_proj <- rslice %>%
      project(tmpl)
    rslice_proj[is.na(land_sea_mask)] <- NA

    yi <- rslice %>% time() %>% year() %>% max()
    out_path <- paste0(out_dir, "/", source_var, "_", yi, ".tif")
    writeRaster(rslice_proj, out_path)
    unlink(paste0(out_path, ".aux.json"))

    message(paste0("wrote ", out_path, " at ", Sys.time()))

  }
  return (NULL)
}

files <- c("Data/ERA5-Land/total_prec.nc",
           "Data/ERA5-Land/e_total.nc")

message("preparing ERA5-Land data...")

dummy <- files %>%
  lapply(.process_era5_land)

# Preprocess Huang et al. water use data ---------------------------------------

message("preparing water use data...")

sctr <- c("mining", "manufacturing", "livestock", "irrigation", "electricity",
          "domestic")
sctr_abb <- c("min", "mfg", "liv", "irr", "elec", "dom")
ghm <- c("watergap", "pcr", "jpjml", "h08")

for (i in 1:length(sctr)) {

  if (sctr[i] == "irrigation") {
    wu_files <- c()
    for (g in 1:length(ghm)) {
      nf <- paste0("Data/wateruse/", sctr[i] ," water use v2/cons_",
                   sctr_abb[i], "_", ghm[g], ".nc")
      wu_files <- c(wu_files, nf)
    }
  } else {
    wu_files <- paste0("Data/wateruse/", sctr[i] ," water use v2/cons_",
                       sctr_abb[i], ".nc")
  }

  irrigation_wu <- list()
  for (fileno in 1:length(wu_files)) {

    conn <- nc_open(wu_files[fileno])
    vals <- ncvar_get(conn, paste0("cons_", sctr_abb[i]))
    lon <- ncvar_get(conn, "lon")
    lat <- ncvar_get(conn, "lat")
    nmonths <- ncvar_get(conn, "month") %>% length()
    nc_close(conn)

    l <- vector(mode = "list", length = nmonths)
    for (colno in 1:nmonths) {

      grd_month <- cbind(lon, lat, vals[,colno])
      r_month <- rast(grd_month, type = "xyz", crs = "+proj=lonlat") %>%
        extend(ext(-180, 180, -90, 90))
      l[[colno]] <- r_month

    }
    if (sctr[i] == "irrigation") {
      irrigation_wu[[fileno]] <- l
    }

  }

  nm_mnths <- rep(seq(1971, 2010), 12) %>%
    sort() %>%
    paste0(".", str_pad(seq(1:12), 2, pad = "0")) %>%
    paste0(sctr_abb[i], "Wu_", .)

  if (length(irrigation_wu) > 0) {
    iwu <- irrigation_wu %>%
      lapply(rast)
    l <- (iwu[[1]] + iwu[[2]] + iwu[[3]] + iwu[[4]]) / 4
    wcons_sctr <- l %>%
      setNames(nm_mnths)
  } else {
    wcons_sctr <- rast(l) %>%
      setNames(nm_mnths)
  }

  out_dir <- "output/temp_drivers_tif_gpkg/wateruse"
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  out_path <- paste0(out_dir, "/wateruse_", sctr[i], ".tif")
  writeRaster(wcons_sctr, out_path)
  message(paste0("wrote ", out_path, " at ", Sys.time()))

}

# Preprocess dam data ----------------------------------------------------------

message("preparing dam data...")

dams <- read_sf("Data/GeoDAR_ICOLD2023_dams.geojson") %>% # proprietary data file
  select(id_v11, `Name of the dam`, Country, River, `Year of Completion`,
         `Reservoir Capacity`, `Catchment area`, `Dam Type`, Purposes) %>%
  setNames(c("id", "name", "country", "river", "year", "reservoir_capacity",
             "catchment_area", "dam_type", "purposes", "geometry")) %>%
  filter(!is.na(reservoir_capacity) & !is.na(year))

st_write(dams, "output/temp_drivers_tif_gpkg/dams.gpkg", delete_dsn = TRUE)

