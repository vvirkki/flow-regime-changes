# ------------------------------------------------------------------------------
# 
# 08_determine_trends
# 
# Estimate annual metrics and Theil-Sen slopes for streamflow and drivers in all
# catchments of the selected catchment sample.
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
library(zoo)
library(doParallel)
library(RobustLinearReg)
options(dplyr.summarise.inform = FALSE)

# Parameters -------------------------------------------------------------------

params <- readRDS("output/temp_rds_gpkg/params.rds")

# Functions --------------------------------------------------------------------

.get_trend_theil_sen <- function(x, p) {

  tsdata <- data.frame(1:length(x), x) %>%
    setNames(c("x", "y"))
  tsmodel <- theil_sen_regression(y ~ x, tsdata)

  mdl <- summary(tsmodel)
  xest <- mdl$coefficients["x",][["Estimate"]]

  pval <- mdl$coefficients["x",][["Pr(>|t|)"]]
  return (tibble(xest, pval, pval < p) %>% setNames(c("slope", "pval", "significant")))

}

# Setup ------------------------------------------------------------------------

signal_type <- c("mean", "sd", "quantile05", "quantile95")

message(paste0("creating streamflow trends at ", Sys.time(), "..."))

streamflow <- readRDS("output/temp_rds_gpkg/05_monthly_streamflow.rds") %>%
  mutate(timestep = make_date(year, month),
         variable = "streamflow")

# Annual streamflow trends -----------------------------------------------------

annual_values <- streamflow %>%
  group_by(id, year) %>%
  summarise(annual_mean = mean(streamflow),
            annual_sd = sd(streamflow),
            annual_quantile05 = quantile(streamflow, 0.05),
            annual_quantile95 = quantile(streamflow, 0.95)) %>%
  ungroup()

# test for linear trends using Theil-Sen
ids <- unique(streamflow$id)
annual_trends_lst <- vector(mode = "list", length = length(ids))

for (i in 1:length(ids)) {

  ctm_trends <- annual_values %>%
    filter(id == ids[i]) %>%
    select(-c(id, year)) %>%
    lapply(.get_trend_theil_sen, params$pval)
  
  # par(mfrow = c(4,1))
  # plot(ctm_trends$annual_mean$tsdata$x, ctm_trends$annual_mean$tsdata$y, type = "l",
  #      xlab = "", ylab = "")
  # abline(ctm_trends$annual_mean$tsmodel, col = "red")
  # plot(ctm_trends$annual_sd$tsdata$x, ctm_trends$annual_sd$tsdata$y, type = "l",
  #      xlab = "", ylab = "")
  # abline(ctm_trends$annual_sd$tsmodel, col = "red")
  # plot(ctm_trends$annual_quantile95$tsdata$x, ctm_trends$annual_quantile95$tsdata$y, type = "l",
  #      xlab = "", ylab = "")
  # abline(ctm_trends$annual_quantile95$tsmodel, col = "red")
  # plot(ctm_trends$annual_quantile05$tsdata$x, ctm_trends$annual_quantile05$tsdata$y, type = "l",
  #      xlab = "", ylab = "")
  # abline(ctm_trends$annual_quantile05$tsmodel, col = "red")

  annual_trends_lst[[i]] <- ctm_trends %>%
    bind_rows() %>%
    mutate(id = ids[i],
           variable = "streamflow",
           signal_type = str_replace(names(ctm_trends), "annual_", ""),
           testvar = "annual",
           method = "theil_sen",
           unit = "mmyr") %>%
    select(id, variable, signal_type, testvar, method,
           slope, unit, pval, significant)

}


saveRDS(list(annual_values, bind_rows(annual_trends_lst)) %>%
               setNames(c("annual_values", "trends")),
        "output/temp_rds_gpkg/08_streamflow_annual_values_trends_theil_sen.rds")

# Parallel driver trend building -----------------------------------------------

drivers_continuous <- c("e_total", "total_prec", "wateruse_all_sectors")
driver_signal_setups <- expand_grid(driver = drivers_continuous, signal_type)

message(paste0("creating driver trends in parallel at ", Sys.time(), "..."))

registerDoParallel(cores = 8)
driver_signals <- foreach(i = 1:nrow(driver_signal_setups),
                          .errorhandling = "pass",
                          .verbose = FALSE) %dopar%
  {

    setup <- driver_signal_setups[i,]
    driver <- readRDS(paste0("output/temp_rds_gpkg/07_", setup$driver, ".rds")) %>%
      mutate(timestep = make_date(year, month))

    driver_annual_out <- paste0("output/temp_rds_gpkg/08_",
                                setup$driver, "_annual_values_trends_theil_sen.rds")
    if (!file.exists(driver_annual_out)) {

      # annual values
      annual_values <- driver %>%
        group_by(id, year) %>%
        summarise(annual_mean = mean(value),
                  annual_sd = sd(value),
                  annual_quantile05 = quantile(value, 0.05),
                  annual_quantile95 = quantile(value, 0.95)) %>%
        ungroup()

      # test for linear trends using Theil-Sen
      ids <- unique(driver$id)
      annual_trends_lst <- vector(mode = "list", length = length(ids))
      for (i in 1:length(ids)) {

        ctm_trends <- annual_values %>%
          filter(id == ids[i]) %>%
          select(-c(id, year)) %>%
          lapply(.get_trend_theil_sen, params$pval)

        annual_trends_lst[[i]] <- ctm_trends %>%
          bind_rows() %>%
          mutate(id = ids[i],
                 variable = setup$driver,
                 signal_type = str_replace(names(ctm_trends), "annual_", ""),
                 testvar = "annual",
                 method = "theil_sen",
                 unit = "mmyr") %>%
          select(id, variable, signal_type, testvar, method,
                 slope, unit, pval, significant)

      }
      saveRDS(list(annual_values, bind_rows(annual_trends_lst)) %>%
                setNames(c("annual_values", "trends")),
              driver_annual_out)

    }

  }

stopImplicitCluster()
