# ------------------------------------------------------------------------------
# 
# res01_assign_analyse_frcs
# 
# Based on Theil-Sen slopes, assign catchments to flow regime change classes
# (FRCs) and aggregate trends and changes in drivers within FRC groups.
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
options(dplyr.summarise.inform = FALSE)

# FUNCTIONS --------------------------------------------------------------------

.get_trends <- function(trends_files_in, var) {

  trends_in <- trends_files_in %>%
    lapply(function(x) {readRDS(x)$trends}) %>%
    bind_rows() %>%
    filter(testvar == var) %>%
    mutate(trend_dir = case_when(
      slope > 0 ~ "a_incrTrend",
      slope < 0 ~ "b_decrTrend",
      slope == 0 ~ "zeroSlope"
    ))

  trends_combined <- trends_in %>%
    select(id, signal_type, trend_dir) %>%
    pivot_wider(id_cols = "id",
                names_from = "signal_type",
                values_from = "trend_dir")

  trends_all <- trends_in %>%
    select(-trend_dir)

  return (list(trends_combined, trends_all) %>%
            setNames(c("combined", "all")))

}

.plot_signal <- function(signals, ctm_id, testvar) {

  ctm <- signals %>%
    filter(id == ctm_id) %>%
    pivot_longer(-c(id, year), names_to = "signal_type") %>%
    rename(timestep = year) %>%
    mutate(signal_type = case_when(
      signal_type == "annual_mean" ~ "a_mean",
      signal_type == "annual_quantile05" ~ "d_quantile05",
      signal_type == "annual_quantile95" ~ "c_quantile95",
      signal_type == "annual_sd" ~ "b_sd",
    ))

  plt <- ggplot(ctm, aes(x = timestep, y = value)) +
    geom_line() +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ signal_type, scales = "free", nrow = 4) +
    ggtitle(ctm_id) +
    theme(axis.title.y = element_blank())

  return (plt)

}

.get_driver_signals_trends <- function(dr, var, ids) {

  fdr <- list.files("output/temp_rds_gpkg/", full.names = TRUE)
  drq <- paste0("_", dr, "_")
  dr_files_in <- fdr[grepl(dr, fdr) & grepl("theil_sen", fdr)]

  dr_signals_in <- dr_files_in %>%
    lapply(function(x) {readRDS(x)$annual_values}) %>%
    bind_rows() %>%
    filter(id %in% ids)

  dr_trends <- .get_trends(dr_files_in, var) %>%
    lapply(function(x) { x %>% filter(id %in% ids) })

  return (list(dr_signals_in, dr_trends) %>%
            setNames(c("signals", "trends")))

}

.get_driver <- function(requested_driver, frc_assignments) {

  driver <- .get_driver_signals_trends(requested_driver, "annual", frc_assignments$id)
  # # optional plotting
  # .plot_signal(driver$signals, "FI_0000096", "annual")
  return (list(trends = driver$trends))

}

.clean_sig <- function(x) {

  return (case_when(is.na(x) ~ "",
                    grepl("nonsignificant", x) ~ "",
                    TRUE ~ str_replace(x, "_a_significant", "")))

}

.bin_to_classes <- function(tbl, grp, bins) {

  coll <- tibble()
  for (i in 1:length(grp)) {
    a <- tbl %>%
      st_drop_geometry() %>%
      mutate(bin = cut((tbl[,grp[i]]) %>% pull(1), bins)) %>%
      group_by(bin) %>%
      summarise(nn = n()) %>%
      mutate(freq = nn / sum(nn),
             var = grp[i])
    coll <- coll %>% bind_rows(a)
  }
  return (coll)

}

.find_diverse_subctm <- function(tbl, row) {

  subctm_cases <- tbl %>%
    filter(grepl(tbl[row,]$children, id)) %>%
    group_by(case) %>%
    summarise(nn = n()) %>%
    mutate(ntot = sum(nn),
           freq = nn / ntot,
           id = tbl[row,]$id) %>%
    select(id, everything())
  return (subctm_cases)

}

.prep_prec_trend_plt <- function(prec_trend_data, outlier_limits, metric) {

  group_prep <- prec_trend_data %>%
    left_join(outlier_limits, by = c("case", "area_class")) %>%
    select(id, case, area_class, matches(metric)) %>%
    setNames(str_replace_all(colnames(.), paste0("_", metric), "")) %>%
    filter(value > minval & value < maxval)

  iqr_by_trend <- group_prep %>%
    group_by(case, prec_trend) %>%
    summarise(nn = n(),
              q25_value = quantile(value, 0.25),
              q75_value = quantile(value, 0.75)) %>%
    mutate(signal_type = metric) %>%
    ungroup()

  iqr_all <- group_prep %>%
    group_by(case) %>%
    summarise(nn = n(),
              q25_value = quantile(value, 0.25),
              q75_value = quantile(value, 0.75)) %>%
    mutate(signal_type = metric,
           prec_trend = "all") %>%
    select(colnames(iqr_by_trend)) %>%
    ungroup()

  group_mean_by_trend <- group_prep %>%
    group_by(case, prec_trend) %>%
    summarise(nn = n(),
              mean_value = mean(value)) %>%
    mutate(signal_type = metric) %>%
    ungroup()

  group_mean_all <- group_prep %>%
    group_by(case) %>%
    summarise(nn = n(),
              mean_value = mean(value)) %>%
    mutate(signal_type = metric,
           prec_trend = "all") %>%
    select(colnames(group_mean_by_trend)) %>%
    ungroup()

  ret <- list(iqr = bind_rows(iqr_by_trend, iqr_all),
              group_mean = bind_rows(group_mean_by_trend, group_mean_all))
  return (ret)

}

.prep_scatter <- function(trend_data, sel_var, metric) {

  trd <- trend_data %>%
    filter(variable == sel_var & signal_type == metric) %>%
    mutate(min_val = mean(value) - 3  * sd(value),
           max_val = mean(value) + 3 * sd(value)) %>%
    filter(value > min_val & value < max_val) %>%
    select(-c(min_val, max_val))

  return (trd)

}

# STREAMFLOW DATA IN -----------------------------------------------------------

params <- readRDS("output/temp_rds_gpkg/params.rds")
if ( !dir.exists("output/figures_tables") ) { dir.create("output/figures_tables") }

ctm_records <- readRDS("output/temp_rds_gpkg/01_catchment_records.rds")

sfin <- list.files("output/temp_rds_gpkg", full.names = TRUE)
sfin <- sfin[grepl("streamflow", sfin) & grepl("theil_sen", sfin)]

streamflow_signals_in <- sfin %>%
  lapply(function(x) { readRDS(x)$annual_values }) %>%
  bind_rows()

# # optional plotting
# .plot_signal(streamflow_signals_in, "FI_0000096", "annual")

# STREAMFLOW TRENDS ------------------------------------------------------------

gg <- read_sf("output/temp_rds_gpkg/03_catchment_geoms.gpkg") %>%
  arrange(-area_km2) %>%
  mutate(render_order = row_number()) %>%
  arrange(id)

figS1_data <- gg %>%
  left_join(ctm_records, by = "id") %>%
  select(id, area_km2, record_max_length, missing_obs_fraction, render_order)

st_write(figS1_data, "output/figures_tables/figS1_sample_characteristics.gpkg",
         delete_dsn = TRUE)

streamflow_trends <- .get_trends(sfin, "annual")

streamflow_component_trends <- streamflow_trends$all %>%
  pivot_wider(id_cols = "id", names_from = "signal_type", values_from = "slope") %>%
  left_join(x = gg, y = ., by = "id")

st_write(streamflow_component_trends, "output/figures_tables/fig2_streamflow_trend_components.gpkg",
         delete_dsn = TRUE)

# for numbers/statements in Section 3.1.
streamflow_components_in_bins <- .bin_to_classes(streamflow_component_trends,
                                                 c("mean", "sd", "quantile05", "quantile95"),
                                                 c(-Inf, 0, Inf))

streamflow_component_sig <- streamflow_trends$all %>%
  pivot_wider(id_cols = "id", names_from = "signal_type", values_from = "significant") %>%
  mutate(across(c(mean, sd, quantile05, quantile95), ~ ifelse(is.na(.x), 0, as.numeric(.x))),
         nmbr_sig = mean + sd + quantile05 + quantile95) %>%
  left_join(x = gg, y = ., by = "id")

st_write(streamflow_component_sig, "output/figures_tables/figS3_streamflow_component_significance.gpkg",
         delete_dsn = TRUE)

# FRC ASSIGNMENT ---------------------------------------------------------------

frc_down <- streamflow_trends$combined %>%
  filter(mean == "b_decrTrend" & quantile05 == "b_decrTrend" & quantile95 == "b_decrTrend") %>%
  mutate(case = "aa_shift_down")

frc_up <- streamflow_trends$combined %>%
  filter(mean == "a_incrTrend" & quantile05 == "a_incrTrend" & quantile95 == "a_incrTrend") %>%
  mutate(case = "ab_shift_up")

frc_shrink <- streamflow_trends$combined %>%
  filter(sd == "b_decrTrend" & quantile05 == "a_incrTrend" & quantile95 == "b_decrTrend") %>%
  mutate(case = "ac_shrink")

frc_expand <- streamflow_trends$combined %>%
  filter(sd == "a_incrTrend" & quantile05 == "b_decrTrend" & quantile95 == "a_incrTrend") %>%
  mutate(case = "ad_expand")

frc_assignments <- bind_rows(frc_shrink, frc_expand, frc_down, frc_up) %>%
  select(id, case) %>%
  arrange(id)

cases_geoms <- gg %>%
  left_join(frc_assignments, by = "id") %>%
  mutate(case = ifelse(is.na(case), "ae_other", case))

st_write(cases_geoms, "output/figures_tables/fig3_frc_assignments.gpkg", delete_dsn = TRUE)

# for Table S1
frc_trend_counts <- cases_geoms %>%
  st_drop_geometry() %>%
  select(id, case) %>%
  left_join(streamflow_trends$combined, by = "id") %>%
  pivot_longer(-c(id, case), names_to = "metric") %>%
  group_by(case, metric, value) %>%
  summarise(nn = n())

write_csv(frc_trend_counts, "output/figures_tables/tableS1.csv")

# finding diverse FRC assignments in catchments for Section 3.1.
check_subctm <- cases_geoms %>%
  st_drop_geometry() %>%
  filter(case != "ae_other")

diverse_subctm <- vector(length = nrow(check_subctm), mode = "list")
for (i in 1:nrow(check_subctm)) {
  diverse_subctm[[i]] <- .find_diverse_subctm(check_subctm, i)
}

nms <- readRDS("output/temp_rds_gpkg/02_catchment_attributes.rds")
large_basins_many_subctm <- diverse_subctm %>%
  bind_rows() %>%
  left_join(nms, by = "id") %>%
  filter(ntot > 40)

ctm_examples_1 <- large_basins_many_subctm %>%
  group_by(id) %>%
  arrange(-freq, .by_group = TRUE) %>%
  slice(1) %>%
  arrange(freq) %>%
  left_join(nms, by = "id")

ctm_examples_2 <- large_basins_many_subctm %>%
  group_by(id) %>%
  summarise(range = max(freq) - min(freq)) %>%
  left_join(nms, by = "id")

# FRC COUNTS IN CATCHMENT AREA CLASSES -----------------------------------------

cases_nogeom <- cases_geoms %>%
  st_drop_geometry() %>%
  mutate(area_class = case_when(
    area_km2 < 2500 ~ "a_small",
    area_km2 < 10000 ~ "b_medium",
    area_km2 > 10000 ~ "c_large"
  ))

counts <- cases_nogeom %>%
  group_by(case, area_class) %>%
  summarise(nn = n(),
            area_class = unique(area_class)) %>%
  ungroup()

ctm_size_breakdown <- ggplot(counts, aes(x = case, y = nn, fill = case)) +
  geom_col() +
  facet_grid(cols = vars(area_class), space = "free_x") +
  scale_y_continuous(limits = c(0,450)) +
  scale_fill_manual(values = c("#f4c40f", "#1f6e9d", "#aa7aa1", "#92c051", "#dfd7b2")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ylab("number of catchments") +
  theme(axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank())

ggsave("output/figures_tables/fig3_ctm_size_breakdown.pdf", ctm_size_breakdown,
       width = 135, height = 60, units = "mm")

# DRIVER DATA IN ---------------------------------------------------------------

all_drivers <- c("e_total", "total_prec", "wateruse_all_sectors")

ctm_frc <- cases_nogeom %>%
  select(id, case, area_class)

driver_lst <- lapply(all_drivers, .get_driver,
                     frc_assignments = ctm_frc) %>%
  setNames(all_drivers)

driver_trends <- driver_lst %>%
  lapply('[[', 1) %>%
  lapply(function(x) {x$all}) %>%
  bind_rows() %>%
  filter(variable %in% c("e_total", "total_prec", "wateruse_all_sectors")) %>%
  select(id, variable, signal_type, slope, significant) %>%
  rename(value = slope)

# for section 3.2
clim_variable_trends <- driver_trends %>%
  filter(variable %in% c("e_total", "total_prec") & signal_type %in% c("mean", "sd")) %>%
  mutate(direction = case_when(
    value < 0 ~ "decreasing",
    value > 0 ~ "increasing",
    TRUE ~ NA
  )) %>%
  group_by(variable, signal_type, direction) %>%
  summarise(nn = n())

fig_s4_data <- tibble() %>%
  bind_rows(.prep_scatter(driver_trends, "e_total", "mean"),
            .prep_scatter(driver_trends, "e_total", "sd"),
            .prep_scatter(driver_trends, "total_prec", "mean"),
            .prep_scatter(driver_trends, "total_prec", "sd")) %>%
  select(id, variable, signal_type, value) %>%
  pivot_wider(id_cols = c(id, signal_type), names_from = variable, values_from = value) %>%
  na.omit()

fig_s4 <- ggplot(fig_s4_data, aes(x = total_prec, y = e_total)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~ signal_type) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  theme_bw()

ggsave("output/figures_tables/figS4.pdf", fig_s4,
       width = 150, height = 75, units = "mm")

# correlation coefficients
R_mean <- cor.test(fig_s4_data %>% filter(signal_type == "mean") %>% pull(total_prec),
                   fig_s4_data %>% filter(signal_type == "mean") %>% pull(e_total))

R_sd <- cor.test(fig_s4_data %>% filter(signal_type == "sd") %>% pull(total_prec),
                 fig_s4_data %>% filter(signal_type == "sd") %>% pull(e_total))

# regression lines
fig_s4_mean_model <- lm(e_total ~ total_prec, fig_s4_data %>% filter(signal_type == "mean")) %>%
  summary()

fig_s4_sd_model <- lm(e_total ~ total_prec, fig_s4_data %>% filter(signal_type == "sd")) %>%
  summary()

dams <- readRDS("output/temp_rds_gpkg/07_dor.rds") %>%
  group_by(id, variable) %>%
  summarise(value = max(value) - min(value)) %>%
  mutate(significant = NA,
         signal_type = "abs_change") %>%
  ungroup()

# numbers for Section 3.2
damming_increases <- dams %>%
  left_join(cases_nogeom %>% select(id, area_km2), by = "id") %>%
  mutate(incr = ifelse(value == 0, FALSE, TRUE),
         islarge = ifelse(area_km2 > 10000, TRUE, FALSE)) %>%
  group_by(islarge, incr) %>%
  summarise(nn = n()) %>%
  mutate(freq = nn / sum(nn)) %>%
  ungroup()

# PLOT PRECIPITATION TRENDS ----------------------------------------------------

prec_trend <- driver_trends %>%
  filter(variable == "total_prec" & signal_type %in% c("mean", "sd")) %>%
  mutate(prec_trend = case_when(
    value < 0 ~ "decr",
    value > 0 ~ "incr"
  )) %>%
  rename(prec_significant = significant) %>%
  select(-variable) %>%
  pivot_wider(id_cols = id, names_from = signal_type,
              values_from = c(prec_trend, prec_significant, value))

prec_trend_magnitudes <- prec_trend %>%
  left_join(ctm_frc, by = "id")

prec_outlier_limits <- prec_trend_magnitudes %>%
  group_by(case, area_class) %>%
  summarise(minval_mean = mean(value_mean) - 2 * sd(value_mean),
            maxval_mean = mean(value_mean) + 2 * sd(value_mean),
            minval_sd = mean(value_sd) - 2 * sd(value_sd),
            maxval_sd = mean(value_sd) + 2 * sd(value_sd))

plt_data_mean <- .prep_prec_trend_plt(prec_trend_magnitudes, prec_outlier_limits, "mean")
plt_data_sd <- .prep_prec_trend_plt(prec_trend_magnitudes, prec_outlier_limits, "sd")

plt_data_iqrs <- bind_rows(plt_data_mean$iqr, plt_data_sd$iqr) %>%
  mutate(prec_trend = paste0(prec_trend, "_", signal_type))

plt_data_group_means <- bind_rows(plt_data_mean$group_mean, plt_data_sd$group_mean) %>%
  mutate(prec_trend = paste0(prec_trend, "_", signal_type))

write_csv(plt_data_iqrs %>% select(case, prec_trend, nn),
          "output/figures_tables/fig5_sample_sizes.csv")

plt_prec_magnitudes <- ggplot() +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_crossbar(data = plt_data_iqrs,
                aes(x = case, y = q75_value,
                    ymin = q25_value, ymax = q75_value, fill = prec_trend),
                position = "dodge", linetype = 0) +
  geom_text(data = plt_data_iqrs,
            aes(x = case, y = 0, label = nn, color = prec_trend)) +
  facet_grid(rows = vars(signal_type),
             space = "free", switch = "both") +
  scale_fill_manual(name = "precipitation trend",
                    values = c("grey65", "grey65", "#fdb462", "#bebada",
                               "#80b1d3", "#8dd3c7")) +
  geom_crossbar(data = plt_data_group_means,
                aes(x = case, y = mean_value, ymin = mean_value, ymax = mean_value),
                color = "grey40") +
  scale_color_manual(name = "precipitation trend",
                     values = c("grey65", "grey65", "#fdb462", "#bebada",
                                        "#80b1d3", "#8dd3c7")) +
  scale_x_discrete(name = "streamflow regime changeclass",
                   position = "top",
                   labels = ~ str_replace_all(.x, "aa_|ab_|ac_|ad_|ae_", "")) +
  ylab("Theil-Sen slope (mm/month)") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50", linewidth = 0.25),
        panel.grid.minor.y = element_line(color = "grey50", linewidth = 0.25),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.x.top = element_text())

ggsave("output/figures_tables/fig5.pdf", plt_prec_magnitudes,
       width = 150, height = 100, units = "mm")

# ADD TRENDS IN ET, WATERUSE AND DOR -------------------------------------------

all_trends <- driver_trends %>%
  filter((variable == "e_total" & signal_type %in% c("mean", "sd")) |
         (variable == "total_prec" & signal_type %in% c("mean", "sd")) |
         (variable == "wateruse_all_sectors" & signal_type == "mean")) %>%
  bind_rows(dams) %>%
  left_join(ctm_frc, by = "id")

# output driver trends for mapping
all_trends %>%
  mutate(trend_component = paste0(variable, "_", signal_type)) %>%
  select(id, case, area_class, trend_component, value) %>%
  pivot_wider(id_cols = c(id, case, area_class),
              names_from = trend_component,
              values_from = value) %>%
  left_join(gg %>% select(id, render_order), by = "id") %>%
  st_as_sf() %>%
  st_write("output/figures_tables/fig4_driver_trend_components.gpkg", delete_dsn = TRUE)

frc_classes <- sort(unique(ctm_frc$case))
driver_variables <- unique(all_trends$variable)
coll <- tibble()

for (i in 1:length(frc_classes)) {

  ctm_in_frc <- cases_nogeom %>%
    filter(case == frc_classes[i])

  trend_partition <- all_trends %>%
    filter(id %in% ctm_in_frc$id) %>%
    left_join(prec_trend %>% select(-starts_with("value")), by = "id") %>%
    select(-starts_with("prec_significant"))

  outlier_limits <- trend_partition %>%
    group_by(variable, area_class, signal_type) %>%
    summarise(minlim = mean(value) - 2 * sd(value),
              maxlim = mean(value) + 2 * sd(value)) %>%
    ungroup()

  trend_partition_olrm <- trend_partition %>%
    left_join(outlier_limits, by = c("variable", "area_class", "signal_type")) %>%
    filter(value > minlim & value < maxlim)

  trend_all <- trend_partition_olrm %>%
    group_by(variable, area_class, signal_type) %>%
    summarise(nn = n(),
              mean_value = mean(value),
              median_value = median(value),
              q25_value = quantile(value, 0.25),
              q75_value = quantile(value, 0.75)) %>%
    mutate(frc_class = frc_classes[i],
           prec_trend = "all") %>%
    select(frc_class, variable, everything()) %>%
    ungroup()

  trend_prec_mean <- trend_partition_olrm %>%
    select(-prec_trend_sd) %>%
    rename(prec_trend = prec_trend_mean) %>%
    filter(!(variable == "e_total" & signal_type == "sd")) %>%
    group_by(variable, prec_trend, area_class) %>%
    summarise(nn = n(),
              mean_value = mean(value),
              median_value = median(value),
              q25_value = quantile(value, 0.25),
              q75_value = quantile(value, 0.75)) %>%
    mutate(frc_class = frc_classes[i],
           prec_trend = paste0("mean_prec", prec_trend)) %>%
    select(frc_class, variable, everything()) %>%
    ungroup()

  trend_prec_sd <- trend_partition_olrm %>%
    select(-prec_trend_mean) %>%
    rename(prec_trend = prec_trend_sd) %>%
    filter(!(variable == "e_total" & signal_type == "mean")) %>%
    group_by(variable, prec_trend, area_class) %>%
    summarise(nn = n(),
              mean_value = mean(value),
              median_value = median(value),
              q25_value = quantile(value, 0.25),
              q75_value = quantile(value, 0.75)) %>%
    mutate(frc_class = frc_classes[i],
           prec_trend = paste0("sd_prec", prec_trend)) %>%
    select(frc_class, variable, everything()) %>%
    ungroup()

  coll <- coll %>%
    bind_rows(trend_all,
              trend_prec_mean,
              trend_prec_sd)
}

# PLOT ALL TRENDS IN ONE -------------------------------------------------------

plt_coll <- coll %>%
  arrange(frc_class, variable, prec_trend) %>%
  mutate(signal_type = ifelse(is.na(signal_type), "dummy", signal_type)) %>%
  select(frc_class, variable, area_class, prec_trend, everything())

p1 <- plt_coll %>%
  filter(frc_class %in% c("aa_shift_down", "ab_shift_up", "ae_other") & grepl("all|mean_", prec_trend)) %>%
  filter(!(variable == "e_total" & signal_type == "sd")) %>%
  filter(!(variable == "total_prec" & signal_type == "sd"))

p2 <- plt_coll %>%
  filter(frc_class %in% c("ac_shrink", "ad_expand") & grepl("all|sd_", prec_trend)) %>%
  filter(!(variable == "e_total" & signal_type == "mean")) %>%
  filter(!(variable == "total_prec" & signal_type == "mean"))

plt_coll <- bind_rows(p1, p2)
write_csv(plt_coll %>% select(-signal_type), "output/figures_tables/fig6_sample_sizes_long.csv")

plt_coll %>%
  select(frc_class, variable, area_class, prec_trend, nn) %>%
  pivot_wider(id_cols = c(frc_class, area_class, prec_trend),
              names_from = variable,
              values_from = nn) %>%
  arrange(area_class, frc_class, prec_trend) %>%
  select(area_class, frc_class, prec_trend, e_total, wateruse_all_sectors, degree_of_regulation) %>%
  write_csv("output/figures_tables/fig6_tableS2_sample_sizes_wide.csv")

plt_coll_frc_only <- plt_coll %>%
  filter(prec_trend == "all") %>%
  mutate(variable = case_when(
    variable == "total_prec" ~ "a_precipitation",
    variable == "e_total" ~ "b_evapotranspiration",
    variable == "wateruse_all_sectors" ~ "c_wateruse",
    variable == "degree_of_regulation" ~ "d_degree_of_regulation"
  ))

scale_fixer <- plt_coll_frc_only %>%
  select(frc_class, variable, area_class) %>%
  mutate(ymin_dummy = case_when(
    variable == "a_precipitation" ~ -0.35,
    variable == "b_evapotranspiration" ~ -0.12,
    variable == "c_wateruse" ~ -0.003,
    variable == "d_degree_of_regulation" ~ 0
  )) %>%
  mutate(ymax_dummy = case_when(
    variable == "a_precipitation" ~ 0.42,
    variable == "b_evapotranspiration" ~ 0.15,
    variable == "c_wateruse" ~ 0.0055,
    variable == "d_degree_of_regulation" ~ 55
  )) %>%
  pivot_longer(cols = c(ymin_dummy, ymax_dummy),
               names_to = "name",
               values_to = "scale_fixer") %>%
  select(-name)

plt_all_trends <- ggplot(plt_coll_frc_only) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_crossbar(aes(x = frc_class, y = mean_value,
                    ymin = q25_value, ymax = q75_value, fill = frc_class),
                position = "dodge", linetype = 0) +
  geom_crossbar(aes(x = frc_class, y = mean_value,
                    ymin = mean_value, ymax = mean_value),
                color = "grey60",
                position = "dodge") +
  geom_point(data = scale_fixer,
             aes(x = frc_class, y = scale_fixer), color = "red") +
  facet_grid(rows = vars(variable),
             cols = vars(area_class),
             scales = "free", space = "free_x", switch = "both") +
  scale_x_discrete(name = "streamflow regime changeclass",
                   position = "top",
                   labels = ~ str_replace_all(.x, "aa_|ab_|ac_|ad_|ae_|b_nonconf_|mix_", "")) +
  scale_y_continuous(position = "right") +
  scale_fill_manual(name = "FRC class",
                    values = c("#f4c40f", "#1f6e9d", "#aa7aa1", "#92c051", "#dfd7b2"),
                    labels = c("shift down", "shift up", "shrink", "expand", "other")) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50", linewidth = 0.25),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.x.top = element_text())

ggsave("output/figures_tables/fig6.pdf", plt_all_trends)


plt_coll_separate_prec <- plt_coll %>%
  mutate(variable = case_when(
    variable == "total_prec" ~ "a_precipitation",
    variable == "e_total" ~ "b_evapotranspiration",
    variable == "wateruse_all_sectors" ~ "c_wateruse",
    variable == "degree_of_regulation" ~ "d_degree_of_regulation"
  ))


plt_all_trends_separate_prec <- ggplot(plt_coll_separate_prec) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_crossbar(aes(x = frc_class, y = mean_value,
                    ymin = q25_value, ymax = q75_value, fill = prec_trend),
                position = "dodge", linetype = 0) +
  geom_crossbar(aes(x = frc_class, y = mean_value,
                    ymin = mean_value, ymax = mean_value, color = prec_trend),
                position = "dodge") +
  facet_grid(rows = vars(variable),
             cols = vars(area_class),
             scales = "free", space = "free_x", switch = "both") +
  scale_x_discrete(name = "streamflow regime changeclass",
                   position = "top",
                   labels = ~ str_replace_all(.x, "aa_|ab_|ac_|ad_|ae_", "")) +
  scale_y_continuous(position = "right") +
  scale_fill_manual(name = "grouping by precipitation trend",
                    values = c("grey65", "#fdb462", "#80b1d3", "#bebada", "#8dd3c7"),
                    labels = c("all", "mean P decreases", "mean P increases",
                               "sd of P decreases", "sd of P increases")) +
  scale_color_manual(name = NULL,
                     values = rep("grey40", 5)) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50", linewidth = 0.25),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.x.top = element_text())

ggsave("output/figures_tables/figS5.pdf", plt_all_trends_separate_prec)

